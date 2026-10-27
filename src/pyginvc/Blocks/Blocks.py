import geopandas as gpd
import numpy as np
class Block:
    def __init__(self, block_id, geometry, name=None):
        self.id       = block_id
        self.name     = name
        self.geometry = geometry
        
        self.wx = 0
        self.wy = 0
        self.wz = 0
        
        self.gps = []
        self.sar = []
        
    def contains(self, point):
        return self.geometry.contains(point)
    
    @staticmethod
    def buildG(lon, lat):
        """
        Rotation matrix mapping an Euler vector to surface velocity.

        Rows are (east, north).  No Earth radius is applied, so omega comes out
        in whatever units make B @ omega match the observed velocities -- i.e.
        the same convention the inversion fitted.  Do NOT mix this with Euler
        poles given in rad/yr.
        """
        rlon = np.deg2rad(lon)
        rlat = np.deg2rad(lat)
        G    = np.array([[-np.sin(rlat)*np.cos(rlon), -np.sin(rlat)*np.sin(rlon), np.cos(rlat)],
               [np.sin(rlon), -np.cos(rlon), 0]])
        return G

class BlockModel:
    def __init__(self, geojsonfile):
        gpf = gpd.read_file(geojsonfile)
        self.blocks = []
        for i, row in gpf.iterrows():
            block = Block(i, row.geometry, None)
            self.blocks.append(block)
        self.nblock = len(self.blocks)
        self.gdf    = gpf
    
    def assign_point(self, points, datatype="gps", fallback="contains"):
        """
        fallback : 'contains' -- a point inside no polygon gets id -1, and
                                 GeoData.process_data then DROPS it silently.
                   'nearest'  -- such a point is given the id of the closest
                                 polygon instead.  Use it to salvage stations
                                 that fall in the gap between two independently
                                 digitised blocks, but note it changes which
                                 stations enter the inversion.
        """
        if fallback not in ("contains", "nearest"):
            raise ValueError("fallback must be 'contains' or 'nearest', got %r" % fallback)

        BlockID = np.full(len(points), -1, dtype=int)
        for i in range(len(points)):
            block = self.find_block(points[i])
            if block is None and fallback == "nearest":
                block = self.find_block_nearest(points[i])

            if block is not None:
                BlockID[i] = block.id
                if datatype == "gps":
                    block.gps.append(points[i])
                elif datatype == 'sar':
                    block.sar.append(points[i])
        return BlockID
        
    
    def find_block(self, point):
        from shapely.geometry import Point
        pt = Point(point[0], point[1])
        for i, row in self.gdf.iterrows():
            if row.geometry.contains(pt):
                return self.blocks[i]
        return None
    
    def buildG_at(self, lon, lat):
        """Rotation matrix at an arbitrary point, same convention as buildG_block."""
        return Block.buildG(lon, lat)

    def find_block_nearest(self, point):
        """
        Nearest block polygon.

        Fallback for points that fall outside every polygon -- either in a gap
        between two independently digitised blocks, or genuinely outside the
        tessellation.  Always returns something; callers that need to know
        whether the point was really inside should use find_block() first.
        """
        from shapely.geometry import Point
        pt = Point(point[0], point[1])
        best, bestd = None, np.inf
        for i, row in self.gdf.iterrows():
            d = row.geometry.distance(pt)
            if d < bestd:
                bestd, best = d, self.blocks[i]
        return best

    def flanking_blocks(self, lon, lat, az_deg, offsets_km=(5., 10., 20., 40.)):
        """
        The two blocks sitting on either side of a point on a fault trace.

        az_deg     : fault strike azimuth, degrees clockwise from north.
        offsets_km : trial distances along the fault normal.  The smallest offset
                     that puts both trial points inside *different* polygons wins,
                     so the sample stays as close to the trace as the digitised
                     geometry allows.

        Returns (id_left, id_right, exact, offset_km):
            id_left  : block on the side 90 deg COUNTER-clockwise from the strike
                       azimuth.  This is the +y side of Okada's fault-local frame,
                       which is what the sign of a strike-slip Green's function is
                       referenced to -- see Coupling.local_basis.
            id_right : the clockwise side.
            exact  : True if both trial points were inside a polygon.  False means
                     the nearest-polygon fallback was used.
            offset : the offset that succeeded, or nan.  A large value is not by
                     itself wrong -- if two blocks are the boundary pair, their
                     relative rotation is the long-term rate whether or not the
                     trace is drawn exactly on their shared edge -- but it does
                     mean the sample was taken far from the fault, so it is worth
                     reporting.  It DOES become wrong if the offset jumped across
                     a third block; check against a block map before trusting
                     anything with offset_km > 20.
        Returns (-1, -1, False, nan) if no usable pair could be found either way.
        """
        for off in offsets_km:
            a, b = self._sample_pair(lon, lat, az_deg, off, nearest=False)
            if a is not None and b is not None and a.id != b.id:
                return a.id, b.id, True, float(off)

        a, b = self._sample_pair(lon, lat, az_deg, offsets_km[0], nearest=True)
        if a is None or b is None or a.id == b.id:
            return -1, -1, False, np.nan
        return a.id, b.id, False, float(offsets_km[0])

    def _sample_pair(self, lon, lat, az_deg, offset_km, nearest=False):
        """
        Two points at +/- offset_km along the normal to a strike of az_deg.
        Returns (block_left, block_right); azimuths run clockwise, so the
        counter-clockwise (left) normal is at az_deg - 90.
        """
        n    = np.deg2rad(az_deg - 90.0)
        dlat = offset_km / 111.0
        dlon = offset_km / (111.0 * max(np.cos(np.deg2rad(lat)), 1e-6))
        de, dn = dlon*np.sin(n), dlat*np.cos(n)

        out = []
        for sgn in (1.0, -1.0):
            b = self.find_block((lon + sgn*de, lat + sgn*dn))
            if b is None and nearest:
                b = self.find_block_nearest((lon + sgn*de, lat + sgn*dn))
            out.append(b)
        return out[0], out[1]

    def buildG_block(self, ndim, block_type="gps"):
        G_list = []
        
        for block in self.blocks:
            if block_type == 'gps':
                pts = block.gps
            G_block = np.zeros((ndim*len(pts), 3))
            for i, pt in enumerate(pts):
                lon, lat = pt
                G_block[ndim*i:ndim*i+ndim,:] = block.buildG(lon, lat)
            G_list.append(G_block)
        return G_list
