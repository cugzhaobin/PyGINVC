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
    
    def buildG(self, lon, lat):
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
    
    def assign_point(self, points, datatype="gps"):
        BlockID = np.full(len(points), -1, dtype=int)
        for i in range(len(points)):
            block = self.find_block(points[i])
            
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
