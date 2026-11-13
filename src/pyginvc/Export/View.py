#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 4 2016

import h5py, os, logging
import numpy as np
import matplotlib.pyplot as plt

from pyginvc.libs import geotools as gt
from pyginvc.Geometry.Patch import Fault
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def _data_extent(lon, lat, margin=0.2):
    '''
    Return [lon_min, lon_max, lat_min, lat_max] with a margin for map extent.
    '''
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    return [lon.min()-margin, lon.max()+margin,
            lat.min()-margin, lat.max()+margin]

def _auto_quiver_scale(lon, lat, *uv, fraction=0.1, margin=0.2):
    '''
    Return a quiver ``scale`` (used together with ``scale_units='xy'``) so
    that the longest arrow spans about ``fraction`` of the shorter map side.

    Parameters
    ----------
    lon, lat : array_like
        Station coordinates in degrees.
    *uv : tuple of (u, v) array_like
        One or more vector fields; all of them share the returned scale so
        observation, prediction and residual are directly comparable.
    fraction : float, optional
        Target length of the longest arrow as a fraction of the map side.
    margin : float, optional
        Same extent margin as ``_data_extent`` (degrees).

    Returns
    -------
    float or None
        Quiver scale, or None when it cannot be determined.
    '''
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    if lon.size == 0:
        return None
    span = min(np.ptp(lon), np.ptp(lat)) + 2.0 * margin
    mags = [np.hypot(np.asarray(u, dtype=float), np.asarray(v, dtype=float)).max()
            for u, v in uv]
    mags = [m for m in mags if np.isfinite(m)]
    maxmag = max(mags) if mags else 0.0
    if span <= 0 or maxmag <= 0:
        return None
    return maxmag / (fraction * span)

def _map_axes(fig, position, extent, ccrs, cfeature):
    '''
    Add a PlateCarree GeoAxes with coastlines, borders and labelled gridlines.
    '''
    ax = fig.add_subplot(position, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':')
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray',
                      alpha=0.4, linestyle='--')
    gl.top_labels   = False
    gl.right_labels = False
    return ax

def _sar_counts(h5):
    '''
    Number of points of each InSAR file (track).

    Falls back to one group for solution files written before ``obs/n_sar``
    existed, or when the stored counts do not match the stored points.
    '''
    n_total = len(h5['obs']['llh_sar'][()])
    if n_total == 0:
        return np.array([], dtype=int)
    counts = np.array([n_total])
    if 'n_sar' in h5['obs'].keys():
        stored = np.atleast_1d(h5['obs']['n_sar'][()]).astype(int)
        if len(stored) > 0 and stored.sum() == n_total:
            counts = stored
    return counts

def plot_obs_mod(solutionfile, scale=None):
    '''
    Plot observation and prediction for GPS and InSAR data

    Parameters
    ----------
    solutionfile : string
        HDF5 file output from PyGINV.Inverse
    scale : float, optional
        Scale of the GPS vectors (``scale_units='xy'``). None means the
        arrows are scaled automatically so that the longest one spans about
        10% of the map side; observation, prediction and residual share the
        same scale for direct comparison.

    Returns
    -------
    None.

    '''

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        logging.warning('cartopy is required for plotting observation and '
                        'prediction, skipping. Install it with '
                        '"pip install cartopy".')
        return None

    with h5py.File(solutionfile, 'r') as h5:
        #
        # For GPS stations
        #
        if 'obs' not in h5.keys():
            print('No group named obs.')
            return None
        else:
            if 'llh_gps' not in h5['obs'].keys():
                print('No dataset named llh_gps.')
                return None
        llh_gps = h5['obs']['llh_gps'][()]
        origin  = h5['flt']['origin'][()]
        enu     = np.zeros((len(llh_gps),2)) 
        for i in range(len(llh_gps)):
            enu[i] = gt.llh2localxy(llh_gps[i], origin)
        if len(llh_gps) == 0:
            print('No GPS data.')
            ndim  = 0
        else:
            ndim  = h5['obs']['ndim'][()]
            d_gps = h5['obs']['d_gps'][()]
            d_gps = d_gps.reshape(len(llh_gps), ndim)
            m_gps = h5['sol']['dhat'][0:ndim*len(llh_gps)]
            m_gps = m_gps.reshape(len(llh_gps), ndim)
            trans = ccrs.PlateCarree()
            ext   = _data_extent(llh_gps[:,1], llh_gps[:,0])
            # 箭头长度自适应：观测/预测/残差共用同一比例，便于直接比较
            qscale = scale
            if qscale is None:
                qscale = _auto_quiver_scale(
                    llh_gps[:,1], llh_gps[:,0],
                    (d_gps[:,0], d_gps[:,1]), (m_gps[:,0], m_gps[:,1]))
            fig   = plt.figure(figsize=(12,5), layout='constrained')
            ax1   = _map_axes(fig, 121, ext, ccrs, cfeature)
            ax2   = _map_axes(fig, 122, ext, ccrs, cfeature)

            ax1.quiver(llh_gps[:,1], llh_gps[:,0], d_gps[:,0], d_gps[:,1],
                       transform=trans, color='b', scale=qscale,
                       scale_units='xy')
            ax1.quiver(llh_gps[:,1], llh_gps[:,0], m_gps[:,0], m_gps[:,1],
                       transform=trans, color='r', scale=qscale,
                       scale_units='xy')
            if ndim == 3:
                ax1.scatter(llh_gps[:,1], llh_gps[:,0], transform=trans,
                            c=d_gps[:,2], cmap=plt.cm.jet, s=50,
                            linewidths=0.1, edgecolors='black')
                ax1.scatter(llh_gps[:,1], llh_gps[:,0], transform=trans,
                            c=m_gps[:,2], cmap=plt.cm.jet, s=15,
                            linewidths=0.1, edgecolors='black')
            ax1.set_title('Observation (blue) VS Prediction (red)')

            ax2.quiver(llh_gps[:,1], llh_gps[:,0],
                       d_gps[:,0]-m_gps[:,0], d_gps[:,1]-m_gps[:,1],
                       transform=trans, color='b', scale=qscale,
                       scale_units='xy')
            if ndim == 3:
                ax2.scatter(llh_gps[:,1], llh_gps[:,0], transform=trans,
                            c=d_gps[:,2]-m_gps[:,2], cmap=plt.cm.jet,
                            linewidths=0.1, edgecolors='black')
            ax2.set_title('Residual')

            fig.savefig('gps_mod_obs.pdf', format='pdf', bbox_inches='tight')
            plt.show()

        # 
        # For InSAR
        #
        if 'llh_sar' in h5['obs'].keys():
            llh_sar = h5['obs']['llh_sar'][()]
            if len(llh_sar) == 0:
                print(' No InSAR points')
            else:
                d_sar = h5['obs']['d_sar'][()]
                # dhat 的行序为 GPS, 水准, InSAR，水准数据夹在中间
                n_lev = len(h5['obs']['llh_lev'][()]) if 'llh_lev' in h5['obs'] else 0
                offset = ndim*len(llh_gps) + n_lev
                m_sar = h5['sol']['dhat'][offset:offset+len(llh_sar)]

                # 多个 InSAR 文件时逐个文件单独绘图
                counts = _sar_counts(h5)
                multi  = len(counts) > 1

                # 逐 track 绘图：每个文件一组观测/预测/残差
                start = 0
                for itrack, count in enumerate(counts):
                    end = start + count
                    if count == 0:
                        start = end
                        continue

                    ob   = d_sar[start:end]
                    mb   = m_sar[start:end]
                    rb   = ob - mb
                    lon  = llh_sar[start:end, 1]
                    lat  = llh_sar[start:end, 0]
                    tag  = ' Track %d' % (itrack + 1) if multi else ''
                    trans = ccrs.PlateCarree()
                    ext   = _data_extent(lon, lat)
                    # 观测和预测共用色标范围，便于直接比较
                    vmin  = min(ob.min(), mb.min())
                    vmax  = max(ob.max(), mb.max())
                    # 残差用以 0 为中心的独立对称色标，避免被观测幅值掩盖
                    rmax  = np.abs(rb).max()
                    if rmax == 0:
                        rmax = 1.0
                    fig   = plt.figure(figsize=(18,5), layout='constrained')
                    ax1   = _map_axes(fig, 131, ext, ccrs, cfeature)
                    ax2   = _map_axes(fig, 132, ext, ccrs, cfeature)
                    ax3   = _map_axes(fig, 133, ext, ccrs, cfeature)

                    s = ax1.scatter(lon, lat, c=ob,
                                    cmap=plt.cm.jet, s=2, vmin=vmin, vmax=vmax,
                                    transform=trans)
                    ax1.set_title('InSAR%s LOS Observation' % tag)
                    cb = fig.colorbar(s, ax=ax1, shrink=0.8)
                    cb.set_label('Line of sight changes (mm)')

                    s = ax2.scatter(lon, lat, c=mb,
                                    cmap=plt.cm.jet, s=2, vmin=vmin, vmax=vmax,
                                    transform=trans)
                    ax2.set_title('InSAR%s LOS Prediction' % tag)
                    cb = fig.colorbar(s, ax=ax2, shrink=0.8)
                    cb.set_label('Line of sight changes (mm)')

                    s = ax3.scatter(lon, lat, c=rb,
                                    cmap=plt.cm.jet, s=2, vmin=-rmax, vmax=rmax,
                                    transform=trans)
                    ax3.set_title('InSAR%s LOS Residual' % tag)
                    cb = fig.colorbar(s, ax=ax3, shrink=0.8)
                    cb.set_label('Line of sight changes (mm)')

                    if multi:
                        fname = 'sar_mod_obs_%d.pdf' % (itrack + 1)
                    else:
                        fname = 'sar_mod_obs.pdf'
                    fig.savefig(fname, format='pdf', bbox_inches='tight')
                    plt.show()
                    start = end


def plot_sar_ramp(solutionfile, show=True):
    '''
    Plot InSAR original data, the estimated ramp and the ramp-removed data.

    One figure per InSAR file (track) that estimated a ramp; tracks without
    ramp are skipped, as is the whole call when no InSAR ramp was estimated.

    Parameters
    ----------
    solutionfile : string
        HDF5 file output from PyGINV program
    show : bool, optional
        Call ``plt.show()`` when True. The default is True.

    Returns
    -------
    None.

    '''
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        logging.warning('cartopy is required for plotting InSAR ramp, '
                        'skipping. Install it with "pip install cartopy".')
        return None

    with h5py.File(solutionfile, 'r') as h5:
        if 'd_sar_ramp' not in h5['obs'].keys():
            print('No InSAR ramp in %s, skip ramp plot.' % solutionfile)
            return None

        llh_sar = h5['obs']['llh_sar'][()]
        d_sar   = h5['obs']['d_sar'][()]
        d_ramp  = h5['obs']['d_sar_ramp'][()]
        if len(d_sar) == 0 or len(d_ramp) != len(d_sar):
            print('No InSAR ramp data to plot, skip.')
            return None

        counts = _sar_counts(h5)
        if len(counts) == 0:
            print('No InSAR points, skip ramp plot.')
            return None

        # 逐 track 是否估计了 ramp；旧文件未记录时按全部估计处理
        flags = np.ones(len(counts), dtype=bool)
        if 'sar_ramp' in h5['obs'].keys():
            stored = np.atleast_1d(h5['obs']['sar_ramp'][()]).astype(bool)
            if len(stored) == len(counts):
                flags = stored
        if not flags.any():
            print('No InSAR ramp was estimated, skip ramp plot.')
            return None

        multi = len(counts) > 1
        start = 0
        for itrack, count in enumerate(counts):
            end = start + count
            if count == 0 or not flags[itrack]:
                start = end
                continue

            ob   = d_sar[start:end]
            rp   = d_ramp[start:end]
            cor  = ob - rp
            lon  = llh_sar[start:end, 1]
            lat  = llh_sar[start:end, 0]
            tag  = ' Track %d' % (itrack + 1) if multi else ''
            trans = ccrs.PlateCarree()
            ext   = _data_extent(lon, lat)
            # 三幅图共用一个色标，原始/ramp/去ramp 可直接比较
            vmin = min(ob.min(), rp.min(), cor.min())
            vmax = max(ob.max(), rp.max(), cor.max())
            if vmax <= vmin:
                vmin, vmax = vmin - 1.0, vmax + 1.0

            fig = plt.figure(figsize=(18, 5), layout='constrained')
            for k, (title, value) in enumerate(
                    (('Original', ob), ('Ramp', rp),
                     ('Ramp removed', cor)), start=1):
                ax = _map_axes(fig, 130 + k, ext, ccrs, cfeature)
                s  = ax.scatter(lon, lat, c=value, cmap=plt.cm.jet, s=2,
                                vmin=vmin, vmax=vmax, transform=trans)
                ax.set_title('InSAR%s LOS %s' % (tag, title))
                cb = fig.colorbar(s, ax=ax, shrink=0.8)
                cb.set_label('Line of sight changes (mm)')

            fname = 'sar_ramp_%d.pdf' % (itrack + 1) if multi else 'sar_ramp.pdf'
            fig.savefig(fname, format='pdf', bbox_inches='tight')
            if show:
                plt.show()
            else:
                plt.close(fig)
            start = end
    return None


def _set_box_aspect(ax, x, y, z, coordtype='llh', depth_exaggeration=2.0):
    '''
    Set the 3D box aspect from true horizontal distances and exaggerated depth.

    matplotlib's 3D axes default to a cube regardless of data ranges, which
    makes longitude/latitude and depth look mismatched.  Here the east and
    north edges are scaled by their true ground distances (km) and the depth
    edge is multiplied by ``depth_exaggeration``.

    Parameters
    ----------
    ax : mpl_toolkits.mplot3d.Axes3D
    x, y, z : array_like
        Coordinates of the plotted vertices.
    coordtype : string, optional
        'llh' (x=lon deg, y=lat deg) or 'enu' (x/y in km).
    depth_exaggeration : float, optional
        Vertical exaggeration of the depth axis. The default is 2.0.
    '''
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    z = np.asarray(z, dtype=float).ravel()
    dx, dy, dz = np.ptp(x), np.ptp(y), np.ptp(z)

    if coordtype == 'llh':
        # 1 deg latitude ~ 111.32 km; 1 deg longitude shrinks with cos(lat)
        lat0 = np.deg2rad(np.nanmean(y))
        dx = dx * 111.32 * max(np.cos(lat0), 1e-3)
        dy = dy * 111.32

    # 退化方向（如恰好为一条线/面）给一个极小比例，避免 3D 轴报错
    eps = max(dx, dy, dz, 1.0) * 1e-3
    ax.set_box_aspect((max(dx, eps), max(dy, eps),
                       max(dz * depth_exaggeration, eps)))


def _slip_arrow_segments(felem, coordtype, slip_scale=0.6):
    '''
    Build in-plane slip-arrow segments for rectangular fault patches.

    Every arrow starts at the patch centre, its length is proportional to
    sqrt(strike_slip**2 + dip_slip**2) and the opening component is ignored.
    The shaft is assembled from the patch edges, so it lies in the fault
    plane, and the two barbs are built from the same in-plane basis.  This
    is the same construction as ``Fault.plot_faultgeom``, and it avoids the
    out-of-plane barbs produced by ``Axes3D.quiver`` on steeply dipping
    faults.

    Parameters
    ----------
    felem : pandas.DataFrame
        Output of ``Fault.FaultGeom2AllVertex`` with the solution slip
        stored in its strike_slip/dip_slip columns.
    coordtype : {'llh', 'enu'}
        Coordinate system of the 3D axes.
    slip_scale : float, optional
        Longest arrow as a fraction of the median patch length.

    Returns
    -------
    list of (2, 3) arrays
        Line segments [p0, p1] forming the shafts and arrow heads.
    '''
    if coordtype == 'llh':
        comps = ('lon', 'lat', 'dep')
    elif coordtype == 'enu':
        comps = ('e', 'n', 'u')
    else:
        return []

    def corner(tag):
        return np.column_stack(
            [np.asarray(felem['%s_%s' % (tag, c)], dtype=float)
             for c in comps])

    lt, rt = corner('lt'), corner('rt')
    lb, rb = corner('lb'), corner('rb')

    ss  = np.nan_to_num(np.asarray(felem['strike_slip'], dtype=float))
    ds  = np.nan_to_num(np.asarray(felem['dip_slip'], dtype=float))
    mag = np.hypot(ss, ds)
    if mag.size == 0 or mag.max() <= 0:
        return []

    # Physical edge lengths (km) set the slip -> plot-length scale.
    e_s = np.column_stack([np.asarray(felem['rt_e'] - felem['lt_e'], float),
                           np.asarray(felem['rt_n'] - felem['lt_n'], float),
                           np.asarray(felem['rt_u'] - felem['lt_u'], float)])
    e_d = np.column_stack([np.asarray(felem['lb_e'] - felem['lt_e'], float),
                           np.asarray(felem['lb_n'] - felem['lt_n'], float),
                           np.asarray(felem['lb_u'] - felem['lt_u'], float)])
    L_s = np.linalg.norm(e_s, axis=1)
    L_d = np.linalg.norm(e_d, axis=1)
    good = np.isfinite(L_s) & np.isfinite(L_d) & (L_s > 0) & (L_d > 0)
    if not np.any(good):
        return []
    arrow_scale = slip_scale * np.median(L_s[good]) / mag.max()

    p_s    = rt - lt                  # right-hand strike edge
    p_d    = lb - lt                  # up-dip edge
    centre = 0.25 * (lt + rt + lb + rb)

    # Unit edge vectors in physical length.  In lon/lat/depth coordinates the
    # horizontal components are degrees, so a geometric normal mixes units
    # and would make the barbs far too long.  Building the in-plane basis
    # from the patch edges keeps shaft and barbs physically consistent.
    L_s = np.where(good, L_s, 1.0)
    L_d = np.where(good, L_d, 1.0)
    u_s = p_s / L_s[:, None]
    u_d = p_d / L_d[:, None]

    head_frac = 0.28                  # barb length as a fraction of shaft
    head_ang  = np.deg2rad(22.0)
    floor     = 1e-9 * mag.max()      # drop round-off-sized slip
    segments  = []
    for i in range(len(felem)):
        if mag[i] <= floor or not good[i]:
            continue
        theta = np.arctan2(ds[i], ss[i])
        v_hat = np.cos(theta)*u_s[i] + np.sin(theta)*u_d[i]
        w_hat = -np.sin(theta)*u_s[i] + np.cos(theta)*u_d[i]
        vlen  = mag[i] * arrow_scale
        if vlen <= 0:
            continue

        tip  = centre[i] + vlen * v_hat
        head = head_frac * vlen
        back = tip - head * np.cos(head_ang) * v_hat
        side = head * np.sin(head_ang)
        segments.append([centre[i], tip])
        segments.append([tip, back + side * w_hat])
        segments.append([tip, back - side * w_hat])
    return segments


def plot_slip_3d(solutionfile, elevation=40, azimuth=-81, coordtype='llh',
                 show=True, depth_exaggeration=2.0, show_slip=True,
                 slip_scale=0.6):
    '''
    Plot 3D finite slip distribution.

    Parameters
    ----------
    solutionfile : string
        HDF5 file output from PyGINV program
    elevation : float, optional
        look angle
    azimuth : float, optional
        look angle. The default is -81.
    coordtype : string, optional
        'llh/enu'. The default is 'llh'.
    depth_exaggeration : float, optional
        Vertical exaggeration of the depth axis. The default is 2.0.
    show_slip : bool, optional
        Draw the slip vector of every rectangular patch, as in the 2D
        plane-view figures: the tail is at the patch centre, the length is
        proportional to sqrt(strike_slip**2 + dip_slip**2) and the whole
        arrow (shaft and head) lies in the fault plane.  The default is True.
    slip_scale : float, optional
        Longest slip arrow as a fraction of the median patch length.
        The default is 0.6.

    Notes
    -----
    cartopy 的 GeoAxes 是二维投影，不能作为 mplot3d 的 3D 坐标轴使用；
    这里通过按真实地面距离设置 box aspect 来保证经纬度与深度的比例一致。

    Returns
    -------
    None.

    '''

    with h5py.File(solutionfile, 'r') as h5:
        if 'flt' not in h5.keys():
            print('No fault information for plotting.')
            return None
        else:
            fig      = plt.figure()
            ax       = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=elevation, azim=azimuth)

            # h5py 读取字符串数据集返回 bytes，需要解码后再比较
            flttype  = h5['flt/flttype'][()]
            if isinstance(flttype, bytes):
                flttype = flttype.decode()

            if flttype == 'rectangle':
                geom     = h5['flt/geom'][()]
                dis_geom = h5['flt/dis_geom'][()]
                slip     = h5['sol/slip'][()]
                # 多个平滑因子时取最后一个（最终模型），与 GMT/模型输出一致
                if slip.ndim > 1:
                    slip = slip[-1]
                slip     = slip[0:len(geom)*3]
                slip     = slip.reshape(len(geom), 3)
                ts       = np.sqrt(slip[:,0]**2+slip[:,1]**2)
                tmax     = ts.max() if ts.size else 0.0
                ts1      = ts / tmax if tmax > 0 else np.zeros_like(ts)
                # Build a Fault view from the stored arrays.  Patch.Fault only
                # fills felem in load_fault(), which reads from a file, so set
                # the four attributes FaultGeom2AllVertex() depends on.
                flt      = Fault.__new__(Fault)
                flt.nf             = len(geom)
                flt.geom_grid      = np.asarray(geom)
                flt.dis_geom_grid  = np.array(dis_geom, dtype=float)
                # Store the solution slip in the geometry table so that
                # FaultGeom2AllVertex() returns strike_slip/dip_slip in the
                # same plotting convention as Fault.plot_faultgeom.
                flt.dis_geom_grid[:, 7:10] = slip.reshape(len(geom), 3)
                flt.origin         = np.array([np.mean(flt.geom_grid[:,[3,5]]),
                                               np.mean(flt.geom_grid[:,[4,6]])])
                felem_df = flt.FaultGeom2AllVertex()
                felem    = felem_df.values
    
                if coordtype == 'llh':
                    idx  = np.array([9,12,15,18])
                    x    = felem[:,idx]
                    y    = felem[:,idx+1]
                    z    = felem[:,idx+2]
                elif coordtype == 'enu':
                    idx  = np.array([21,24,27,30])
                    x    = felem[:,idx]
                    y    = felem[:,idx+1]
                    z    = felem[:,idx+2]

                # 经纬度按真实地面距离、深度按 depth_exaggeration 放大
                _set_box_aspect(ax, x, y, z, coordtype, depth_exaggeration)

                for i in range(len(felem)):
                    verts= [list(zip(x[i], y[i], z[i]))]
                    subfault = Poly3DCollection(verts, facecolors='w', linewidths=1, alpha=0.3)
                    subfault.set_facecolor(plt.cm.jet(ts1[i]))
                    subfault.set_linewidth(0.05)
                    ax.add_collection3d(subfault)

                # Slip arrows, drawn as explicit in-plane segments.  The
                # arrowheads of Axes3D.quiver are rotated about a horizontal
                # axis and would leave a steeply dipping fault plane, so the
                # shafts and barbs are assembled from the patch edges instead.
                if show_slip:
                    segments = _slip_arrow_segments(felem_df, coordtype,
                                                    slip_scale)
                    if segments:
                        import matplotlib.patheffects as pe
                        from mpl_toolkits.mplot3d.art3d import Line3DCollection
                        arrows = Line3DCollection(segments, colors='w',
                                                  linewidths=0.9)
                        # Dark rim keeps the light arrows legible on both the
                        # dark blue and the yellow/red ends of the slip map.
                        arrows.set_path_effects([
                            pe.Stroke(linewidth=1.9, foreground='0.25'),
                            pe.Normal()])
                        ax.add_collection3d(arrows)
                s    = ax.scatter(geom[:,4], geom[:,3], c=ts/1000., cmap=plt.cm.jet, s=0.00001, lw=0)
                cb   = plt.colorbar(s, shrink=0.8, pad=-0.07)
                cb.set_label('Total slip(m)')
                if coordtype == 'llh':
                    ax.set_xlabel('Longitude')
                    ax.set_ylabel('Latitude')
                elif coordtype == 'enu':
                    ax.set_xlabel('East (km)')
                    ax.set_ylabel('North (km)')
                ax.set_zlabel('Depth (km)')
                plt.show()
            elif flttype == 'triangle':
                element  = h5['flt/element'][()]-1
                slip     = h5['sol/slip'][()]
                # 多个平滑因子时取最后一个（最终模型），与 GMT/模型输出一致
                if slip.ndim > 1:
                    slip = slip[-1]
                slip     = slip[0:len(element)*3]
                slip     = slip.reshape(len(element), 3)
                ts       = np.sqrt(slip[:,0]**2+slip[:,1]**2)
                tmax     = ts.max() if ts.size else 0.0
                ts1      = ts / tmax if tmax > 0 else np.zeros_like(ts)

                if coordtype == 'llh':
                    # lat, lon, down
                    vertex_llh = h5['flt/vertex_llh'][()]
                    _set_box_aspect(ax, vertex_llh[:,1], vertex_llh[:,0],
                                    -vertex_llh[:,2], 'llh',
                                    depth_exaggeration)
                    idx = element[:,0]
                    ax.scatter(vertex_llh[:,1],
                               vertex_llh[:,0],
                              -vertex_llh[:,2],
                               s=0.1)
                    for i in range(len(element)):
                        indx = element[i,:]
                        lat  = np.array(vertex_llh[indx, 0])
                        lon  = np.array(vertex_llh[indx, 1])
                        dep  = np.array(-vertex_llh[indx, 2])
                        verts= [list(zip(lon, lat, dep))]
                        subfault = Poly3DCollection(verts)
                        subfault.set_facecolor(plt.cm.jet(ts1[i]))
                        subfault.set_linewidth(1.5)
                        ax.add_collection3d(subfault)

                    idx = element[:,0]
                    s   = ax.scatter(vertex_llh[idx,1], 
                                     vertex_llh[idx,0], 
                                     c=ts, cmap=plt.cm.jet, s=1e-5, lw=0)
                    cb = plt.colorbar(s, shrink=0.8, pad=-0.07)
                    cb.set_label('Total slip (mm)')
        
                    plt.xlabel('Longitude')
                    plt.ylabel('Latitude')
                    ax.set_zlabel('Depth (km)')
            
                elif coordtype in {'neu', 'enu', 'NEU', 'ENU'}:
                    vertex_enu = h5['flt/vertex_enu'][()]
                    _set_box_aspect(ax, vertex_enu[:,1], vertex_enu[:,0],
                                    -vertex_enu[:,2], 'enu',
                                    depth_exaggeration)
                    idx = element[:,0]
                    s   = ax.scatter(vertex_enu[idx,1],
                                     vertex_enu[idx,0],
                                     -vertex_enu[idx,2],
                                     c=ts, cmap=plt.cm.jet, s=1e-5, lw=0)

                    for i in range(len(element)):
                        indx = element[i,:]
                        e    = np.array(vertex_enu[indx, 1])
                        n    = np.array(vertex_enu[indx, 0])
                        u    = np.array(-vertex_enu[indx, 2])
                        verts= [list(zip(e, n, u))]
                        subfault = Poly3DCollection(verts)
                        subfault.set_facecolor(plt.cm.jet(ts1[i]))
                        subfault.set_linewidth(1.5)
                        ax.add_collection3d(subfault)
                    # plot color bar
                    cb = plt.colorbar(s, shrink=0.8, pad=-0.07)
                    cb.set_label('Total slip (m)')
                    
                    
                    plt.xlabel('East (m)')
                    plt.ylabel('North (m)')
                    ax.set_zlabel('Depth (m)')
                plt.show()


def plot_lcurve(solutionfile, curvefile=None, show=True):
    '''
    Plot the trade-off (L-) curve between model roughness and data misfit.

    The curve is only drawn when at least two smoothing factors were tested,
    i.e. when the ``ruff_misfit`` table exists.

    Parameters
    ----------
    solutionfile : string
        HDF5 solution file; used to locate ``ruff_misfit`` in the same folder.
    curvefile : string, optional
        Trade-off table written by ``OutputUtils.write_lcurve``. None means
        ``ruff_misfit`` next to ``solutionfile``.
    show : bool, optional
        Call ``plt.show()`` when True. The default is True.

    Returns
    -------
    None.

    '''
    if curvefile is None:
        curvefile = os.path.join(
            os.path.dirname(os.path.abspath(solutionfile)), 'ruff_misfit')
    if not os.path.isfile(curvefile):
        print('No trade-off table (%s), skip L-curve.' % curvefile)
        return None

    table = np.atleast_2d(np.loadtxt(curvefile, comments='#'))
    if table.size == 0 or table.shape[0] < 2:
        print('Need at least two smoothing factors for an L-curve, skip.')
        return None

    ruff   = table[:, 0]
    misfit = table[:, 1]
    smo    = table[:, 2]

    fig, ax = plt.subplots(figsize=(6, 5), layout='constrained')
    # 粗糙度和残差都为正时用双对数坐标，L 形最清楚
    loglog = np.all(ruff > 0) and np.all(misfit > 0)
    if loglog:
        ax.loglog(ruff, misfit, 'o-', color='b', lw=1, ms=5)
    else:
        ax.plot(ruff, misfit, 'o-', color='b', lw=1, ms=5)

    # 标注每个点对应的平滑因子
    for r, m, s in zip(ruff, misfit, smo):
        ax.annotate('%g' % s, (r, m), textcoords='offset points',
                    xytext=(4, 4), fontsize=8, color='0.3')

    # 取点到首尾连线垂距最大者为拐点，即折衷解
    if loglog and len(ruff) >= 3:
        lr      = np.log10(ruff)
        lm      = np.log10(misfit)
        order   = np.argsort(lr)
        scale_x = np.ptp(lr) or 1.0
        scale_y = np.ptp(lm) or 1.0
        x       = (lr[order] - lr.min()) / scale_x
        y       = (lm[order] - lm.min()) / scale_y
        dist    = np.abs((y[-1] - y[0]) * (x - x[0])
                         - (x[-1] - x[0]) * (y - y[0]))
        icor    = order[np.argmax(dist)]
        ax.plot(ruff[icor], misfit[icor], 'r*', ms=12,
                label='corner: smoothing = %g' % smo[icor])
        ax.legend(loc='best', fontsize=8)

    ax.set_xlabel('Model roughness')
    ax.set_ylabel('Data misfit')
    ax.set_title('Trade-off (L-) curve')
    ax.grid(True, ls='--', lw=0.3, alpha=0.5)

    outpdf = os.path.join(os.path.dirname(os.path.abspath(curvefile)),
                          'lcurve.pdf')
    fig.savefig(outpdf, format='pdf', bbox_inches='tight')
    if show:
        plt.show()
    return None