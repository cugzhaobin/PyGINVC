"""
Plot Helpers - Visualization functions for the Streamlit GUI.

Each function returns a ``matplotlib.figure.Figure`` object instead of calling
``plt.show()``, so it can be rendered with ``st.pyplot(fig)`` in Streamlit.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for Streamlit
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from pyginvc.libs import geotools as gt


def plot_gps_obs_mod(data, sol, green=None, scale=200) -> Figure:
    """
    Plot GPS observation vs prediction as quiver maps.

    Parameters
    ----------
    data : GeoData
    sol : GeoInversion (with .dhat attribute)
    green : optional, unused (kept for API consistency)
    scale : float, quiver scale factor

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(12, 5))
    fig.set_tight_layout(True)

    if data.llh_gps.size == 0:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No GPS data", ha="center", va="center", fontsize=16)
        return fig

    llh_gps = data.llh_gps
    ndim = data.ndim
    d_gps = data.d_gps.reshape(len(llh_gps), ndim)

    # predicted GPS
    dhat = sol.dhat
    n_gps_comp = ndim * len(llh_gps)
    m_gps = dhat[:n_gps_comp].reshape(len(llh_gps), ndim)

    # Observation vs Prediction
    ax1 = fig.add_subplot(121)
    ax1.quiver(
        llh_gps[:, 1], llh_gps[:, 0],
        d_gps[:, 0], d_gps[:, 1],
        color="b", scale=scale, label="Observed"
    )
    ax1.quiver(
        llh_gps[:, 1], llh_gps[:, 0],
        m_gps[:, 0], m_gps[:, 1],
        color="r", scale=scale, label="Predicted"
    )
    if ndim == 3:
        sc = ax1.scatter(
            llh_gps[:, 1], llh_gps[:, 0],
            c=d_gps[:, 2], cmap=plt.cm.jet, s=50,
            linewidths=0.1, edgecolors="black"
        )
        fig.colorbar(sc, ax=ax1, label="Vertical (mm)")
    ax1.set_title("GPS Observation vs Prediction")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.legend(loc="upper right", fontsize=8)

    # Residual
    ax2 = fig.add_subplot(122)
    r_e = d_gps[:, 0] - m_gps[:, 0]
    r_n = d_gps[:, 1] - m_gps[:, 1]
    ax2.quiver(
        llh_gps[:, 1], llh_gps[:, 0],
        r_e, r_n,
        color="b", scale=scale
    )
    if ndim == 3:
        r_u = d_gps[:, 2] - m_gps[:, 2]
        sc2 = ax2.scatter(
            llh_gps[:, 1], llh_gps[:, 0],
            c=r_u, cmap=plt.cm.jet, linewidths=0.1, edgecolors="black"
        )
        fig.colorbar(sc2, ax=ax2, label="Vertical residual (mm)")
    ax2.set_title("GPS Residual")
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")

    return fig


def plot_sar_obs_mod(data, sol, green=None) -> Figure:
    """
    Plot InSAR observation vs prediction as scatter maps.

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(12, 5))
    fig.set_tight_layout(True)

    if data.llh_sar.size == 0:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No InSAR data", ha="center", va="center", fontsize=16)
        return fig

    llh_sar = data.llh_sar
    d_sar = data.d_sar
    ndim = data.ndim

    # offset into dhat for SAR data
    n_gps_comp = ndim * (len(data.llh_gps) if data.llh_gps.size > 0 else 0)
    n_lev_comp = len(data.d_lev) if data.d_lev.size > 0 else 0
    offset = n_gps_comp + n_lev_comp
    m_sar = sol.dhat[offset : offset + len(llh_sar)]

    # Observation
    ax1 = fig.add_subplot(121)
    sc1 = ax1.scatter(llh_sar[:, 1], llh_sar[:, 0], c=d_sar, cmap=plt.cm.jet, s=2)
    fig.colorbar(sc1, ax=ax1, label="LOS changes (mm)")
    ax1.set_title("InSAR LOS Observation")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")

    # Prediction
    ax2 = fig.add_subplot(122)
    sc2 = ax2.scatter(llh_sar[:, 1], llh_sar[:, 0], c=m_sar, cmap=plt.cm.jet, s=2)
    fig.colorbar(sc2, ax=ax2, label="LOS changes (mm)")
    ax2.set_title("InSAR LOS Prediction")
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")

    return fig


def plot_residuals(data, sol, scale=200) -> Figure:
    """
    Plot residual analysis for GPS and InSAR data.

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(12, 5))
    fig.set_tight_layout(True)

    has_gps = data.llh_gps.size > 0
    has_sar = data.llh_sar.size > 0

    if not has_gps and not has_sar:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No data for residual analysis", ha="center", va="center")
        return fig

    dhat = sol.dhat

    if has_gps:
        llh_gps = data.llh_gps
        ndim = data.ndim
        d_gps = data.d_gps.reshape(len(llh_gps), ndim)
        m_gps = dhat[: ndim * len(llh_gps)].reshape(len(llh_gps), ndim)
        r_e = d_gps[:, 0] - m_gps[:, 0]
        r_n = d_gps[:, 1] - m_gps[:, 1]

        ax1 = fig.add_subplot(121 if has_sar else 111)
        ax1.quiver(llh_gps[:, 1], llh_gps[:, 0], r_e, r_n, color="b", scale=scale)
        if ndim == 3:
            r_u = d_gps[:, 2] - m_gps[:, 2]
            sc = ax1.scatter(
                llh_gps[:, 1], llh_gps[:, 0],
                c=r_u, cmap=plt.cm.jet, linewidths=0.1, edgecolors="black"
            )
            fig.colorbar(sc, ax=ax1, label="Vert. residual (mm)")
        ax1.set_title("GPS Residual")
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")

    if has_sar:
        llh_sar = data.llh_sar
        d_sar = data.d_sar
        ndim = data.ndim
        n_gps_comp = ndim * (len(data.llh_gps) if data.llh_gps.size > 0 else 0)
        n_lev_comp = len(data.d_lev) if data.d_lev.size > 0 else 0
        offset = n_gps_comp + n_lev_comp
        m_sar = dhat[offset : offset + len(llh_sar)]
        r_sar = d_sar - m_sar

        ax2 = fig.add_subplot(122)
        sc2 = ax2.scatter(llh_sar[:, 1], llh_sar[:, 0], c=r_sar, cmap=plt.cm.RdBu_r, s=2)
        fig.colorbar(sc2, ax=ax2, label="LOS residual (mm)")
        ax2.set_title("InSAR Residual")
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")

    return fig


def plot_slip_3d(flt, sol, elevation=40, azimuth=-81, coordtype="llh") -> Figure:
    """
    Plot 3D slip distribution on fault patches.

    Parameters
    ----------
    flt : Fault (with .FaultGeom2AllVertex() method)
    sol : GeoInversion (with .slip attribute)
    elevation : float, view elevation angle
    azimuth : float, view azimuth angle
    coordtype : str, 'llh' or 'enu'

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev=elevation, azim=azimuth)

    # Get fault element vertices
    felem = flt.FaultGeom2AllVertex()
    slip = sol.slip[-1, : flt.nf * 3].reshape(flt.nf, 3)
    ts = np.sqrt(slip[:, 0] ** 2 + slip[:, 1] ** 2)
    ts_max = ts.max() if ts.max() > 0 else 1.0
    ts_norm = ts / ts_max

    if coordtype == "llh":
        x_cols = ["lt_lon", "rt_lon", "rb_lon", "lb_lon"]
        y_cols = ["lt_lat", "rt_lat", "rb_lat", "lb_lat"]
        z_cols = ["lt_dep", "rt_dep", "rb_dep", "lb_dep"]
    else:
        x_cols = ["lt_e", "rt_e", "rb_e", "lb_e"]
        y_cols = ["lt_n", "rt_n", "rb_n", "lb_n"]
        z_cols = ["lt_u", "rt_u", "rb_u", "lb_u"]

    x = felem[x_cols].values
    y = felem[y_cols].values
    z = felem[z_cols].values

    for i in range(len(felem)):
        verts = [list(zip(x[i], y[i], z[i]))]
        subfault = Poly3DCollection(verts, facecolors="w", linewidths=1, alpha=0.3)
        subfault.set_facecolor(plt.cm.jet(ts_norm[i]))
        subfault.set_linewidth(0.5)
        ax.add_collection3d(subfault)

    # Scatter for colorbar reference
    if coordtype == "llh":
        sc = ax.scatter(
            felem["lt_lon"], felem["lt_lat"],
            c=ts / 1000.0, cmap=plt.cm.jet, s=0.01, lw=0
        )
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    else:
        sc = ax.scatter(
            felem["lt_e"], felem["lt_n"],
            c=ts / 1000.0, cmap=plt.cm.jet, s=0.01, lw=0
        )
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")

    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
    cb.set_label("Total slip (m)")
    ax.set_zlabel("Depth (km)")
    ax.set_title("3D Slip Distribution")

    return fig


def plot_l_curve(sol) -> Figure:
    """
    Plot L-curve (roughness vs misfit) and individual misfit components.

    Parameters
    ----------
    sol : GeoInversion (with .smoothness, .misfit, .misfit_gps, .misfit_sar, .smo_facts)

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(12, 5))
    fig.set_tight_layout(True)

    smoothness = sol.smoothness
    misfit = sol.misfit
    smo_facts = sol.smo_facts
    misfit_gps = sol.misfit_gps
    misfit_sar = sol.misfit_sar

    if len(smo_facts) < 1:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "Not enough smoothing factors for L-curve",
                ha="center", va="center")
        return fig

    # L-Curve
    ax1 = fig.add_subplot(121)
    ax1.plot(smoothness, misfit, "bo-", markersize=6)
    if len(smo_facts) > 1:
        for i, sf in enumerate(smo_facts):
            ax1.annotate(f"{sf:.2f}", (smoothness[i], misfit[i]),
                         fontsize=7, textcoords="offset points", xytext=(5, 5))
    ax1.set_xlabel("Roughness")
    ax1.set_ylabel("Misfit (WRSS)")
    ax1.set_title("L-Curve")
    ax1.grid(True, alpha=0.3)

    # Misfit components
    ax2 = fig.add_subplot(122)
    width = 0.25
    x_pos = np.arange(len(smo_facts))
    has_gps = misfit_gps is not None and np.any(misfit_gps > 0)
    has_sar = misfit_sar is not None and np.any(misfit_sar > 0)

    bars = [ax2.bar(x_pos, misfit, width, label="Total")]
    if has_gps:
        bars.append(ax2.bar(x_pos + width, misfit_gps, width, label="GPS"))
    if has_sar:
        bars.append(ax2.bar(x_pos + 2 * width, misfit_sar, width, label="SAR"))

    ax2.set_xticks(x_pos + width)
    ax2.set_xticklabels([f"{sf:.2f}" for sf in smo_facts], fontsize=8)
    ax2.set_xlabel("Smoothing Factor")
    ax2.set_ylabel("WRSS")
    ax2.set_title("Misfit vs Smoothing Factor")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3, axis="y")

    return fig


def plot_data_summary(data) -> Figure:
    """
    Plot a summary of loaded observation data.

    Returns
    -------
    Figure
    """
    fig = Figure(figsize=(12, 5))
    fig.set_tight_layout(True)

    has_gps = data.llh_gps.size > 0
    has_sar = data.llh_sar.size > 0

    if not has_gps and not has_sar:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No data loaded", ha="center", va="center", fontsize=16)
        return fig

    n_plots = int(has_gps) + int(has_sar)
    plot_idx = 1

    if has_gps:
        ax1 = fig.add_subplot(1, n_plots, plot_idx)
        plot_idx += 1
        ndim = data.ndim
        d_gps = data.d_gps.reshape(len(data.llh_gps), ndim)
        ax1.quiver(
            data.llh_gps[:, 1], data.llh_gps[:, 0],
            d_gps[:, 0], d_gps[:, 1],
            color="b", scale=200
        )
        if ndim == 3:
            sc = ax1.scatter(
                data.llh_gps[:, 1], data.llh_gps[:, 0],
                c=d_gps[:, 2], cmap=plt.cm.jet, s=30,
                linewidths=0.1, edgecolors="black"
            )
            fig.colorbar(sc, ax=ax1, label="Vertical (mm)")
        ax1.set_title(f"GPS Data ({len(data.llh_gps)} stations)")
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")

    if has_sar:
        ax2 = fig.add_subplot(1, n_plots, plot_idx)
        sc2 = ax2.scatter(
            data.llh_sar[:, 1], data.llh_sar[:, 0],
            c=data.d_sar, cmap=plt.cm.jet, s=2
        )
        fig.colorbar(sc2, ax=ax2, label="LOS (mm)")
        ax2.set_title(f"InSAR Data ({len(data.llh_sar)} points)")
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")

    return fig
