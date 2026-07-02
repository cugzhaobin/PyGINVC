"""
PyGINVC GUI - Streamlit-based web application for fault slip inversion.

This is the main entry point. Launch with:

    streamlit run src/pyginvc/gui/app.py

Or after ``pip install -e .``:

    pyginvc-gui
"""
import sys
import os
import subprocess

import numpy as np
import streamlit as st

from pyginvc.gui.sidebar import render_sidebar, validate
from pyginvc.gui.utils.inversion_runner import run_load_data, run_inversion
from pyginvc.gui.utils.log_handler import setup_streamlit_logging
from pyginvc.gui.utils.plot_helpers import (
    plot_gps_obs_mod,
    plot_sar_obs_mod,
    plot_slip_3d,
    plot_residuals,
    plot_l_curve,
    plot_data_summary,
)


def _init_session_state():
    """Initialize session_state with default values on first run."""
    defaults = {
        "cfg": {},
        "data": None,
        "flt": None,
        "green": None,
        "lap": None,
        "sol": None,
        "log_messages": [],
        "inversion_done": False,
        "data_loaded": False,
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def _render_main_area():
    """Render the main content area with tabs for visualization."""

    if st.session_state.get("inversion_done"):
        _render_inversion_results()
    elif st.session_state.get("data_loaded"):
        _render_data_preview()
    else:
        _render_welcome()


def _render_welcome():
    """Render the welcome/landing page."""
    st.title("PyGINVC - Fault Slip Inversion")
    st.markdown(
        """
        Welcome to the **PyGINVC** graphical interface for geodetic fault slip inversion.

        ### How to use:
        1. **Upload data files** in the sidebar — GPS, InSAR, and/or Leveling
        2. **Upload a fault geometry file** and configure patch subdivision
        3. Adjust **Green's function**, **smoothing**, and **boundary** parameters as needed
        4. Click **"Load Data"** to preview your observation data and fault geometry
        5. Click **"Run Inversion"** to execute the slip inversion and view results

        ### Supported data types:
        - **GNSS/GPS**: 2D or 3D velocity/displacement vectors (GMT, IOS formats)
        - **InSAR**: Line-of-sight displacement
        - **Leveling**: Vertical displacement

        ### Supported models:
        - Rectangular fault patches (Okada or layered Earth Green's functions)
        - Constrained or unconstrained linear inversion
        - Laplacian smoothing with adjustable factors
        """
    )


def _render_data_preview():
    """Render data preview after loading."""
    data = st.session_state.get("data")
    flt = st.session_state.get("flt")

    st.header("Data Preview")

    # ---- Stats cards ----
    n_gps = len(data.llh_gps) if data.llh_gps.size > 0 else 0
    n_sar = len(data.llh_sar) if data.llh_sar.size > 0 else 0
    n_lev = len(data.llh_lev) if data.llh_lev.size > 0 else 0
    n_patches = flt.nf if flt is not None else 0

    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("GPS Stations", n_gps)
    col2.metric("InSAR Points", n_sar)
    col3.metric("Leveling", n_lev)
    col4.metric("Fault Patches", n_patches)
    col5.metric("GPS Dim", f"{data.ndim}D" if data else "-")

    st.divider()

    # ---- Main content: plots on left, info on right ----
    left, right = st.columns([3, 1])

    with left:
        if data is not None:
            fig = plot_data_summary(data)
            st.pyplot(fig)
        else:
            st.info("No data loaded.")

    with right:
        if flt is not None:
            st.subheader("Fault Geometry")
            st.write(f"**Patches:** {flt.nf}")
            st.write(f"**Grid:** {flt.nsegs} (strike) x {flt.ndeps} (dip)")
            if flt.origin is not None and (isinstance(flt.origin, (list, tuple)) and len(flt.origin) > 0 or hasattr(flt.origin, 'size') and flt.origin.size > 0):
                st.write(f"**Origin:** {flt.origin}")

        st.subheader("Loaded Files")
        cfg = st.session_state.get("cfg", {})
        dd = cfg.get("dict_data", {})
        gpsfiles = dd.get("gpsfile", [])
        sarfiles = dd.get("sarfile", [])
        levfiles = dd.get("levfile", [])
        faultfile = cfg.get("dict_fault", {}).get("faultfile", "")

        if faultfile:
            st.caption(f"**Fault:** {os.path.basename(faultfile)}")
        if isinstance(gpsfiles, list) and gpsfiles:
            st.caption(f"**GPS ({len(gpsfiles)}):**")
            for p in gpsfiles:
                st.caption(f"  {os.path.basename(p)}")
        elif gpsfiles:
            st.caption(f"**GPS:** {os.path.basename(str(gpsfiles))}")
        if isinstance(sarfiles, list) and sarfiles:
            st.caption(f"**InSAR ({len(sarfiles)}):**")
            for p in sarfiles:
                st.caption(f"  {os.path.basename(p)}")
        elif sarfiles:
            st.caption(f"**InSAR:** {os.path.basename(str(sarfiles))}")
        if isinstance(levfiles, list) and levfiles:
            st.caption(f"**LEV ({len(levfiles)}):**")
            for p in levfiles:
                st.caption(f"  {os.path.basename(p)}")
        elif levfiles:
            st.caption(f"**LEV:** {os.path.basename(str(levfiles))}")

    # ---- File content preview (full width below plots) ----
    st.divider()
    _render_file_contents(gpsfiles, sarfiles, levfiles, faultfile)


def _render_file_contents(gpsfiles, sarfiles, levfiles, faultfile):
    """Render collapsible file content previews for uploaded data files."""
    st.subheader("File Contents")

    def _show_file(label, path, max_lines=50):
        """Display file content in an expander, truncated to max_lines."""
        if not path or not os.path.isfile(path):
            return
        try:
            with open(path, "r") as f:
                lines = f.readlines()
            total = len(lines)
            preview = "".join(lines[:max_lines])
            if total > max_lines:
                preview += f"\n... ({total - max_lines} more lines, {total} total)"
            with st.expander(f"{label} ({total} lines)", expanded=False):
                st.code(preview, language=None)
        except Exception as e:
            with st.expander(f"{label}", expanded=False):
                st.caption(f"Cannot read: {e}")

    # Fault geometry file
    if faultfile:
        _show_file(f"Fault: {os.path.basename(faultfile)}", faultfile)

    # GPS files
    if isinstance(gpsfiles, list):
        for p in gpsfiles:
            _show_file(f"GPS: {os.path.basename(p)}", p)
    elif gpsfiles:
        _show_file(f"GPS: {os.path.basename(str(gpsfiles))}", str(gpsfiles))

    # InSAR files
    if isinstance(sarfiles, list):
        for p in sarfiles:
            _show_file(f"InSAR: {os.path.basename(p)}", p)
    elif sarfiles:
        _show_file(f"InSAR: {os.path.basename(str(sarfiles))}", str(sarfiles))

    # Leveling files
    if isinstance(levfiles, list):
        for p in levfiles:
            _show_file(f"LEV: {os.path.basename(p)}", p)
    elif levfiles:
        _show_file(f"LEV: {os.path.basename(str(levfiles))}", str(levfiles))


def _render_inversion_results():
    """Render inversion results with dual-column layout."""
    data = st.session_state.get("data")
    flt = st.session_state.get("flt")
    sol = st.session_state.get("sol")

    st.header("Inversion Results")

    # ---- Enhanced metrics bar ----
    cols = st.columns(6)
    if hasattr(sol, "misfit") and sol.misfit is not None:
        cols[0].metric("Final Misfit", f"{sol.misfit[-1]:.2f}")
    if hasattr(sol, "moment") and sol.moment is not None:
        Mo, Mw = sol.moment[-1]
        cols[1].metric("Mw", f"{Mw:.2f}")
        cols[2].metric("M0 (N\u00b7m)", f"{Mo:.2e}")
    if hasattr(sol, "smo_facts"):
        cols[3].metric("Smoothing Factors", len(sol.smo_facts))
        # Best smoothing factor (lowest misfit)
        if len(sol.misfit) > 0:
            best_idx = int(np.argmin(sol.misfit))
            cols[4].metric("Best Smooth.", f"{sol.smo_facts[best_idx]:.3f}")
    if hasattr(sol, "misfit_gps") and sol.misfit_gps is not None and len(sol.misfit_gps) > 0:
        cols[5].metric("GPS Residual", f"{sol.misfit_gps[-1]:.2f}")

    st.divider()

    # ---- Dual-column layout ----
    row1_left, row1_right = st.columns(2)

    with row1_left:
        st.subheader("GPS: Observation vs Model")
        fig = plot_gps_obs_mod(data, sol)
        st.pyplot(fig)

    with row1_right:
        st.subheader("InSAR: Observation vs Model")
        fig = plot_sar_obs_mod(data, sol)
        st.pyplot(fig)

    st.divider()

    row2_left, row2_right = st.columns(2)

    with row2_left:
        st.subheader("3D Slip Distribution")
        col_a, col_b = st.columns(2)
        elev = col_a.slider("Elevation", 0, 90, 40, key="slip_elev")
        azim = col_b.slider("Azimuth", -180, 180, -81, key="slip_azim")
        coordtype = st.radio("Coordinates", ["llh", "enu"], horizontal=True, key="slip_coord")
        fig = plot_slip_3d(flt, sol, elevation=elev, azimuth=azim, coordtype=coordtype)
        st.pyplot(fig)

    with row2_right:
        st.subheader("L-Curve")
        fig = plot_l_curve(sol)
        st.pyplot(fig)

    st.divider()

    # ---- Full-width residuals ----
    st.subheader("Residual Analysis")
    fig = plot_residuals(data, sol)
    st.pyplot(fig)

    st.divider()

    # ---- Download + Log ----
    col_a, col_b = st.columns(2)
    with col_a:
        outpath = st.session_state.get("cfg", {}).get("dict_export", {}).get("outpath", "")
        sol_h5 = os.path.join(outpath, "solutions.h5") if outpath else ""
        if sol_h5 and os.path.isfile(sol_h5):
            with open(sol_h5, "rb") as f:
                st.download_button(
                    "Download solutions.h5",
                    data=f,
                    file_name="solutions.h5",
                    mime="application/octet-stream",
                )
        else:
            st.caption("Results saved to output directory.")
    with col_b:
        with st.expander("Log Messages", expanded=False):
            _render_log_tab()


def _render_log_tab():
    """Render the log messages tab."""
    logs = st.session_state.get("log_messages", [])
    if logs:
        st.code("\n".join(logs), language="log")
        st.download_button(
            "Download Log",
            data="\n".join(logs),
            file_name="pyginvc_log.txt",
            mime="text/plain",
        )
    else:
        st.info("No log messages yet.")


def _render_summary_stats():
    """Render summary statistics of loaded data and fault."""
    data = st.session_state.get("data")
    flt = st.session_state.get("flt")

    if data is not None:
        st.subheader("Data Summary")
        n_gps = len(data.llh_gps) if data.llh_gps.size > 0 else 0
        n_sar = len(data.llh_sar) if data.llh_sar.size > 0 else 0
        n_lev = len(data.llh_lev) if data.llh_lev.size > 0 else 0
        st.write(f"- GPS stations: **{n_gps}**")
        st.write(f"- InSAR points: **{n_sar}**")
        st.write(f"- Leveling stations: **{n_lev}**")
        st.write(f"- GPS dimension: **{data.ndim}**")

    if flt is not None:
        st.subheader("Fault Summary")
        st.write(f"- Total patches: **{flt.nf}**")
        st.write(f"- Grid: **{flt.nsegs}** x **{flt.ndeps}**")


def _handle_load_data(cfg: dict):
    """Execute the Load Data workflow."""
    valid, errors = validate(cfg)
    if not valid:
        for e in errors:
            st.error(e)
        return

    with st.status("Loading data...", expanded=True) as status:
        try:
            result = run_load_data(cfg, progress_callback=st.write)
            st.session_state["data"] = result["data"]
            st.session_state["flt"] = result["flt"]
            st.session_state["data_loaded"] = True
            status.update(label="Data loaded!", state="complete")
        except Exception as e:
            status.update(label="Error loading data", state="error")
            st.exception(e)


def _handle_run_inversion(cfg: dict):
    """Execute the Run Inversion workflow."""
    valid, errors = validate(cfg)
    if not valid:
        for e in errors:
            st.error(e)
        return

    data = st.session_state.get("data")
    flt = st.session_state.get("flt")

    # If data not loaded yet, load it first
    if data is None or flt is None:
        with st.status("Loading data...", expanded=True) as status:
            try:
                result = run_load_data(cfg, progress_callback=st.write)
                data = result["data"]
                flt = result["flt"]
                st.session_state["data"] = data
                st.session_state["flt"] = flt
                st.session_state["data_loaded"] = True
                status.update(label="Data loaded!", state="complete")
            except Exception as e:
                status.update(label="Error loading data", state="error")
                st.exception(e)
                return

    # Run inversion
    with st.status("Running inversion...", expanded=True) as status:
        try:
            result = run_inversion(cfg, data, flt, progress_callback=st.write)
            st.session_state["green"] = result["green"]
            st.session_state["lap"] = result["lap"]
            st.session_state["sol"] = result["sol"]
            st.session_state["inversion_done"] = True
            status.update(label="Inversion complete!", state="complete")
        except Exception as e:
            status.update(label="Inversion failed", state="error")
            st.exception(e)


def main():
    """Main entry point for the Streamlit application."""
    st.set_page_config(
        page_title="PyGINVC - Fault Slip Inversion",
        page_icon="🌍",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    _init_session_state()

    # Setup logging
    setup_streamlit_logging(st.session_state["log_messages"])

    # Render sidebar and get config
    cfg = render_sidebar()

    # Handle actions triggered by form submission
    if st.session_state.get("_load_data", False) and cfg is not None:
        st.session_state["_load_data"] = False
        _handle_load_data(cfg)

    if st.session_state.get("_run_inversion", False) and cfg is not None:
        st.session_state["_run_inversion"] = False
        _handle_run_inversion(cfg)

    # Render main content
    _render_main_area()


if __name__ == "__main__":
    # Allow running with: python -m pyginvc.gui.app
    # or: streamlit run src/pyginvc/gui/app.py
    import streamlit.runtime

    if streamlit.runtime.exists():
        main()
    else:
        # Not inside streamlit runtime; launch streamlit as a subprocess
        script_path = os.path.abspath(__file__)
        sys.exit(subprocess.run([sys.executable, "-m", "streamlit", "run", script_path]).returncode)
