"""
PyGINVC GUI Sidebar - Configuration panels for all inversion parameters.

This module renders the sidebar form controls that correspond to the YAML
configuration sections: dict_data, dict_fault, dict_green, dict_weight,
dict_bound, dict_export.

Data files are uploaded via file selectors (outside the form for immediate
processing). Uploaded files are saved to a temporary directory and their
paths are passed to the inversion runner.
"""
import os
import tempfile
import streamlit as st
import yaml


def _temp_dir():
    """Get or create the session-scoped temporary directory for uploaded files."""
    if "_temp_dir" not in st.session_state:
        st.session_state["_temp_dir"] = tempfile.mkdtemp(prefix="pyginvc_")
    return st.session_state["_temp_dir"]


def _save_uploaded_file(uploaded_file, subdir="") -> str:
    """
    Save a single UploadedFile to the session temp directory.

    Returns the absolute path of the saved file.
    """
    base = _temp_dir()
    if subdir:
        base = os.path.join(base, subdir)
        os.makedirs(base, exist_ok=True)
    # Keep original filename, avoid collisions
    fname = uploaded_file.name
    dst = os.path.join(base, fname)
    # If a file with the same name was already saved this session, overwrite
    try:
        with open(dst, "wb") as f:
            f.write(uploaded_file.getbuffer())
    except Exception:
        # getbuffer() can fail for large files; fall back to read()
        with open(dst, "wb") as f:
            f.write(uploaded_file.read())
    return dst


def _save_uploaded_files(uploaded_files, subdir="") -> list:
    """
    Save multiple UploadedFile objects and return a list of absolute paths.
    """
    if not uploaded_files:
        return []
    paths = []
    for uf in uploaded_files:
        path = _save_uploaded_file(uf, subdir)
        paths.append(path)
    return paths


def render_sidebar():
    """
    Render the complete sidebar with all configuration sections.
    Returns the assembled config dict when the form is submitted.
    """
    with st.sidebar:
        st.title("PyGINVC")

        # ---- All widgets inside ONE form ----
        with st.form("inversion_form"):
            # File uploaders (inside form for reliable state management)
            gps_files_val, sar_files_val, lev_files_val, gfiletype_val = _render_data_uploaders()
            fault_file_val, nsegs_val, ndeps_val, doSubFault_val = _render_fault_uploader()

            st.divider()

            cfg = st.session_state.get("cfg", {})
            dict_green = _render_green_section(cfg.get("dict_green", {}))
            dict_weight = _render_weight_section(cfg.get("dict_weight", {}))
            dict_bound = _render_bound_section(cfg.get("dict_bound", {}))
            dict_export = _render_export_section(cfg.get("dict_export", {}))

            st.divider()
            col1, col2 = st.columns(2)
            load_data = col1.form_submit_button("Load Data", use_container_width=True)
            run_inv = col2.form_submit_button(
                "Run Inversion", use_container_width=True, type="primary"
            )

        # Download config (outside form, uses session-stored cfg)
        saved_cfg = st.session_state.get("cfg", {})
        if saved_cfg:
            yaml_str = yaml.dump(saved_cfg, default_flow_style=False, sort_keys=False)
            st.download_button(
                "Export Config (YAML)",
                data=yaml_str,
                file_name="pyginvc_config.yaml",
                mime="text/yaml",
                use_container_width=True,
            )

    # Store actions
    st.session_state["_load_data"] = load_data
    st.session_state["_run_inversion"] = run_inv

    if load_data or run_inv:
        # Save uploaded files to temp on form submission
        gps_paths = _save_uploaded_files(gps_files_val, "gps")
        sar_paths = _save_uploaded_files(sar_files_val, "sar")
        lev_paths = _save_uploaded_files(lev_files_val, "lev")
        fault_path = ""
        if fault_file_val is not None:
            fault_path = _save_uploaded_file(fault_file_val, "fault")

        dict_data = {
            "gpsfile": gps_paths,
            "sarfile": sar_paths,
            "levfile": lev_paths,
            "gfiletype": gfiletype_val,
        }
        dict_fault = {
            "faultfile": fault_path,
            "nsegs": nsegs_val,
            "ndeps": ndeps_val,
            "doSubFault": doSubFault_val,
        }
        assembled = {
            "dict_data": dict_data,
            "dict_fault": dict_fault,
            "dict_green": dict_green,
            "dict_weight": dict_weight,
            "dict_bound": dict_bound,
            "dict_export": dict_export,
        }
        st.session_state["cfg"] = assembled
        return assembled
    return st.session_state.get("cfg")


def _render_data_uploaders():
    """Render GPS/SAR/LEV file uploaders (inside form).

    Returns (gps_files, sar_files, lev_files, gfiletype)
    where gps_files etc. are lists of UploadedFile or None/[] if none selected.
    """
    st.subheader("Data Files")

    gfiletype = st.selectbox(
        "GPS Data Format",
        options=["GMT2D", "GMT3D", "IOS2D", "IOS3D"],
        index=0,
        key="gfiletype",
    )

    # GPS files (multi)
    gps_files = st.file_uploader(
        "GPS Data Files",
        type=["gmtvec", "dat", "txt", "vec"],
        accept_multiple_files=True,
        key="gps_uploader",
        help="Upload one or more GPS velocity/displacement files.",
    )
    if gps_files:
        with st.expander(f"GPS: {len(gps_files)} file(s) selected", expanded=False):
            for uf in gps_files:
                st.caption(uf.name)

    # InSAR files (multi)
    sar_files = st.file_uploader(
        "InSAR Data Files",
        type=["dat", "txt", "grd"],
        accept_multiple_files=True,
        key="sar_uploader",
        help="Upload one or more InSAR LOS displacement files.",
    )
    if sar_files:
        with st.expander(f"InSAR: {len(sar_files)} file(s) selected", expanded=False):
            for uf in sar_files:
                st.caption(uf.name)

    # Leveling files (multi)
    lev_files = st.file_uploader(
        "Leveling Data Files",
        type=["dat", "txt"],
        accept_multiple_files=True,
        key="lev_uploader",
        help="Upload one or more leveling data files.",
    )
    if lev_files:
        with st.expander(f"Leveling: {len(lev_files)} file(s) selected", expanded=False):
            for uf in lev_files:
                st.caption(uf.name)

    return gps_files, sar_files, lev_files, gfiletype


def _render_fault_uploader():
    """Render fault geometry file uploader and grid params (inside form).

    Returns (fault_file, nsegs, ndeps, doSubFault)
    where fault_file is an UploadedFile or None.
    """
    st.subheader("Fault Geometry")

    fault_file = st.file_uploader(
        "Fault Geometry File",
        type=[],  # accept any extension
        accept_multiple_files=False,
        key="fault_uploader",
        help="Upload the fault geometry definition file (FaultGeom format).",
    )
    if fault_file is not None:
        st.caption(f"Selected: {fault_file.name}")

    col1, col2 = st.columns(2)
    nsegs = col1.number_input(
        "nsegs (strike)",
        min_value=1,
        max_value=500,
        value=30,
        step=1,
        key="nsegs",
    )
    ndeps = col2.number_input(
        "ndeps (dip)",
        min_value=1,
        max_value=500,
        value=20,
        step=1,
        key="ndeps",
    )
    doSubFault = st.checkbox(
        "Subdivide Fault into Patches",
        value=True,
        key="doSubFault",
    )

    return fault_file, nsegs, ndeps, doSubFault


def _render_green_section(defaults: dict) -> dict:
    """Render dict_green configuration section."""
    st.subheader("Green's Function")

    grnmethod = st.selectbox(
        "Method",
        options=["okada", "wang"],
        index=["okada", "wang"].index(defaults.get("grnmethod", "okada")),
        key="grnmethod",
    )

    col1, col2 = st.columns(2)
    nu = col1.slider(
        "Poisson Ratio (nu)",
        min_value=0.0,
        max_value=0.5,
        value=float(defaults.get("nu", 0.25)),
        step=0.01,
        key="nu",
    )
    modulus = col2.number_input(
        "Shear Modulus (Pa)",
        value=float(defaults.get("modulus", 3e10)),
        format="%.1e",
        key="modulus",
    )

    st.write("Green's Function Components:")
    greentype_default = defaults.get("greentype", [1, 0, 0])
    if not isinstance(greentype_default, list) or len(greentype_default) != 3:
        greentype_default = [1, 0, 0]

    col1, col2, col3 = st.columns(3)
    ss_comp = col1.checkbox("Strike-Slip", value=bool(greentype_default[0]), key="gt_ss")
    ds_comp = col2.checkbox("Dip-Slip", value=bool(greentype_default[1]), key="gt_ds")
    op_comp = col3.checkbox("Opening", value=bool(greentype_default[2]), key="gt_op")

    st.write("Boundary Conditions [top, bottom, end, start]:")
    bcs_default = defaults.get("bcs", [0, 0, 0, 0])
    if not isinstance(bcs_default, list) or len(bcs_default) != 4:
        bcs_default = [0, 0, 0, 0]

    bcs_cols = st.columns(4)
    bcs_labels = ["Top", "Bottom", "End", "Start"]
    bcs = []
    for i, (col, label) in enumerate(zip(bcs_cols, bcs_labels)):
        val = col.selectbox(
            label,
            options=[0, 1],
            index=int(bcs_default[i]) if bcs_default[i] in [0, 1] else 0,
            key=f"bcs_{i}",
        )
        bcs.append(val)

    greenfile = st.selectbox(
        "Green's Function File",
        options=["SAVE", "Load from file"],
        index=0 if defaults.get("greenfile", "SAVE") == "SAVE" else 1,
        key="greenfile_mode",
    )
    if greenfile == "Load from file":
        greenfile_path = st.text_input(
            "Green's Function Path",
            value=defaults.get("greenfile", ""),
            key="greenfile_path",
        )
    else:
        greenfile_path = "SAVE"

    col1, col2 = st.columns(2)
    gps_ramp = col1.checkbox(
        "GPS Ramp", value=defaults.get("gps_ramp", True), key="gps_ramp"
    )
    sar_ramp = col2.checkbox(
        "SAR Ramp", value=defaults.get("sar_ramp", True), key="sar_ramp"
    )
    rake_beta = st.number_input(
        "Rake Beta",
        value=float(defaults.get("rake_beta", 0)),
        key="rake_beta",
    )

    return {
        "grnmethod": grnmethod,
        "nu": nu,
        "modulus": modulus,
        "greentype": [int(ss_comp), int(ds_comp), int(op_comp)],
        "bcs": bcs,
        "greenfile": greenfile_path,
        "gps_ramp": gps_ramp,
        "sar_ramp": sar_ramp,
        "rake_beta": rake_beta,
    }


def _render_weight_section(defaults: dict) -> dict:
    """Render dict_weight configuration section."""
    st.subheader("Weight & Smoothing")

    sf_default = defaults.get("smoothfactor", [0.1, 0.5, 0.1])
    if not isinstance(sf_default, list) or len(sf_default) != 3:
        sf_default = [0.1, 0.5, 0.1]

    st.write("Smoothing Factor [start, end, step]:")
    sf_cols = st.columns(3)
    sf_start = sf_cols[0].number_input(
        "Start", value=float(sf_default[0]), step=0.01, key="sf_start"
    )
    sf_end = sf_cols[1].number_input(
        "End", value=float(sf_default[1]), step=0.01, key="sf_end"
    )
    sf_step = sf_cols[2].number_input(
        "Step", value=float(sf_default[2]), step=0.01, key="sf_step"
    )

    st.write("Data Weights:")
    wgps_default = defaults.get("wgps", [1])
    wsar_default = defaults.get("wsar", [1])
    wlev_default = defaults.get("wlev", [1])
    wgps_val = wgps_default[0] if isinstance(wgps_default, list) else wgps_default
    wsar_val = wsar_default[0] if isinstance(wsar_default, list) else wsar_default
    wlev_val = wlev_default[0] if isinstance(wlev_default, list) else wlev_default

    w_cols = st.columns(3)
    wgps = w_cols[0].number_input(
        "GPS Weight", value=float(wgps_val), step=0.1, key="wgps"
    )
    wsar = w_cols[1].number_input(
        "SAR Weight", value=float(wsar_val), step=0.1, key="wsar"
    )
    wlev = w_cols[2].number_input(
        "LEV Weight", value=float(wlev_val), step=0.1, key="wlev"
    )

    return {
        "smoothfactor": [sf_start, sf_end, sf_step],
        "wgps": [wgps],
        "wsar": [wsar],
        "wlev": [wlev],
    }


def _render_bound_section(defaults: dict) -> dict:
    """Render dict_bound configuration section."""
    st.subheader("Boundary Conditions")

    bound_switch = st.checkbox(
        "Enable Bounds", value=defaults.get("bound_switch", True), key="bound_switch"
    )

    if bound_switch:
        ss_default = defaults.get("ss_range", [-100, 100])
        ds_default = defaults.get("ds_range", [0, 2000])
        op_default = defaults.get("op_range", [0, 0.01])

        if not isinstance(ss_default, list) or len(ss_default) != 2:
            ss_default = [-100, 100]
        if not isinstance(ds_default, list) or len(ds_default) != 2:
            ds_default = [0, 2000]
        if not isinstance(op_default, list) or len(op_default) != 2:
            op_default = [0, 0.01]

        st.write("Slip Ranges [min, max] (mm):")
        st.write("**Strike-Slip:**")
        ss_cols = st.columns(2)
        ss_min = ss_cols[0].number_input(
            "Min", value=float(ss_default[0]), key="ss_min"
        )
        ss_max = ss_cols[1].number_input(
            "Max", value=float(ss_default[1]), key="ss_max"
        )

        st.write("**Dip-Slip:**")
        ds_cols = st.columns(2)
        ds_min = ds_cols[0].number_input(
            "Min", value=float(ds_default[0]), key="ds_min"
        )
        ds_max = ds_cols[1].number_input(
            "Max", value=float(ds_default[1]), key="ds_max"
        )

        st.write("**Opening:**")
        op_cols = st.columns(2)
        op_min = op_cols[0].number_input(
            "Min", value=float(op_default[0]), key="op_min"
        )
        op_max = op_cols[1].number_input(
            "Max", value=float(op_default[1]), key="op_max"
        )

        slip_lb = st.text_input(
            "Lower Bound File (optional)",
            value=defaults.get("slip_lb", ""),
            key="slip_lb",
        )
        slip_ub = st.text_input(
            "Upper Bound File (optional)",
            value=defaults.get("slip_ub", ""),
            key="slip_ub",
        )
    else:
        ss_min, ss_max = -100, 100
        ds_min, ds_max = 0, 2000
        op_min, op_max = 0, 0.01
        slip_lb, slip_ub = "", ""

    return {
        "bound_switch": bound_switch,
        "ss_range": [ss_min, ss_max],
        "ds_range": [ds_min, ds_max],
        "op_range": [op_min, op_max],
        "slip_lb": slip_lb,
        "slip_ub": slip_ub,
        "sar_plane_range": defaults.get("sar_plane_range", []),
    }


def _render_export_section(defaults: dict) -> dict:
    """Render dict_export configuration section."""
    st.subheader("Export Settings")

    outpath = st.text_input(
        "Output Directory",
        value=defaults.get("outpath", "./output"),
        key="outpath",
    )
    view = st.checkbox(
        "Show Results After Inversion",
        value=defaults.get("view", True),
        key="view",
    )
    vecscale = st.number_input(
        "Vector Scale",
        value=int(defaults.get("vecscale", 5000)),
        step=100,
        key="vecscale",
    )

    return {
        "outpath": outpath,
        "view": view,
        "vecscale": vecscale,
    }


def validate(cfg: dict) -> tuple:
    """
    Validate the assembled configuration dictionary.

    Returns
    -------
    tuple : (is_valid: bool, errors: list[str])
    """
    errors = []

    # Data files (now lists of paths)
    dd = cfg.get("dict_data", {})
    gpsfiles = dd.get("gpsfile", [])
    sarfiles = dd.get("sarfile", [])
    levfiles = dd.get("levfile", [])

    # Normalize: ensure lists
    if isinstance(gpsfiles, str):
        gpsfiles = [gpsfiles] if gpsfiles else []
    if isinstance(sarfiles, str):
        sarfiles = [sarfiles] if sarfiles else []
    if isinstance(levfiles, str):
        levfiles = [levfiles] if levfiles else []

    # Check each GPS file exists
    for f in gpsfiles:
        if f and not os.path.isfile(f):
            errors.append(f"GPS file does not exist: {f}")
    for f in sarfiles:
        if f and not os.path.isfile(f):
            errors.append(f"SAR file does not exist: {f}")
    for f in levfiles:
        if f and not os.path.isfile(f):
            errors.append(f"Leveling file does not exist: {f}")

    # At least one data source
    has_gps = any(f for f in gpsfiles)
    has_sar = any(f for f in sarfiles)
    has_lev = any(f for f in levfiles)
    if not has_gps and not has_sar and not has_lev:
        errors.append("At least one data file (GPS, SAR, or LEV) must be uploaded.")

    # Fault
    df = cfg.get("dict_fault", {})
    faultfile = df.get("faultfile", "")
    if not faultfile or not os.path.isfile(faultfile):
        errors.append(f"Fault geometry file is required: {faultfile}")

    nsegs = df.get("nsegs", 0)
    ndeps = df.get("ndeps", 0)
    if nsegs < 1 or ndeps < 1:
        errors.append("nsegs and ndeps must be >= 1")

    # Green's function
    dg = cfg.get("dict_green", {})
    nu = dg.get("nu", 0)
    if not (0 < nu < 0.5):
        errors.append(f"Poisson ratio nu must be in (0, 0.5), got {nu}")

    modulus = dg.get("modulus", 0)
    if modulus <= 0:
        errors.append(f"Shear modulus must be > 0, got {modulus}")

    # Smoothing factor
    dw = cfg.get("dict_weight", {})
    sf = dw.get("smoothfactor", [0, 0, 0])
    if len(sf) == 3 and sf[2] <= 0:
        errors.append(f"Smoothing factor step must be > 0, got {sf[2]}")

    # Boundary
    db = cfg.get("dict_bound", {})
    if db.get("bound_switch", False):
        ss = db.get("ss_range", [0, 0])
        if len(ss) == 2 and ss[0] > ss[1]:
            errors.append(f"Strike-slip range min ({ss[0]}) > max ({ss[1]})")
        ds = db.get("ds_range", [0, 0])
        if len(ds) == 2 and ds[0] > ds[1]:
            errors.append(f"Dip-slip range min ({ds[0]}) > max ({ds[1]})")
        op = db.get("op_range", [0, 0])
        if len(op) == 2 and op[0] > op[1]:
            errors.append(f"Opening range min ({op[0]}) > max ({op[1]})")

    return (len(errors) == 0, errors)
