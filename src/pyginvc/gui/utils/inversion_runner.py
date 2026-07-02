"""
Inversion Runner - Encapsulates the slip inversion workflow as a callable function.

This module mirrors the logic from ``pyginvc.scripts.SlipInversion.run_inv()``
and adapts it for use in the Streamlit GUI, providing a progress callback so that
``st.status`` can display step-by-step feedback.
"""
import os
import logging

from pyginvc.GeoData.GeoData import GeoData
from pyginvc.Geometry.Patch import Fault
from pyginvc.Greens.Okada import Okada
from pyginvc.Laplacian.RectPlaneLap import RectPlaneLap
from pyginvc.Inversion.GeoInversion import GeoInversion
from pyginvc.Export.Output import Output

logger = logging.getLogger(__name__)


def run_load_data(cfg: dict, progress_callback=None) -> dict:
    """
    Load geodetic data and fault geometry (Step 1 of the workflow).

    Parameters
    ----------
    cfg : dict
        Full configuration dictionary with keys dict_data, dict_fault, etc.
    progress_callback : callable, optional
        Called with ``progress_callback(message: str)`` to report progress.

    Returns
    -------
    dict : {"data": GeoData, "flt": Fault}
    """

    def report(msg):
        logger.info(msg)
        if progress_callback:
            progress_callback(msg)

    dict_data = cfg["dict_data"]
    dict_fault = cfg["dict_fault"]

    # ---- GeoData ----
    report("Loading geodetic data...")
    gpsfiles = dict_data.get("gpsfile", [])
    sarfiles = dict_data.get("sarfile", [])
    levfiles = dict_data.get("levfile", [])
    gfiletype = dict_data.get("gfiletype", "GMT2D")

    # Normalize to lists (supports single string or list from config)
    if isinstance(gpsfiles, str):
        gpsfiles = [gpsfiles] if gpsfiles else []
    if isinstance(sarfiles, str):
        sarfiles = [sarfiles] if sarfiles else []
    if isinstance(levfiles, str):
        levfiles = [levfiles] if levfiles else []

    data = GeoData(gpsfiles, sarfiles, levfiles, gfiletype)
    data.load_data()

    n_gps = len(data.llh_gps) if data.llh_gps.size > 0 else 0
    n_sar = len(data.llh_sar) if data.llh_sar.size > 0 else 0
    n_lev = len(data.llh_lev) if data.llh_lev.size > 0 else 0
    report(f"Data loaded: {n_gps} GPS, {n_sar} InSAR, {n_lev} Leveling stations.")

    # ---- Fault ----
    report("Loading fault geometry...")
    faultfile = dict_fault["faultfile"]
    nsegs = dict_fault["nsegs"]
    ndeps = dict_fault["ndeps"]
    doSubFault = dict_fault.get("doSubFault", True)

    flt = Fault(faultfile, nsegs, ndeps, doSubFault, origin=[])
    flt.load_fault()
    report(f"Fault loaded: {flt.nf} patches ({nsegs} x {ndeps}).")

    return {"data": data, "flt": flt}


def run_inversion(cfg: dict, data, flt, progress_callback=None) -> dict:
    """
    Run the full inversion pipeline (Step 2 of the workflow).

    Parameters
    ----------
    cfg : dict
        Full configuration dictionary.
    data : GeoData
        Loaded geodetic data (from run_load_data).
    flt : Fault
        Loaded fault geometry (from run_load_data).
    progress_callback : callable, optional
        Called with ``progress_callback(message: str)`` to report progress.

    Returns
    -------
    dict : {"green": green, "lap": lap, "sol": sol, "output": output}
    """

    def report(msg):
        logger.info(msg)
        if progress_callback:
            progress_callback(msg)

    dict_green = cfg["dict_green"]
    dict_weight = cfg["dict_weight"]
    dict_bound = cfg["dict_bound"]
    dict_export = cfg["dict_export"]
    dict_data = cfg["dict_data"]

    lapmethod = str(dict_data.get("lapmethod", "1"))

    # ---- Green's Function ----
    report("Building Green's functions...")
    grnmethod = dict_green.get("grnmethod", "okada")
    if grnmethod == "wang":
        from pyginvc.Greens.EDCMP import EDCMP
        green = EDCMP(flt, data, dict_green)
    else:
        green = Okada(flt, data, dict_green)
    green.build_greens()
    report("Green's functions built successfully.")

    # ---- Laplacian ----
    report("Building Laplacian smoothing matrix...")
    bcs = dict_green.get("bcs", [0, 0, 0, 0])
    lap = RectPlaneLap(flt.nsegs, flt.ndeps, bcs, method=lapmethod)
    lap.build_laps()
    report("Laplacian matrix built.")

    # ---- Inversion ----
    report("Running slip inversion...")
    sol = GeoInversion(flt, data, green, lap, dict_weight, dict_bound)
    sol.run_inversion_clean()

    n_smooth = len(sol.smo_facts)
    report(f"Inversion complete. {n_smooth} smoothing factor(s) evaluated.")

    # ---- Output ----
    outpath = dict_export.get("outpath", "./output")
    report(f"Saving results to {outpath}...")

    # Ensure output directory exists
    if outpath and not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True)

    out = Output(flt, data, green, sol=sol, archdir=outpath)
    out.OutputSolution()
    out.archive_outfile()
    report("Results saved successfully.")

    # Log summary
    if hasattr(sol, "misfit") and sol.misfit is not None:
        report(f"Final misfit: {sol.misfit[-1]:.4f}")
    if hasattr(sol, "moment") and sol.moment is not None:
        Mo, Mw = sol.moment[-1]
        report(f"Geodetic moment: M0 = {Mo:.3e} N·m, Mw = {Mw:.2f}")

    return {"green": green, "lap": lap, "sol": sol, "output": out}
