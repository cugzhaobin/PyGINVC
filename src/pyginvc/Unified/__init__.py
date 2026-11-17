#!/usr/bin/env python
"""Self-contained output and plotting helpers for unified slip inversion.

These modules intentionally live outside ``pyginvc.Export``.  The existing
export classes assume a single flat fault type; the unified model can mix
rectangular and triangular faults, which requires a different HDF5 layout
and a viewer that understands it.  Old result files remain readable through
the legacy branch of :class:`pyginvc.Unified.UnifiedView`.

Import the concrete classes from their modules, e.g.
``from pyginvc.Unified.UnifiedOutput import UnifiedOutput``.  Keeping this
package initialiser empty avoids importing matplotlib for headless runs.
"""
