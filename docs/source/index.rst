Welcome to 'deblend-sofia-detections' documentation!
=================================
**deblend-sofia-detections** is a python package implementation and enhanced version of 
the watershed methods published in https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract

This started as an adapted copy of the accompanying ipython notebook of that paper and evolved from there
and as such Qifeng Huang deserves credit for the original implementation.

Changes from upstream
---------------------

This documentation describes an enhanced fork of
`PeterKamphuis/deblend_sofia_detections <https://github.com/PeterKamphuis/deblend_sofia_detections>`_.
The documented comparison is pinned between upstream ``v0.0.4`` (``a6daef3``)
and this fork's ``v1.0.0`` (``d78fa1b``). The fork adds selected-source runs,
source-level failure isolation and reporting, optical/H I/catalogue QA overlays,
multi-dimensional optical FITS support and a conversion utility, manual-only
watershed marker control, expanded documentation, and regression tests. The core
watershed method and the requirement to work on a copy of the SoFiA products are
unchanged.

See :doc:`Fork_Differences` for the detailed behavior and compatibility notes.

Contents
--------
.. toctree::



  Changes from upstream <Fork_Differences.rst>
  Introduction, Installation & Easy use <readme.md>
  The yaml File Input <Advanced.rst>
  Copyright <license.rst>
