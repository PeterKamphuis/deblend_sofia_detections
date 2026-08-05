Welcome to 'deblend-sofia-detections' documentation!
=================================
**deblend-sofia-detections** is a Python package implementing the watershed
methods published in https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract.

This started as an adapted copy of the accompanying ipython notebook of that paper and evolved from there
and as such Qifeng Huang deserves credit for the original implementation.

Project lineage and acknowledgement
-----------------------------------

This documentation accompanies a fork of
`PeterKamphuis/deblend_sofia_detections <https://github.com/PeterKamphuis/deblend_sofia_detections>`_.
We sincerely thank Peter Kamphuis for creating, maintaining, and openly sharing
the original package on which this work is based. The core watershed workflow,
SoFiA integration, and project structure originate in that project and in the
method developed by Qifeng Huang and collaborators.

For version clarity, the fork-specific notes are pinned between Peter's ``v0.0.4``
release (``a6daef3``) and this fork's ``v1.0.0`` release (``d78fa1b``). They cover
selected-source runs, source-level failure reporting, QA overlays,
multi-dimensional optical FITS support, a conversion utility, manual-only marker
control, and supporting tests. These notes identify which repository version owns
each option and are intended as provenance, not as an evaluation of relative merit.

See :doc:`Fork_Differences` for detailed attribution, additions, continuity, and
compatibility notes. See :doc:`Citing` for paper-ready references, BibTeX, and
suggested acknowledgement wording.

Contents
--------
.. toctree::



  Project lineage and fork additions <Fork_Differences.rst>
  Introduction, Installation & Easy use <readme.md>
  The yaml File Input <Advanced.rst>
  Citing and acknowledging the software <Citing.rst>
  Copyright <license.rst>
