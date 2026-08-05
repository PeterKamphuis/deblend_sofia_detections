Welcome to 'deblend-sofia-detections' documentation!
=====================================================
**deblend-sofia-detections** is a Python package implementing the watershed
methods published in https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract.

This started as an adapted copy of the accompanying ipython notebook of that paper and evolved from there
and as such Qifeng Huang deserves credit for the original implementation.

.. warning::

   This is beta research software. Its output masks and catalogues are candidate
   deblending solutions, not confirmation that multiple astrophysical sources are
   present. An astronomer must visually inspect every proposed split in the H I
   cube, channel maps, moment maps, spectra, position-velocity structure, and
   counterpart data before using it scientifically. Passing the program's internal
   checks, or conserving the parent flux, does not establish that emission was
   assigned to the correct galaxy. Rejecting a split or classifying it as unresolved
   is appropriate when the available evidence is ambiguous.

Project lineage and acknowledgement
-----------------------------------

This documentation accompanies a fork of
`PeterKamphuis/deblend_sofia_detections <https://github.com/PeterKamphuis/deblend_sofia_detections>`_.
We sincerely thank Peter Kamphuis for creating, maintaining, and openly sharing
the original package on which this work is based. The core watershed workflow,
SoFiA integration, and project structure originate in that project and in the
method developed by Qifeng Huang and collaborators.

For version clarity, the fork-specific notes are pinned from Peter's ``v0.0.4``
release (``a6daef3``), through this fork's ``v1.0.0`` release (``d78fa1b``), to
the ``v1.1.0`` release and current documented development commit ``b215ca0``.
They cover
selected-source runs, source-level failure reporting, QA overlays,
multi-dimensional optical FITS support, a conversion utility, manual-only marker
control, automatic Legacy Surveys DR10 positional counterparts, optional
moment-0 support filtering, targeted optical-region deblending, Gaia-mask
diagnostics, position-velocity QA projections, and supporting tests. These notes
identify which repository version owns each option and are intended as provenance,
not as an evaluation of relative merit.

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
