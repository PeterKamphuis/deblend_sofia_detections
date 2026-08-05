.. _citing:

Citing and acknowledging the software
======================================

Thank you for crediting the people whose work made this software possible. This
fork has two distinct parts of its software provenance: the additions maintained
by Eric Maina, and the original ``deblend_sofia_detections`` package created and
openly shared by Peter Kamphuis. The scientific watershed workflow follows Huang
et al. (2025). Keeping these contributions separate makes the credit accurate
without implying that Peter authored later fork-specific changes.

Scientific-use caution
----------------------

Citation does not replace validation. This package is beta research software, and
its output is a proposed segmentation rather than proof that the components are
distinct astrophysical sources. An astronomer must visually examine every changed
source in the H I cube, channel maps, moment maps, spectra, position-velocity
structure, and counterpart data before including the result in a scientific
analysis. A split passing the program's internal checks, or conserving the parent
flux, can still assign emission to the wrong galaxy. Ambiguous results should be
reported as rejected or unresolved rather than as secure deblends.

Recommended citations for a scientific paper
--------------------------------------------

If this fork contributes to a published analysis, please include all three of the
following references:

#. **The exact release of this fork that was used.** Cite the version or commit so
   readers can reproduce the software environment. For release ``v1.0.0``:

   Maina, E. (2026). *deblend-sofia-detections* (Version 1.0.0) [Computer
   software]. https://github.com/3rico/deblend_sofia_detections

   The automatic DR10 counterpart path, optional moment-0 support filter,
   targeted optical-region deblending, and per-source Gaia-mask diagnostics were
   introduced after ``v1.0.0``. The documented implementation containing all of
   them is commit ``44532f95bab9d7259a4f081666139907853a7f4b``. Until a later
   release contains them, cite that full commit hash rather than identifying such
   an analysis as unmodified ``v1.0.0``.

#. **Peter Kamphuis's original package.** This acknowledges the software
   foundation from which the fork developed:

   Kamphuis, P. (2026). *deblend_sofia_detections* (Version 0.0.4) [Computer
   software]. https://github.com/PeterKamphuis/deblend_sofia_detections

#. **The scientific method paper.** Cite this when the watershed-deblending
   workflow supports the analysis:

   Huang, Q., Wang, J., Lin, X., et al. (2025). WALLABY Pilot Survey: Star
   Formation Enhancement and Suppression in Gas-rich Galaxy Pairs. *The
   Astrophysical Journal, 980*, 157.
   https://doi.org/10.3847/1538-4357/ad9579

Please adapt punctuation and author-list truncation to the journal's reference
style, while retaining the software title, author, year, version, and repository
URL. If an archival DOI is assigned to a later software release, prefer that DOI
to the repository URL.

Suggested wording
-----------------

A methods section could say:

   H I detections were deblended with ``deblend-sofia-detections`` v1.0.0
   (Maina 2026), a fork of Peter Kamphuis's ``deblend_sofia_detections`` v0.0.4
   (Kamphuis 2026), using the watershed-deblending approach described by Huang et
   al. (2025).

Replace ``v1.0.0`` with the full Git commit when post-release functionality was
used.

An acknowledgements section may additionally say:

   We thank Peter Kamphuis for creating and openly sharing the original
   ``deblend_sofia_detections`` package on which this fork is based.

These sentences are templates, not mandatory wording. They are intended to make
the lineage explicit while preserving each contributor's role.

BibTeX
------

The following entries can be copied into a bibliography database and adjusted to
the target journal's conventions:

.. code-block:: bibtex

   @software{maina_deblend_sofia_detections_2026,
     author  = {Maina, Eric},
     title   = {deblend-sofia-detections},
     version = {1.0.0},
     year    = {2026},
     url     = {https://github.com/3rico/deblend_sofia_detections},
     note    = {Computer software}
   }

   @software{maina_deblend_sofia_detections_44532f9_2026,
     author  = {Maina, Eric},
     title   = {deblend-sofia-detections},
     version = {Git commit 44532f95bab9d7259a4f081666139907853a7f4b},
     year    = {2026},
     url     = {https://github.com/3rico/deblend_sofia_detections/tree/44532f95bab9d7259a4f081666139907853a7f4b},
     note    = {Computer software; post-v1.0.0 DR10, targeted optical-deblending, and Gaia-diagnostic implementation}
   }

   @software{kamphuis_deblend_sofia_detections_2026,
     author  = {Kamphuis, Peter},
     title   = {deblend_sofia_detections},
     version = {0.0.4},
     year    = {2026},
     url     = {https://github.com/PeterKamphuis/deblend_sofia_detections},
     note    = {Computer software}
   }

   @article{huang_wallaby_pairs_2025,
     author  = {Huang, Qifeng and Wang, Jing and Lin, Xuchen and Oh, Se-Heon and
                Chen, Xinkai and Catinella, Barbara and Deg, Nathan and
                Dénes, Helga and For, Bi-Qing and Koribalski, Baerbel and
                Lee-Waddell, Karen and Rhee, Jonghwan and Shen, Austin and
                Shao, Li and Spekkens, Kristine and Staveley-Smith, Lister and
                Westmeier, Tobias and Wong, O. Ivy and Bosma, Albert},
     title   = {{WALLABY Pilot Survey: Star Formation Enhancement and
                 Suppression in Gas-rich Galaxy Pairs}},
     journal = {The Astrophysical Journal},
     volume  = {980},
     pages   = {157},
     year    = {2025},
     doi     = {10.3847/1538-4357/ad9579}
   }

Machine-readable citation metadata
----------------------------------

The repository root contains ``CITATION.cff``. GitHub and compatible reference
managers can use this file to generate a citation for the fork. Its ``references``
also record Peter Kamphuis's original project and the Huang et al. method paper.

For an unreleased checkout, replace the release version with the output of
``deblend -v`` when it identifies the build, and record the full Git commit hash.
The release or commit should also appear in the paper's reproducibility records.
