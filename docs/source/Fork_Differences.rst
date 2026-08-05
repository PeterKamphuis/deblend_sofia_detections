.. _fork-differences:

Changes from upstream
=====================

Scope and provenance
--------------------

This repository is an enhanced fork of
`PeterKamphuis/deblend_sofia_detections <https://github.com/PeterKamphuis/deblend_sofia_detections>`_.
The comparison on this page is intentionally pinned to exact committed states:

* upstream ``main`` at release ``v0.0.4``, commit
  `a6daef3 <https://github.com/PeterKamphuis/deblend_sofia_detections/commit/a6daef3>`_;
* this fork at release ``v1.0.0``, commit
  `d78fa1b <https://github.com/3rico/deblend_sofia_detections/commit/d78fa1b>`_.

Only committed changes through ``v1.0.0`` are described. Work that has not yet
been committed and tagged is excluded so that every statement can be reproduced
from the named revisions.

The fork retains the upstream watershed method, SoFiA input-product conventions,
counterpart-based filtering, master-mask update, and final field
re-parameterisation. Its changes concentrate on controlling which sources run,
surviving source-level failures, handling real-world optical FITS products, and
making the proposed splits easier to audit scientifically.

Summary of differences
----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 36 44

   * - Area
     - Upstream ``v0.0.4``
     - This fork ``v1.0.0``
   * - Source selection
     - Iterates over every SoFiA catalogue ID.
     - Adds a validated ``input.source_ids`` allowlist; an empty list still
       processes every source.
   * - Source exceptions
     - An uncaught exception stops the field run.
     - Continues by default, records detailed failures, and offers an explicit
       fail-fast setting.
   * - Debug QA
     - Primarily intermediate FITS masks and terminal output.
     - Adds optical/H I/catalogue and trial-child PNG overlays plus an ECSV record
       of plotted catalogue coordinates.
   * - Optical FITS dimensions
     - Expects a directly usable 2-D manual image.
     - Normalises RGB and other multi-plane images to 2-D using celestial WCS.
   * - Optical conversion
     - No standalone converter.
     - Adds ``scripts/convert_optical_fits_to_2d.py`` with safe output and several
       collapse methods.
   * - Watershed markers
     - Combines automatically detected optical sources with manual markers.
     - Optionally uses only manual catalogue markers while retaining automatic
       detections as QA context.
   * - Documentation and tests
     - Concise installation and configuration notes.
     - Adds an operational guide and regression tests for fork-specific behavior.

1. Select only the source IDs that need testing
------------------------------------------------

The fork adds ``input.source_ids``. For example::

  input:
    source_ids:
      - "42"
      - "57"
      - "221"

The requested IDs are normalised for comparison with Astropy/OmegaConf scalar
types and checked against the original SoFiA catalogue before optical images or
cubelets are processed. Unknown IDs raise an input error. Catalogue ordering and
catalogue scalar types are preserved, and duplicate requested values are processed
once. With ``source_ids: []`` the behavior remains the upstream all-source loop.

This option limits candidate processing, not the final scientific data model. If
one selected source is split successfully, new component labels are inserted into
the complete copied field mask and SoFiA re-parameterises that complete mask.
Unselected detections remain in the final field products.

2. Continue after one source fails
----------------------------------

The main catalogue loop now isolates exceptions per source. The default is::

  general:
    continue_on_source_error: true

An exception records the SoFiA ID, cubelet path, exception class, message, and full
traceback, then processing continues with the next selected ID. At the end of each
run, ``Watershed_Output/deblend_failures.log`` is overwritten with current-run
requested, successful, and failed counts and all recorded failures. It is also
rewritten after a completely successful run, preventing a stale failure report
from being mistaken for current output.

Set ``continue_on_source_error: false`` to stop at the first exception. The report
is written before that exception is re-raised.

3. Add optical, H I, catalogue, and child-component QA
------------------------------------------------------

With ``general.debug: true``, each processed source can receive these new products:

``optical_hi_catalogue_overlay_source_<ID>.png``
  A WCS-aware grayscale optical cutout with the original parent H I footprint in
  purple, faint moment-0 contours, cyan automatic optical-source outlines and
  centroids, and yellow manual or accepted NED catalogue positions.

``catalogue_positions_source_<ID>.ecsv``
  The catalogue origin, object name, RA, and Dec associated with each plotted
  catalogue position. An empty but valid table is written when no position falls
  inside the cutout.

``optical_hi_components_overlay_source_<ID>.png``
  The same optical and parent-source context with a distinct colour for each child
  returned by the first trial SoFiA parameterisation, plus its measured catalogue
  position. These are proposed children before counterpart merging or rejection;
  they are not automatically accepted sources.

The first QA image is written early and updated after counterpart searches, so it
often survives a later watershed or trial-SoFiA failure. Visualisation errors are
reported as warnings and do not abort the scientific source processing. The fork
adds Matplotlib as a runtime dependency for these products.

4. Handle multi-dimensional manual optical FITS files
-----------------------------------------------------

The fork no longer assumes that every manual optical FITS image is already 2-D.
It inspects WCS axis metadata, maps the two celestial pixel axes onto NumPy array
axes, and averages all non-celestial axes. This supports conventional plane-first
cubes as well as RGB/RGBA arrays whose colour axis is stored elsewhere.

The reduction records the original shape and collapsed axes, retains celestial
WCS for later cutouts, and raises focused input errors when the FITS file has no
usable image HDU, no identifiable celestial WCS, fewer than two dimensions, or no
overlap with the SoFiA field.

The standalone ``scripts/convert_optical_fits_to_2d.py`` utility offers:

* automatic ``<input>_2d.fits`` or explicit output naming;
* mean, median, or first-plane reduction of non-celestial axes;
* image-HDU selection by index or extension name;
* explicit overwrite control while forbidding input self-overwrite;
* float32 output by default;
* retained celestial WCS and useful observing metadata;
* FITS ``HISTORY`` provenance and post-write 2-D verification.

5. Allow a manual catalogue to control watershed markers
---------------------------------------------------------

The fork adds::

  input:
    manual_markers_only: true

In this mode, automatic Photutils source regions are removed from the marker array
before manual catalogue ellipses are added. Only manual markers overlapping the
parent H I mask can seed the optical watershed. Automatically detected regions
remain visible in cyan in the QA plot so the operator can see what was excluded,
and the plot title records the marker mode.

This mode requires at least one manual input table. It is intended for crowded
fields where a vetted counterpart list is safer than every optical segmentation
island. The default is ``false``, preserving combined automatic and manual marker
behavior.

6. Expand operations documentation and regression coverage
-----------------------------------------------------------

The fork replaces the short README with a beginner-oriented guide covering the
cube-mask-catalogue relationship, watershed routes, required SoFiA products,
installation, copied-workspace safety, source allowlists, manual images and
catalogues, outputs, scientific QA, and troubleshooting. The advanced
configuration reference documents each new setting.

Focused tests cover source-ID selection and error reporting, stale failure-log
replacement, WCS-aware QA plots and catalogue records, trial-child overlays,
multi-axis optical reduction, the standalone converter CLI, and manual-only marker
selection.

Compatibility and unchanged limitations
----------------------------------------

The new controls are compatible with the upstream workflow when left at their
defaults:

* ``input.source_ids: []`` processes all catalogue IDs;
* ``input.manual_markers_only: false`` retains automatic optical markers;
* ``general.debug: false`` does not generate the new PNG/ECSV QA products.

The exception policy is intentionally different: this fork continues after a
source-level error by default. Set ``general.continue_on_source_error: false`` for
upstream-style fail-fast behavior.

Scientific limitations remain. The software still divides only voxels already in
the parent SoFiA mask, assigns each voxel wholly to one child, and can mistake
rotating-disc peaks, clumps, projected optical objects, or displaced gas for
separate galaxies. It can overwrite the copied master mask and regenerate enabled
SoFiA products. The new controls and QA plots make decisions more manageable and
auditable; they do not remove the need for scientific review or the requirement to
work on a complete copy.

Commit-by-commit map
--------------------

* ``638f517`` — source-ID filtering, validation, initial full user guide, and
  source-selection tests.
* ``8bf9ae4`` — source-level exception isolation, current-run failure reports, and
  fail-fast control.
* ``89eaafb`` — optical/H I/catalogue QA overlay, catalogue-position ECSV, and
  non-fatal visualisation handling.
* ``c7c7111`` — multi-plane optical FITS normalisation and standalone conversion
  utility.
* ``11124b6`` — trial H I child-component QA overlay.
* ``d78fa1b`` — manual-catalogue-only watershed marker mode.
