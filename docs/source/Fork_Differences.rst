.. _fork-differences:

Project lineage and fork additions
==================================

Acknowledgement and purpose
---------------------------

This repository builds on
`PeterKamphuis/deblend_sofia_detections <https://github.com/PeterKamphuis/deblend_sofia_detections>`_.
We sincerely thank Peter Kamphuis for creating and maintaining the original
package, publishing it as open-source software, and providing the scientific and
technical foundation on which this fork depends. The watershed method also follows
the work of Qifeng Huang and collaborators, who retain credit for the underlying
method and implementation lineage.

This page records provenance and version-specific additions so users can tell
which settings and outputs belong to the code they are running. The description is
descriptive rather than evaluative: the additions reflect the operational needs of
a particular large-field workflow and may not be appropriate for every use case.

For scientific use, :doc:`Citing` gives separate citations for this fork, Peter
Kamphuis's original package, and the Huang et al. method paper, together with
suggested wording that explicitly thanks Peter for openly sharing the original
project.

Version scope
-------------

The notes are intentionally pinned to exact committed states:

* Peter Kamphuis's original project at release ``v0.0.4``, commit
  `a6daef3 <https://github.com/PeterKamphuis/deblend_sofia_detections/commit/a6daef3>`_;
* this fork's last tagged release, ``v1.0.0``, commit
  `d78fa1b <https://github.com/3rico/deblend_sofia_detections/commit/d78fa1b>`_;
* the current documented post-release implementation, commit
  `44532f9 <https://github.com/3rico/deblend_sofia_detections/commit/44532f9>`_.

Only committed additions through ``44532f9`` are described. The automatic DR10,
targeted optical-region deblending, and Gaia-diagnostic additions are newer than
the ``v1.0.0`` tag. Until they are included in a later release, scientific users
should record and cite the full Git commit rather than identifying those runs as
unmodified ``v1.0.0``.

The fork retains the original project's watershed method, SoFiA input-product
conventions, counterpart-based filtering, master-mask update, and final field
re-parameterisation. Its additions concentrate on controlling which sources run,
recording source-level failures, handling multi-plane optical FITS products, and
making proposed splits easier to audit scientifically. Post-``v1.0.0`` additions
also provide opt-in Legacy Surveys DR10 positional counterparts, optional H I
support filtering for those markers, targeted optical-region deblending, and
explicit Gaia-mask provenance products.

Summary of additions maintained in this fork
--------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Area
     - Addition and intended use in this fork
   * - Source selection
     - Adds a validated ``input.source_ids`` allowlist; an empty list still
       retains the original all-source workflow.
   * - Source exceptions
     - Records detailed source-level failures, continues by default, and offers
       an explicit fail-fast setting.
   * - Debug QA
     - Adds optical/H I/catalogue and trial-child PNG overlays plus an ECSV record
       of plotted catalogue coordinates.
   * - Optical FITS dimensions
     - Extends manual-image handling to normalise RGB and other multi-plane data
       to 2-D using celestial WCS.
   * - Optical conversion
     - Adds ``scripts/convert_optical_fits_to_2d.py`` with safe output and several
       collapse methods.
   * - Watershed markers
     - Optionally uses only manual catalogue markers while retaining automatic
       detections as QA context.
   * - Automatic catalogue
     - Optionally downloads and caches a field-limited Legacy Surveys DR10 Tractor
       table when no manual catalogue is supplied.
   * - DR10 marker filtering
     - Optionally retains a selected DR10 marker only when the same optical region
       contains a finite positive, beam-scale parent moment-0 maximum.
   * - Targeted optical deblending
     - Optionally applies Photutils multi-threshold deblending only to an optical
       region containing at least two mapped moment-0 peaks.
   * - Gaia-mask provenance
     - Writes the masked optical background and matching binary mask per source in
       debug mode, including query-success and masked-pixel metadata.
   * - Documentation and tests
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
once. With ``source_ids: []`` the original all-source loop is retained.

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

Set ``continue_on_source_error: false`` when stopping at the first exception is
preferred. The report is written before that exception is re-raised.

3. Add optical, H I, catalogue, and child-component QA
------------------------------------------------------

With ``general.debug: true``, each processed source can receive these new products:

``optical_hi_catalogue_overlay_source_<ID>.png``
  A WCS-aware grayscale optical cutout with the original parent H I footprint in
  purple, faint moment-0 contours, cyan automatic optical-source outlines and
  centroids, yellow manual/selected DR10/accepted NED catalogue positions, and red
  crosses for DR10 positions rejected by the optional moment-0 filter.

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

When the corresponding modes are active, the same source-level debug directory
also contains Gaia masking products and a DR10 moment-0 decision audit, described
below.

4. Handle multi-dimensional manual optical FITS files
-----------------------------------------------------

In addition to the original 2-D manual-image path, the fork accepts
multi-dimensional optical FITS images. It inspects WCS axis metadata, maps the two
celestial pixel axes onto NumPy array axes, and averages all non-celestial axes.
This supports conventional plane-first cubes as well as RGB/RGBA arrays whose
colour axis is stored elsewhere.

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

6. Use an automatic Legacy Surveys DR10 catalogue
--------------------------------------------------

Commit ``7c83d62`` adds these opt-in settings::

  input:
    manual_input_tables: [null]
    auto_query_catalogue: true
    galaxy_types: [REX, EXP, DEV, SER]

The software constructs a field footprint from the full-field moment-0 celestial
WCS, queries ``ls_dr10.tractor`` through the public NOIRLab Data Lab TAP service,
and caches the returned CSV together with the exact query metadata. Cache reuse
requires the table columns, sky bounds, selected morphology types, service URL,
and query text to match. ``input.original_tables: true`` forces a fresh download.

For each parent H I detection, a catalogue position must lie inside both the
projected parent mask and a positive optical-segmentation label. At most one
allowed row is retained per optical region, chosen by the highest finite
``flux_g``. Catalogue coordinates are assigned to the nearest segmentation pixel,
including at region boundaries. The selected position becomes a watershed marker
and may serve as a positional child counterpart.

Any supplied manual catalogue takes precedence and prevents the automatic query.
DR10 matches are positional, not spectroscopic. Neither a Tractor morphology nor
the brightest ``flux_g`` in a region proves that an object lies at the H I redshift
or owns the detected gas.

7. Optionally require parent moment-0 support for DR10 markers
--------------------------------------------------------------

The post-release setting::

  input:
    filter_dr10_markers_by_moment0_peaks: true

searches the parent moment-0 map for finite positive local maxima using an
elliptical footprint derived from ``BMAJ``, ``BMIN``, ``BPA``, and the celestial
pixel-scale matrix. A selected DR10 row remains a marker only when a peak maps
into that row's exact optical-segmentation label. A rejected catalogue-associated
label is removed from seeding, while unrelated automatic optical regions remain.

The filter is independent of ``input.use_peak_deblending`` and never applies to a
manual catalogue. It has no global amplitude or signal-to-noise threshold: a very
faint positive local maximum can pass. It is therefore a topology-based support
test, not astrophysical confirmation.

With ``general.debug: true``, the source directory records:

* ``parent_moment0_used_for_dr10_peak_filter_source_<ID>.fits``;
* ``moment0_peak_map_source_<ID>.fits``; and
* ``dr10_moment0_peak_filter_audit_source_<ID>.ecsv``.

The audit includes catalogue identity, morphology, ``flux_g``, optical label,
accepted/rejected status, matched peak position and value, and rejection reason.

8. Deblend only optical regions containing multiple H I peaks
--------------------------------------------------------------

The opt-in settings::

  input:
    deblend_optical_regions_with_multiple_moment0_peaks: true
    optical_deblend_nlevels: 32
    optical_deblend_contrast: 0.001
    optical_deblend_min_pixels: 20

project the beam-scale moment-0 peak map into the raw cyan optical
segmentation. Photutils multi-threshold deblending is applied only to an
original cyan label containing at least two peaks; the pixel membership of all
other labels is preserved. DR10 type filtering and highest-``flux_g`` selection
then run independently in each resulting sublabel, followed by the moment-0
acceptance filter. This limits the operation to the peak-selected labels rather
than applying it to the complete segmentation.

The option requires ``filter_dr10_markers_by_moment0_peaks: true``. It is
inactive for manual catalogues and when automatic DR10 or optical processing is
inactive. Debug mode writes the exact cyan segmentations before and after the
targeted operation as
``optical_segmentation_before_targeted_deblend_source_<ID>.fits`` and
``optical_segmentation_after_targeted_deblend_source_<ID>.fits``.

Multiple moment-0 peaks are only a trigger for a candidate optical segmentation.
They can arise within one rotating or clumpy galaxy, from noise or residual
foreground structure, or from unsuitable beam metadata. The resulting sublabels
must be inspected in the H I cube, channel maps, spectra, position-velocity
structure, and counterpart data before they are interpreted astrophysically.

9. Record per-source Gaia masking provenance
---------------------------------------------

In debug mode, commit ``7c83d62`` writes:

* ``background_gaia_masked_source_<ID>.fits``, the optical background immediately
  after applying the Gaia mask and before background subtraction or smoothing;
* ``gaia_star_mask_source_<ID>.fits``, a binary mask in which one denotes a masked
  pixel.

Both headers record ``GAIA_OK``, ``MASKNPIX``, and ``MASKFRAC``. This distinguishes
a successful query that found no maskable stars from a failed Gaia query followed
by the package's unmasked fallback. When debug mode is enabled and a cached
cleaned image lacks these products, the Gaia step is rebuilt so the diagnostics
correspond to the source being processed.

10. Expand operations documentation and regression coverage
------------------------------------------------------------

The fork complements the original project documentation with a beginner-oriented
guide covering the cube-mask-catalogue relationship, watershed routes, required
SoFiA products, installation, copied-workspace safety, source allowlists, manual
images and catalogues, outputs, scientific QA, and troubleshooting. The advanced
configuration reference documents each new setting.

Focused tests cover source-ID selection and error reporting, stale failure-log
replacement, WCS-aware QA plots and catalogue records, trial-child overlays,
multi-axis optical reduction, the standalone converter CLI, manual-only marker
selection, DR10 query caching and counterpart selection, WCS-aware moment-0
filtering, nearest-pixel catalogue assignment, debug audit products, and Gaia-mask
provenance. They also cover targeted optical splitting, preservation of already
separated regions, peak-map reuse, parameter validation, and the before/after
segmentation products.

Continuity with the original workflow and scientific limitations
-----------------------------------------------------------------

The original project's broad workflow remains available when the fork-specific
controls are left at their defaults:

* ``input.source_ids: []`` processes all catalogue IDs;
* ``input.manual_markers_only: false`` retains automatic optical markers;
* ``input.auto_query_catalogue: false`` makes no DR10 network request;
* ``input.filter_dr10_markers_by_moment0_peaks: false`` leaves selected DR10
  markers unfiltered when automatic catalogue mode is enabled;
* ``input.deblend_optical_regions_with_multiple_moment0_peaks: false`` leaves the
  initial connected optical segmentation unchanged;
* ``general.debug: false`` does not generate the new PNG/ECSV QA products.

The exception policy is a version-specific choice: this fork continues after a
source-level error by default. Set ``general.continue_on_source_error: false``
when fail-fast behavior is preferred.

Scientific limitations remain. The software still divides only voxels already in
the parent SoFiA mask, assigns each voxel wholly to one child, and can mistake
rotating-disc peaks, clumps, projected optical objects, or displaced gas for
separate galaxies. It can overwrite the copied master mask and regenerate enabled
SoFiA products. The new controls and QA plots make decisions more manageable and
auditable; they do not remove the need for scientific review or the requirement to
work on a complete copy.

.. warning::

   This fork is beta research software. Internal acceptance, a positional DR10
   association, moment-0 peak support, or conservation of parent flux does not
   establish that a split is astrophysically correct. An astronomer must visually
   inspect every changed source in the H I cube, channel maps, moment maps,
   spectra, position-velocity structure, and counterpart data. Ambiguous results
   should be rejected or classified as unresolved.

Fork addition history
---------------------

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
* ``7c83d62`` — automatic DR10 catalogue selection, optional beam-scale moment-0
  support filtering, DR10 decision audits, per-source Gaia-mask diagnostics, and
  nearest-pixel catalogue-to-segmentation assignment.
* ``44532f9`` — opt-in Photutils deblending restricted to optical regions with
  multiple mapped moment-0 peaks, peak-map reuse, before/after segmentation
  diagnostics, parameter validation, and focused regression tests.
