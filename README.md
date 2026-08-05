# deblend-sofia-detections

## Project lineage and additions in this fork

This repository builds on
[`PeterKamphuis/deblend_sofia_detections`](https://github.com/PeterKamphuis/deblend_sofia_detections),
created and maintained by Peter Kamphuis. We are grateful to Peter for developing
the original package, making it openly available, and providing the implementation
and project structure on which this fork depends. The underlying watershed method
also follows the work of Qifeng Huang and collaborators, credited below.

The purpose of this section is attribution and version clarity. It records which
options and support tools were added in this fork so users can choose the correct
documentation for the code they are running. These notes are descriptive rather
than evaluative: the additions reflect the needs of a particular large-field
deblending workflow and are not assumed to be preferable for every use case.

### Version scope

The notes are pinned to reproducible Git revisions:

| Project state | Revision | Role in this documentation |
| --- | --- | --- |
| Peter Kamphuis's original project | [`v0.0.4` / `a6daef3`](https://github.com/PeterKamphuis/deblend_sofia_detections/commit/a6daef3) | Foundation and reference revision from which this fork developed. |
| This fork's previous tagged release | [`v1.0.0` / `d78fa1b`](https://github.com/3rico/deblend_sofia_detections/commit/d78fa1b) | Tagged baseline containing selected-source controls, failure reporting, optical QA, multi-plane FITS support, and manual-only markers. |
| This fork's current release | [`v1.1.0`](https://github.com/3rico/deblend_sofia_detections/tree/v1.1.0) | Adds automatic DR10 counterparts, optional moment-0 support filtering, targeted optical-region deblending, Gaia-mask diagnostics, and corrected fork-maintainer metadata. |

Only committed additions included in `v1.1.0` are described. Scientific analyses
should record and cite this release tag, together with the configuration and
scientific inputs used for the run.

### Additions maintained in this fork

| Area | Addition and intended use |
| --- | --- |
| Source selection | `input.source_ids` provides an optional validated allowlist for controlled trials on known candidate blends. |
| Failure handling | Source-level exceptions can be recorded while later sources continue, with an option to retain fail-fast behavior. |
| Optical/catalogue QA | Per-source PNG overlays and ECSV records connect optical detections, catalogue positions, the parent H I footprint, moment-0 structure, and trial child components. |
| Manual optical FITS input | RGB and other multi-plane FITS data can be reduced to 2-D by identifying celestial axes from WCS and collapsing the remaining axes. |
| Standalone optical conversion | `scripts/convert_optical_fits_to_2d.py` provides mean, median, or first-plane collapse modes with guarded output handling. |
| Marker control | `input.manual_markers_only` allows a vetted catalogue to define watershed seeds while automatic detections remain available as QA context. |
| Automatic catalogue | `input.auto_query_catalogue` can download and cache a field-limited Legacy Surveys DR10 Tractor table when no manual catalogue is supplied. |
| DR10 marker filtering | `input.filter_dr10_markers_by_moment0_peaks` can require a selected DR10 object's exact optical region to contain a positive, beam-scale parent moment-0 maximum. |
| Targeted optical deblending | Photutils multi-threshold deblending can be limited to cyan regions containing two or more moment-0 peaks, while other regions retain their existing pixel membership. |
| Gaia-mask provenance | Debug mode writes the optical cutout immediately after Gaia masking and the matching binary mask, including query-success and masked-pixel metadata. |
| Documentation and tests | The fork includes an operational guide and focused regression tests for its version-specific controls and QA paths. |

### Detailed addition notes

1. **Targeted source-ID runs.** Set `input.source_ids` to quoted SoFiA catalogue
   IDs such as `["42", "57"]`. The selection is validated before optical-image or
   cubelet processing, follows catalogue order, and processes duplicate requested
   IDs once. An empty list preserves the original all-source loop. If a selected
   split is accepted, the final re-parameterisation still operates on the complete
   field mask so unchanged detections remain in the field catalogue.

2. **Source-level failure isolation.** The main catalogue loop catches an exception
   from one cubelet and continues when `general.continue_on_source_error: true`,
   which is the fork default. Each run rewrites
   `Watershed_Output/deblend_failures.log` with requested, successful, and failed
   counts plus the source ID, cubelet path, exception type, reason, and traceback.
   Set the option to `false` when fail-fast behavior is preferred.

3. **Two complementary QA overlays.** With `general.debug: true`,
   `optical_hi_catalogue_overlay_source_<ID>.png` shows the optical background,
   original H I footprint, moment-0 contours, automatic optical detections, and
   catalogue positions. `catalogue_positions_source_<ID>.ecsv` records the plotted
   catalogue coordinates. After trial SoFiA parameterisation,
   `optical_hi_components_overlay_source_<ID>.png` adds separately coloured child
   masks and their measured positions. Plotting failures are warnings and do not
   abort scientific processing.

4. **Robust multi-dimensional optical FITS handling.** Manual images are inspected
   using their WCS. The code maps the celestial pixel axes to NumPy axes, averages
   all other axes to produce a 2-D image, retains celestial WCS information, and
   raises focused errors for missing image data, unusable celestial WCS, or absent
   sky overlap. The standalone converter supports explicit HDU selection,
   `mean`/`median`/`first` collapse methods, overwrite protection, float32 output,
   provenance in FITS `HISTORY`, and post-write dimensionality verification.

5. **Manual-catalogue-only watershed markers.** When
   `input.manual_markers_only: true`, automatic Photutils regions do not seed the
   optical watershed. Only manual catalogue markers that overlap the parent H I
   mask are used. Automatic regions remain visible in cyan in the QA plot, and the
   plot title records the active marker mode. The option requires a manual input
   table and defaults to `false` for compatibility.

6. **Beginner and operations documentation.** The fork explains the cube-mask-
   catalogue relationship, the optical and peak watershed routes, the exact SoFiA
   products required, copied-workspace safety, selected-source and whole-field
   workflows, parameter meanings, output interpretation, scientific validation,
   and troubleshooting. Focused unit tests cover each fork-specific helper and QA
   path.

7. **Automatic Legacy Surveys DR10 counterparts.** When no manual catalogue is
   supplied, `input.auto_query_catalogue: true` downloads a field-limited subset
   of `ls_dr10.tractor` from NOIRLab Data Lab and caches both the CSV and exact
   query metadata. Eligible morphologies are controlled by `input.galaxy_types`.
   Per parent source, at most one eligible row is selected in each optical region,
   using the highest finite `flux_g`. These are positional counterparts, not
   spectroscopic confirmations.

8. **Optional H I support for DR10 markers.** With
   `input.filter_dr10_markers_by_moment0_peaks: true`, a selected DR10 marker is
   retained only when a finite positive local maximum, measured with a
   synthesized-beam footprint, falls inside the same optical-segmentation label.
   The filter uses no global flux or signal-to-noise threshold. Debug mode records
   the exact moment-0 input, binary peak map, and an ECSV decision audit; rejected
   catalogue positions remain visible as red crosses.

9. **Targeted splitting of blended optical regions.** With
   `input.deblend_optical_regions_with_multiple_moment0_peaks: true`, Photutils
   multi-threshold deblending is applied only to a raw cyan region containing at
   least two beam-scale moment-0 peaks. Other regions retain their existing pixel
   membership. The DR10 per-region selection and moment-0 filter then operate on
   the new sublabels.

10. **Per-source Gaia-mask diagnostics.** Debug mode writes the background cutout
   immediately after Gaia masking and a corresponding binary mask under the
   source's `debug_products` directory. FITS headers record whether the Gaia query
   succeeded, the number of masked pixels, and their fraction of the cutout. A
   failed Gaia query is therefore distinguishable from a successful zero-star
   result.

### Foundations retained from the original project

The core scientific intent is unchanged: the program starts from an existing
SoFiA detection, proposes integer child labels inside the parent mask, asks SoFiA
to measure those labels, checks counterparts, updates the master mask when a split
survives, and finally re-parameterises the full field. The fork does not make
deblending scientifically automatic, does not discover unmasked H I emission, and
does not make the workflow non-destructive. The copied-workspace warning below
applies equally to the original project and this fork.

For the original broad all-source workflow on the new controls, leave
`input.source_ids` empty, `input.manual_markers_only` false,
`input.auto_query_catalogue` false, and
`input.filter_dr10_markers_by_moment0_peaks` false, and
`input.deblend_optical_regions_with_multiple_moment0_peaks` false. Set
`general.continue_on_source_error: false` if the run must stop at the first
source-level exception. Debug overlays are only created when `general.debug` is
enabled.

The corresponding documentation page is
[Project lineage and fork additions](https://github.com/3rico/deblend_sofia_detections/blob/need_based_modifications/docs/source/Fork_Differences.rst).

## Package overview

`deblend-sofia-detections` is a post-processing tool for H I source detections
created by [SoFiA-2](https://gitlab.com/SoFiA-Admin/SoFiA-2). It is designed for
the case where SoFiA has labelled emission from two or more galaxies as one
source. The package proposes a division of that existing source mask, asks SoFiA
to measure the proposed components, and uses optical/spectroscopic counterparts
to decide whether more than one component should remain.

The package can process a full survey or cluster-field catalogue containing many
normal and blended detections. By default it checks every SoFiA catalogue ID. Use
the optional `source_ids` setting to run a controlled trial on known candidates.

> [!WARNING]
> **Beta research software—astronomical review is required.** The generated child
> masks and catalogues are candidate deblending solutions, not confirmed
> astrophysical sources. An internal “accepted” split only means that it passed the
> software's procedural checks. Before using a split scientifically, an astronomer
> must visually inspect the masks against the H I cube channel by channel, moment
> maps, spectra, position-velocity structure, and credible optical or spectroscopic
> counterparts. Flux conservation alone is not validation: the total flux can be
> preserved while emission is assigned to the wrong galaxy. Ambiguous cases should
> be classified as rejected or unresolved rather than reported as secure splits.

> [!CAUTION]
> The program can overwrite the master SoFiA mask and then rerun SoFiA, which can
> regenerate the catalogue, cubelets, moments, and other enabled products. Work on
> a complete copy of the original SoFiA run, never on your only production copy.

## Contents

- [What problem does it solve?](#what-problem-does-it-solve)
- [Terminology](#terminology)
- [What the program does](#what-the-program-does)
- [Requirements](#requirements)
- [Required SoFiA input](#required-input-from-the-original-sofia-run)
- [Safe working copy](#make-a-safe-working-copy)
- [Quick start](#quick-start-selected-known-blends)
- [Configuration](#important-configuration-settings)
- [Manual optical images and catalogues](#supplying-a-manual-optical-image)
- [Automatic DR10 catalogue](#automatically-querying-the-legacy-surveys-dr10-catalogue)
- [Outputs](#understanding-the-output)
- [Scientific quality assurance](#scientific-quality-assurance)
- [Troubleshooting](#troubleshooting)

## What problem does it solve?

An H I data cube has two spatial axes and one spectral/velocity axis. SoFiA finds
significant voxels and stores their source membership in an integer mask. A simple
field might look like this:

```text
data cube voxel values  +  source mask labels  ->  catalogue measurements
```

If two galaxies overlap on the sky and in velocity, SoFiA may connect their
emission and assign all affected voxels to one ID:

```text
connected emission from Galaxy A and Galaxy B -> source ID 42
```

The catalogue then contains one row whose flux, position, spectrum, and velocity
width describe the combined detection. This package attempts to replace that one
mask region with multiple labelled regions:

```text
part of the original mask -> source ID 42
remaining component       -> a new source ID above the old catalogue maximum
```

The measured cube values are not altered. Only voxel ownership in the mask is
changed. Each voxel is assigned wholly to one component; the program cannot assign
fractions of one voxel to several galaxies.

The method is based on the workflow described by
[Huang et al. (2025, ApJ, 980, 157)](https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract).

## Terminology

| Term | Meaning in this package |
| --- | --- |
| **Cube** | The original 3-D H I measurement: two sky axes plus a spectral/velocity axis. |
| **Voxel** | One 3-D cube element, equivalent to a pixel in one velocity channel. |
| **Mask** | A 3-D integer array aligned with the cube. Zero is background; each positive integer identifies one SoFiA source. |
| **Catalogue** | The SoFiA table containing one row of measured properties for each mask ID. |
| **Cubelet** | A smaller cube, mask, and set of moments cut around one catalogue source. |
| **Marker** | A seed location representing a possible child galaxy before watershed growth. |
| **Watershed** | A segmentation algorithm that grows labelled regions from markers until they meet. |
| **Counterpart** | An optical or spectroscopic galaxy associated with a proposed H I component. |

The parent detection is the original SoFiA source. The children are the proposed
sources created by dividing that parent's mask.

## What the program does

For each selected SoFiA source, the current pipeline performs the following steps:

1. It opens the source cubelet, source mask, and moment-0 map from the existing
   SoFiA run.
2. It obtains a field-wide optical image. The default is DSS2 Red from SkyView;
   a WCS-aware optical FITS image can be supplied instead.
3. For each source, it makes an optical cutout, estimates and subtracts the
   background, attempts to mask Gaia DR3 stars, smooths the image, and detects
   optical objects inside the H I footprint.
4. It uses the optical objects as markers for a watershed segmentation of either
   the full H I cube (default) or the moment-0 map.
5. If peak deblending is enabled, it finds significant local H I maxima and runs
   a second 3-D watershed. The number of peak markers is capped by
   `maximum_no_peaks`.
6. It chooses a candidate mask in this order: successful peak-based 3-D mask,
   successful optical 3-D mask, then successful optical 2-D mask.
7. It calls the separate `sofia` executable with the S+C finder, linking,
   reliability filtering, and dilation disabled. The bundled trial template also
   has threshold finding disabled, so SoFiA measures the proposed labelled
   components rather than rediscovering sources.
8. It searches NED, any supplied manual counterpart tables, and any selected
   automatic DR10 positions. Proposed components without a suitable counterpart
   may be merged back together.
9. If more than one viable component remains, it inserts the relabelled source mask
   into the master field mask.
10. After all selected IDs have been checked, it reruns SoFiA on the complete field
    mask if at least one split was accepted.

```mermaid
flowchart LR
    A["Existing SoFiA cube, mask, catalogue, and cubelets"] --> B["One selected source ID"]
    B --> C["Optical markers and H I peak markers"]
    C --> D["2-D or 3-D watershed masks"]
    D --> E["SoFiA re-parameterises proposed children"]
    E --> F["NED and manual-counterpart checks"]
    F --> G["Accepted child labels enter the master mask"]
    G --> H["SoFiA re-parameterises the complete field"]
```

This is an automated proposal and filtering pipeline, not a definitive
astrophysical blend classifier. Every accepted split still needs scientific
inspection.

## What it does not do

- It does not run initial H I source finding. You must complete a normal SoFiA run
  first.
- It does not discover emission outside the existing parent source mask.
- It does not guarantee that multiple H I peaks correspond to multiple galaxies.
- It does not guarantee that an optical object owns the nearest H I emission.
- It does not fractionally divide unresolved voxels.
- It does not remove the need to inspect channel maps, spectra, moment maps, and
  counterparts.

A rotating disc can have several H I peaks or a double-horned spectrum. Disturbed
gas, tidal bridges, and ram-pressure-stripped tails can also be offset from their
stellar counterparts. In such cases, “unresolved” or “reject the split” can be the
correct scientific conclusion.

## Requirements

### Python package

- Python 3.11 or later
- The Python dependencies declared in `pyproject.toml`

Install the released package in a virtual environment:

```bash
python3.11 -m venv deblend_venv
source deblend_venv/bin/activate
python -m pip install --upgrade pip
python -m pip install deblend_sofia_detections
```

To use the fork-specific additions, install this repository rather than the PyPI
release of the original project:

```bash
git clone https://github.com/3rico/deblend_sofia_detections.git
cd deblend_sofia_detections
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

An editable install is the right choice when testing changes from this checkout,
including configuration features that have not yet appeared in a packaged release.

Confirm that the command is installed:

```bash
deblend -v
```

### SoFiA-2 executable

SoFiA is a separate program. Installing this Python package does **not** install
SoFiA. The executable must be available when candidate components and the final
field are parameterised.

```bash
command -v sofia
```

If the executable has another name or is not on `PATH`, set it explicitly:

```yaml
internal:
  sofia: /absolute/path/to/sofia
```

On a cluster, SoFiA may first need to be loaded through the local module system.

### Network access

The default workflow accesses:

- SkyView for a DSS2 Red image, unless a manual optical FITS image is supplied;
- Gaia DR3 for foreground-star masking; and
- NED for counterpart matching; and
- NOIRLab Data Lab when `input.auto_query_catalogue: true` and no manual
  catalogue is supplied.

Gaia and NED failures have limited fallbacks, but reliable counterpart information
is central to deciding whether a proposed split survives. For crowded fields,
supplying a vetted manual spectroscopic catalogue is strongly recommended.

## Required input from the original SoFiA run

Use the **exact SoFiA parameter file used for the full field whose detections you
want to deblend**. Do not use the generic `sofia_template.par` as the input for an
existing run.

The deblender reads these values from your original `.par` file:

```text
input.data
output.directory
output.filename
```

Those settings must resolve to the original cube and the corresponding SoFiA
products. If `output.filename = cluster`, the expected layout is approximately:

```text
working_copy/
├── cluster_sofia.par
├── cluster_cube.fits
└── sofia_output/
    ├── cluster_cat.xml          # cluster_cat.txt is also supported
    ├── cluster_mask.fits        # full-field labelled mask
    ├── cluster_mom0.fits        # full-field moment-0 map
    └── cluster_cubelets/
        ├── cluster_1_cube.fits
        ├── cluster_1_mask.fits
        ├── cluster_1_mom0.fits
        ├── cluster_2_cube.fits
        ├── cluster_2_mask.fits
        ├── cluster_2_mom0.fits
        └── ...
```

The original SoFiA run therefore needs catalogue, mask, moment, and cubelet output
enabled. Standard parameterised SoFiA catalogues should contain the fields used by
the deblender, including source ID, sky position, velocity, flux, source bounds,
ellipse measurements, and noise/width measurements.

### Which `.par` file should I use?

| Parameter file | Use it? | Reason |
| --- | --- | --- |
| The `.par` used to detect the complete cluster or survey field | Yes | It points to the master cube and all existing SoFiA products. |
| A `.par` from a deliberately isolated one-source SoFiA run | Only for that isolated test | The deblender will only know about products from that run. |
| The generated `sofia_template.par` | No | It is an internal template for measuring proposed child masks. |
| The generated `deblend_sofia.par` | No | The program creates it when re-parameterising the updated field mask. |

## Make a safe working copy

Copy the cube, original `.par` file, catalogue, mask, moment maps, and cubelet
directory into a trial workspace. Preserve an untouched baseline separately.

```bash
cp -a /data/cluster_sofia_run /data/cluster_deblend_trial
cd /data/cluster_deblend_trial
```

Inspect the copied `.par` file:

```bash
rg -n '^\s*(input\.data|output\.directory|output\.filename)\s*=' cluster_sofia.par
```

Update absolute paths that still point at the production output. A safe copied file
might contain:

```text
input.data = /data/cluster_deblend_trial/cluster_cube.fits
output.directory = /data/cluster_deblend_trial/sofia_output
output.filename = cluster
output.overwrite = true
```

The data cube may remain in a shared read-only location, but `output.directory`
must point to the copied SoFiA products. Relative output paths are resolved from
the input cube directory by the current code, so absolute paths are safest. The
final rerun preserves the original SoFiA output settings; in the working copy,
`output.overwrite` must permit those products to be regenerated.

For an additional before/after comparison:

```bash
mkdir -p baseline
cp sofia_output/cluster_cat.xml baseline/cluster_cat_before.xml
cp sofia_output/cluster_mask.fits baseline/cluster_mask_before.fits
```

## Quick start: selected known blends

The safest first run is a small allowlist of known candidate blends.

Create `cluster_deblend.yml`:

```yaml
input:
  # This is the original full-field SoFiA parameter file in the working copy.
  sofia_parameters: /data/cluster_deblend_trial/cluster_sofia.par

  # Quote IDs. Remove this block or use [] to process the entire catalogue.
  source_ids:
    - "42"
    - "57"
    - "221"

  # Used only when manual_input_tables contains no catalogue path.
  auto_query_catalogue: false
  galaxy_types: [REX, EXP, DEV, SER]
  filter_dr10_markers_by_moment0_peaks: false
  deblend_optical_regions_with_multiple_moment0_peaks: false
  optical_deblend_nlevels: 32
  optical_deblend_contrast: 0.001
  optical_deblend_min_pixels: 20

  use_optical_deblending: true
  use_cube_deblending: true
  use_peak_deblending: true
  maximum_no_peaks: 5
  compactness: 0.0

general:
  verbose: true
  ncpu: 1
  optical_pixel_scale: 5.0
  counterpart_region: Ellipse
  debug: true

directories:
  run_directory: /data/cluster_deblend_trial
  ancillary_directory: /data/cluster_deblend_trial/ancillary_data
  watershed_directory: /data/cluster_deblend_trial/Watershed_Output
```

Run it and keep a log:

```bash
source deblend_venv/bin/activate
cd /data/cluster_deblend_trial
deblend configuration_file=/data/cluster_deblend_trial/cluster_deblend.yml \
  2>&1 | tee cluster_deblend.log
```

Configuration and command-line values use OmegaConf `key=value` syntax, not
traditional `--key value` options. A setting can be overridden for one run:

```bash
deblend configuration_file=cluster_deblend.yml \
  general.debug=true \
  input.maximum_no_peaks=3 \
  input.compactness=0.0
```

The override does not edit the YAML file. Record the exact command used for each
scientific trial.

## Running the complete catalogue

An omitted or empty `source_ids` list preserves the original behavior and attempts
every catalogue detection:

```yaml
input:
  sofia_parameters: /data/cluster_deblend_trial/cluster_sofia.par
  source_ids: []
```

Normal detections are allowed in the same field as blends. The code checks them
one at a time and only updates the master mask when a candidate remains split after
the measurement and counterpart stages. This automated screening can still make
false positive or false negative decisions, so review every changed source.

For a large cluster field, begin with a short allowlist. Once the configuration,
optical image, counterpart table, and QA procedure are validated, expand the list
or remove it.

## Generate the packaged examples

Run the following command in a disposable or configuration directory:

```bash
deblend print_examples=true
```

It writes:

- `deblend_sofia_detections_default.yml`, containing all user-facing settings in
  a recommended starter configuration;
- `sofia_template.par`, the generic internal SoFiA measurement template.

The command exits after writing the files. The template is useful for inspection,
but it is not a replacement for the `.par` file from your original SoFiA run.

The generated starter intentionally enables the automatic DR10, moment-0 support,
and targeted optical-region path, while disabling the separate H I peak watershed:

```yaml
input:
  auto_query_catalogue: true
  filter_dr10_markers_by_moment0_peaks: true
  deblend_optical_regions_with_multiple_moment0_peaks: true
  use_peak_deblending: false
```

These are generated-example choices, not changes to the compatibility defaults
listed below. The user must still replace `input.sofia_parameters` and the
mandatory `directories.data_directory: ???` value before running.

## Important configuration settings

### Input settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `input.sofia_parameters` | `sofia_input.par` | Original `.par` file for the complete SoFiA run being deblended. |
| `input.source_ids` | `[]` | Optional catalogue-ID allowlist. Empty means all IDs. Unknown IDs stop before source processing. |
| `input.manual_input_tables` | `[null]` | Optional `.csv`, `.txt`, or compatible pickled counterpart tables. |
| `input.auto_query_catalogue` | `false` | When `true` and no manual catalogue is supplied, download and cache the field subset of `ls_dr10.tractor` from the public NOIRLab Data Lab TAP service. |
| `input.galaxy_types` | `[REX, EXP, DEV, SER]` | DR10 Tractor morphology types accepted as galaxy counterparts. The comparison is case-insensitive. |
| `input.filter_dr10_markers_by_moment0_peaks` | `false` | When automatic DR10 selection is active, require a finite positive beam-scale moment-0 maximum inside the selected row's exact cyan optical region. Manual catalogues are never filtered. |
| `input.deblend_optical_regions_with_multiple_moment0_peaks` | `false` | Run Photutils multi-threshold deblending only on a raw cyan region containing at least two moment-0 peaks. Requires the DR10 moment-0 filter to be enabled. |
| `input.optical_deblend_nlevels` | `32` | Number of Photutils multi-threshold levels used for targeted optical deblending. |
| `input.optical_deblend_contrast` | `0.001` | Minimum fraction of a region's flux required for a local optical peak to become a separate sublabel. |
| `input.optical_deblend_min_pixels` | `20` | Minimum connected pixels above a Photutils threshold required for a targeted object to be deblended. |
| `input.manual_markers_only` | `false` | When `true`, only manual catalogue positions inside the parent H I mask seed optical watershed runs; automatic optical detections remain QA diagnostics. Requires a manual input table. |
| `input.manual_optical_image` | `[null]` | Optional celestial-WCS optical FITS image. A 2-D image is used directly; multi-plane images are collapsed over non-celestial axes. Only the first image is currently used. |
| `input.original_tables` | `false` | Reread supplied text tables instead of using cached pickle files, and force a fresh automatic DR10 download when that mode is active. This is not a backup/safety option. |
| `input.original_images` | `false` | Recreate cached optical products from the supplied image or SkyView. This is not a backup/safety option. |
| `input.use_optical_deblending` | `true` | Detect optical markers and try an optical watershed split. |
| `input.use_cube_deblending` | `true` | Use the 3-D cube for optical-marker watershed; `false` uses the 2-D moment-0 map. |
| `input.use_peak_deblending` | `true` | Try the H I local-peak 3-D watershed and prefer it when successful. |
| `input.maximum_no_peaks` | `5` | Maximum H I peak markers retained for one parent source. |
| `input.compactness` | `0.0` | Watershed compactness penalty. Start at zero for extended/disturbed H I. |

### General settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `general.verbose` | `true` | Print progress and decisions. Keep enabled during initial trials. |
| `general.ncpu` | available CPUs | Number passed to the peak-finding interface. The current fast local-maximum implementation itself is not parallelised. |
| `general.multiprocessing` | `true` | Reserved multiprocessing preference. The current main source loop does not branch on this setting. |
| `general.continue_on_source_error` | `true` | Record a source-level exception and continue with the next selected ID. |
| `general.optical_pixel_scale` | `5.0` | Requested number of optical pixels across the H I beam; downloaded pixels are capped at 4 arcsec. |
| `general.counterpart_region` | `Ellipse` | Spatial region used for counterpart matching. Supported choices include `Beam`, `3Beam`, `Box`, `Ellipse`, and `Full_Ellipse`. |
| `general.debug` | `false` | Write additional diagnostic products, including a per-source optical/H I/catalogue QA plot. **This is not a dry-run mode.** |

### Directory settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `directories.run_directory` | current directory | Directory in which processing runs. It must already exist. |
| `directories.ancillary_directory` | `ancillary_data` | Cached optical images, cleaned cutouts, tables, and debug material. |
| `directories.watershed_directory` | `Watershed_Output` | Per-source watershed and trial-SoFiA products. |

The `sofia` configuration section is populated from the original `.par` file and
normally should not be set manually. `internal.sofia` is the useful exception when
the SoFiA executable is not simply named `sofia`.

### Source-level failure handling

An exception raised while processing one cubelet is recorded and the run continues
with the next selected SoFiA ID by default. At the end of every run,
`Watershed_Output/deblend_failures.log` contains the source ID, cubelet path,
exception type, reason, and full traceback for every failed attempt. The file is
overwritten on each run, including successful runs, so an old failure report cannot
be mistaken for the current result.

A source that completes normally without being split is counted as successful.
Only an uncaught exception is recorded as a failure. To stop at the first exception,
set:

```yaml
general:
  continue_on_source_error: false
```

## Supplying a manual optical image

Set a FITS image with valid celestial WCS that overlaps the full SoFiA moment-0
field:

```yaml
input:
  manual_optical_image:
    - /data/cluster_deblend_trial/optical/r_band.fits
  original_images: true
```

The program cuts the image to the SoFiA field and stores processed products in the
ancillary directory. After the first successful run, set `original_images: false`
to reuse the cached image. Set it back to `true` after replacing the optical FITS
or when you intentionally want to rebuild all processed cutouts.

A manual image must contain a valid celestial WCS. Two-dimensional FITS images
are used directly. For an RGB, RGBA, or other multi-plane image, the program
identifies the two celestial axes from the WCS and averages every non-celestial
array axis to form a grayscale 2-D image before making the cutout. The log records
the original shape and collapsed axes. If the celestial axes cannot be identified,
the file contains no suitable image HDU, or it does not overlap the SoFiA field,
the run stops with an explicit input error rather than the less informative
Astropy dimension-mismatch traceback.

### Standalone optical FITS converter

The repository includes `scripts/convert_optical_fits_to_2d.py` for converting
an optical image before running the deblender. It accepts any FITS image with at
least two dimensions and a celestial WCS. The default operation averages all
non-celestial planes. Supplying only the input creates `<input>_2d.fits` in the
same directory:

```bash
python scripts/convert_optical_fits_to_2d.py \
  "/data/Abell-548-colour.fits"
```

This writes `/data/Abell-548-colour_2d.fits`. An explicit output path is also
supported:

```bash
python scripts/convert_optical_fits_to_2d.py \
  "/data/Abell-548-colour.fits" \
  "/data/Abell-548-optical-2d.fits"
```

For a more outlier-resistant collapse, use the median:

```bash
python scripts/convert_optical_fits_to_2d.py \
  input.fits output_2d.fits \
  --method median
```

Use `--method first` to select the first plane along every non-celestial axis,
`--hdu 1` (or an extension name) to choose a specific image HDU, and
`--overwrite` to replace an existing output file. The script never allows the
input path itself to be overwritten. It writes a float32 image by default,
retains the celestial WCS and useful observing metadata, records the conversion
in FITS `HISTORY`, and verifies that the written result is two-dimensional.

Point the YAML at the converted file and rebuild the cached optical products:

```yaml
input:
  manual_optical_image:
    - /data/Abell-548-optical-2d.fits
  original_images: true
```

## Supplying a manual counterpart table

Manual counterparts are especially valuable in crowded clusters where automatic
database matching may be incomplete or ambiguous. A CSV must provide a commented
column-name row followed by a commented unit row. For example:

```csv
# Name,RA,DEC,Velocity,PA,sma,smb
# -,deg,deg,km/s,deg,arcsec,arcsec
Galaxy_A,150.123400,-26.123400,14500,45,18,10
Galaxy_B,150.128100,-26.120700,14620,112,13,8
```

Configure it with:

```yaml
input:
  manual_input_tables:
    - /data/cluster_deblend_trial/catalogues/counterparts.csv
  # Optional: exclude automatic Photutils detections from watershed markers.
  manual_markers_only: true
  original_tables: true
```

With `manual_markers_only: true`, the manual table acts as the optical-marker
allowlist. Every catalogue entry whose marker overlaps the parent H I footprint
can seed the watershed; use a table containing only the counterparts intended
for the selected blends. Automatic Photutils detections are still shown as cyan
QA outlines and crosses, but they are not watershed seeds. Catalogue positions
outside the parent H I footprint can still appear as yellow QA dots and are not
used as markers. The QA image title reports the active watershed-marker mode.

The important fields are:

- `Name`: stable counterpart identifier;
- `RA`, `DEC`: sky position, normally in decimal degrees;
- `Velocity`: spectroscopic systemic velocity, normally in km/s;
- `PA`: optional optical position angle;
- `sma`, `smb`: optional semi-major and semi-minor marker sizes.

Accepted aliases exist for several velocity and size columns, but explicit names
and units are safest. Check carefully for redshift versus velocity, m/s versus
km/s, degrees versus sexagesimal coordinates, arcseconds versus arcminutes, and
semi-axis versus full diameter.

The loader caches text/CSV tables as `.pkl` files. After editing a source table,
use `original_tables: true` for the next run so the text is reread. It can then be
returned to `false` to reuse the cache.

## Automatically querying the Legacy Surveys DR10 catalogue

If no manual counterpart table is available, the deblender can query the public
[NOIRLab Data Lab TAP service](https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/CatalogDataAccessTAPSCS/CatalogDataAccessTAPSCS.html)
for the current field's Legacy Surveys DR10 Tractor sources:

```yaml
input:
  manual_input_tables: [null]
  auto_query_catalogue: true
  galaxy_types:
    - REX
    - EXP
    - DEV
    - SER
  filter_dr10_markers_by_moment0_peaks: true
  deblend_optical_regions_with_multiple_moment0_peaks: true
  optical_deblend_nlevels: 32
  optical_deblend_contrast: 0.001
  optical_deblend_min_pixels: 20
```

The query uses `ls_dr10.tractor`, keeps primary-brick rows, and downloads only the
identifier, position, morphology, and `flux_g` columns needed by the deblender.
The field subset and query metadata are cached as
`ancillary_data/catalogues/<SoFiA-basename>_ls_dr10_tractor.csv` and `.json`.
Set `input.original_tables: true` to force a fresh download; otherwise an exact
cache match is reused.

The field query is prepared once, before the per-source loop. The full-field
moment-0 FITS must therefore exist and contain usable celestial WCS. A missing
image, invalid response, network failure, or unavailable TAP service stops the run
before source processing unless a compatible cached query is available.

For each parent H I detection, a DR10 position must fall inside both the purple
H I footprint and a positive cyan optical-detection region. At most one DR10 row
is associated with each cyan region. Rows whose `type` is not in
`input.galaxy_types` are rejected, and if two or more eligible rows occupy the
same region, the row with the highest finite `flux_g` is selected. The selected
position replaces that cyan region's generic optical marker and is also
available as a positional (not spectroscopic) counterpart for the trial H I
child component.

When `input.filter_dr10_markers_by_moment0_peaks: true`, those per-region DR10
selections pass through one additional gate before the watershed. The parent
moment-0 map is searched for finite positive local maxima with an elliptical
maximum-filter footprint derived from `BMAJ`, `BMIN`, `BPA`, and the celestial
WCS. A selected row is retained only when a peak maps into the exact same cyan
segmentation label. No global peak-flux or signal-to-noise threshold is imposed.
A rejected row and its associated cyan region are removed from watershed
seeding; unrelated automatic cyan regions remain unchanged. This option is
independent of `input.use_peak_deblending`, has no effect when automatic DR10 or
optical deblending is inactive, and never filters a manual catalogue. Missing or
invalid parent moment-0 WCS or beam metadata raises a source-level input error.

When
`input.deblend_optical_regions_with_multiple_moment0_peaks: true`, the same
peak map is first projected into the raw cyan segmentation. Photutils
multi-threshold deblending runs only for a cyan label containing at least two
peaks; every other label is left unchanged. DR10 morphology and highest-`flux_g`
selection then run independently in each resulting sublabel, followed by the
moment-0 acceptance filter. This targeted trigger avoids applying Photutils
deblending globally, which can split structure within an already distinct galaxy.
The option requires `filter_dr10_markers_by_moment0_peaks: true`, and has no
effect for manual catalogues or when automatic DR10/optical deblending is inactive.
`optical_deblend_nlevels`, `optical_deblend_contrast`, and
`optical_deblend_min_pixels` control the Photutils operation. Lower contrast is
more aggressive; review the before/after segmentation and QA overlay after any
parameter change.

This trigger is a segmentation aid, not evidence that the peaks belong to
different galaxies. A rotating disc, star-forming clumps, residual foreground
structure, noise, or an imperfect beam model can produce multiple maxima. Treat
every resulting sublabel as a candidate and apply the full channel-map, spectral,
position-velocity, and counterpart review described in the beta warning above.

Any non-null entry in `input.manual_input_tables` takes absolute precedence. In
that case no DR10 request is made, even if `auto_query_catalogue: true`; an
invalid manual path therefore raises an input error rather than silently falling
back to the online catalogue.

> [!WARNING]
> A selected DR10 row is only a positional prior. Tractor morphology, high
> `flux_g`, proximity to H I, or a positive moment-0 maximum does not demonstrate
> that the object is a galaxy at the H I redshift or that it owns the associated
> gas. The optional moment-0 filter deliberately has no global S/N threshold, so
> even a very faint positive maximum can pass. Inspect every retained and rejected
> marker against redshift information, the optical image, channel maps, spectra,
> and position-velocity structure.

## Understanding the output

With default directories, source ID 42 produces a directory similar to:

```text
sofia_output/
└── Watershed_Output/
    └── watershed_source_42/
        ├── deblended_cube_mask_based_on_optical.fits
        ├── deblended_moment0_mask_based_on_optical.fits
        ├── deblended_cube_mask_based_on_peaks.fits
        ├── final_mask.fits
        ├── utilized_mask.fits
        ├── Sofia_Output/
        │   ├── sofia_input.par
        │   ├── ..._cat.xml
        │   ├── ..._mask.fits
        │   └── ..._cubelets/
        └── debug_products/
            ├── background_gaia_masked_source_42.fits
            ├── gaia_star_mask_source_42.fits
            ├── parent_moment0_used_for_dr10_peak_filter_source_42.fits
            ├── moment0_peak_map_source_42.fits
            ├── dr10_moment0_peak_filter_audit_source_42.ecsv
            ├── optical_segmentation_before_targeted_deblend_source_42.fits
            ├── optical_segmentation_after_targeted_deblend_source_42.fits
            ├── optical_hi_catalogue_overlay_source_42.png
            ├── optical_hi_components_overlay_source_42.png
            └── catalogue_positions_source_42.ecsv
```

Not every file appears for every source:

- `deblended_*_based_on_optical.fits` files are candidate optical-marker masks;
- `deblended_cube_mask_based_on_peaks.fits` is the H I peak candidate;
- `final_mask.fits` can be written while unmatched components are merged;
- `utilized_mask.fits` is the relabelled accepted mask inserted into the master
  field mask;
- `Sofia_Output/` contains the temporary SoFiA measurement of proposed children;
- `debug_products/` exists when `general.debug: true` and contains intermediate
  markers, smoothed cubes, watershed diagnostics, and per-source QA images.
  `background_gaia_masked_source_<ID>.fits` is the optical cutout immediately
  after Gaia masking, with masked pixels stored as NaN.
  `gaia_star_mask_source_<ID>.fits` is the matching binary mask (`1` for masked,
  `0` for retained) and records `GAIA_OK`, `MASKNPIX`, and `MASKFRAC` in its FITS
  header. When the DR10 moment-0 peak filter is active,
  `parent_moment0_used_for_dr10_peak_filter_source_<ID>.fits` preserves the
  exact in-memory parent map used for the decision,
  `moment0_peak_map_source_<ID>.fits` stores the binary representative peak
  pixels, and `dr10_moment0_peak_filter_audit_source_<ID>.ecsv` records every
  selected DR10 identity, morphology, `flux_g`, optical label, decision, matched
  peak coordinate/value, and rejection reason. In the overlays, accepted DR10
  positions remain yellow circles and rejected diagnostic-only positions are
  red crosses. When targeted optical deblending is enabled, the
  `optical_segmentation_before_targeted_deblend_source_<ID>.fits` and
  `optical_segmentation_after_targeted_deblend_source_<ID>.fits` files preserve
  the exact cyan label maps on either side of Photutils processing. Their headers
  record the targeted original-label count and deblending parameters. In
  `optical_hi_catalogue_overlay_source_<ID>.png`, the optical cutout is the
  grayscale background, the original H I footprint is purple, faint lavender
  contours show moment-0 intensity, cyan outlines and crosses mark automatically
  detected optical sources, yellow dots mark accepted catalogue positions, and
  red crosses mark DR10 positions rejected by the optional moment-0 filter.
  `catalogue_positions_source_<ID>.ecsv` records the catalogue name, object name,
  RA, and Dec for every plotted catalogue position; detailed filter decisions
  remain in the separate DR10 audit table. If no catalogue position falls inside
  the cutout, the table is empty and the plot says so explicitly.
- `optical_hi_components_overlay_source_<ID>.png` retains the same optical,
  parent-H I, optical-detection, and catalogue context, then adds one distinct
  contour colour for every child from the first trial SoFiA parameterisation.
  A matching coloured circle and label mark the child's SoFiA catalogue RA/Dec.
  These are raw candidate children shown before counterpart-based merging or
  rejection; they are not confirmation of separate galaxies or accepted final
  catalogue sources.

The QA plot is first written immediately after optical-source detection and is
updated after counterpart matching. Consequently, an early plot is normally
still available if a later watershed or trial-SoFiA step fails. Plotting errors
are reported as warnings and do not fail or skip the scientific source
processing.

If at least one split is accepted, the code also:

1. overwrites the copied full-field `<basename>_mask.fits` with the new labels;
2. writes `deblend_sofia.par` next to the original parameter file; and
3. reruns SoFiA on the complete updated mask using the original output settings.

The final catalogue and other enabled products therefore remain field-level SoFiA
outputs, including all unchanged detections and any new child IDs.

## Scientific quality assurance

Finding `utilized_mask.fits` means a split survived the software's internal path;
it does not prove that the split is physically correct.

For every changed source, inspect the original parent and proposed children in a
3-D viewer such as CARTA or a combination of DS9 and spectral/PV tools. At minimum:

1. Display mask IDs using discrete colours or integer contours.
2. Overlay the masks on the optical image and H I moment-0 map.
3. Step through individual velocity channels.
4. Inspect moment-1/velocity maps and a position-velocity diagram.
5. Compare each child spectrum with the original parent spectrum.
6. Check positions and systemic velocities against the proposed counterparts.
7. Repeat reasonable configurations, especially peak limits and the choice of
   peak versus optical-only deblending.

A plausible child should be spatially coherent over consecutive channels, span a
meaningful number of resolution elements, have a sensible spectrum and velocity,
and correspond to a credible galaxy. Warning signs include:

- a component present in only one channel;
- a marker located on a foreground star or unrelated projected object;
- a normal rotating disc split into two horns or bright clumps;
- a boundary that cuts abruptly across continuous velocity structure;
- a tidal or stripped tail changing ownership after a small parameter change;
- large disagreement between optical and peak-based masks; or
- a child with no credible positional and velocity counterpart.

Check approximate integrated-flux conservation:

```text
fractional difference = abs(sum(child fluxes) - parent flux) / abs(parent flux)
```

Good flux agreement is necessary but not sufficient. The total can be conserved
while emission is assigned to the wrong child.

A useful review classification is:

- **Accept**: coherent and well-supported split;
- **Accept with caveat**: useful division with a documented ambiguity;
- **Reject**: likely over-segmentation or incorrect assignment;
- **Unresolved**: available resolution does not support a defensible division.

## Troubleshooting

### `deblend` is not found

Activate the virtual environment and confirm the installation:

```bash
source deblend_venv/bin/activate
python -m pip show deblend_sofia_detections
command -v deblend
```

### `sofia` is not found

SoFiA must be installed or loaded separately:

```bash
command -v sofia
```

If needed, set `internal.sofia` to the executable's absolute path. A failed SoFiA
subprocess writes `sofia_output.txt` in the relevant SoFiA run directory.

### The catalogue or cubelets cannot be found

Check the three path-defining settings in the copied `.par` file:

```bash
rg -n '^\s*(input\.data|output\.directory|output\.filename)\s*=' cluster_sofia.par
```

Then verify the expected products:

```bash
ls sofia_output/cluster_cat.xml
ls sofia_output/cluster_mask.fits
ls sofia_output/cluster_mom0.fits
ls sofia_output/cluster_cubelets/cluster_42_cube.fits
ls sofia_output/cluster_cubelets/cluster_42_mask.fits
ls sofia_output/cluster_cubelets/cluster_42_mom0.fits
```

If these do not exist, regenerate the original SoFiA run with catalogue, mask,
moments, and cubelets enabled.

### A requested source ID is rejected

`source_ids` must match IDs in the original SoFiA catalogue. IDs should be quoted
in YAML. The error lists requested IDs that are absent.

### A known blend is not split

Possible causes include:

- only one usable optical marker was detected;
- too few significant H I peaks survived the size/beam checks;
- proposed children were merged because counterparts were not found;
- the optical image is too shallow, poorly registered, or contaminated;
- the original SoFiA mask does not contain enough emission to separate the objects.

Retry with `general.debug: true`, inspect the marker and candidate-mask products,
and consider a manual optical image or counterpart table. Do not treat increasing
`maximum_no_peaks` as proof that additional components are real.

### Changes to a catalogue or optical image are ignored

Set the corresponding refresh option for one run:

```yaml
input:
  original_tables: true
  original_images: true
```

These settings refresh cached inputs. They do not make the run non-destructive and
do not protect the master SoFiA products.

### The full-field run is slow

Use `source_ids` to validate a small candidate list first. Disable debug products
after the workflow is understood, reuse cached optical/table products, and only
then expand to the complete catalogue.

## Citing and acknowledging the software

Thank you for giving scientific credit to the people whose work made this package
possible. If this fork contributes to a paper, please cite:

1. **The exact version of this fork used in the analysis.** For `v1.1.0`:
   Maina, E. (2026), *deblend-sofia-detections*, Version 1.1.0 [Computer
   software], https://github.com/3rico/deblend_sofia_detections.
2. **Peter Kamphuis's original package**, which provided the software foundation:
   Kamphuis, P. (2026), *deblend_sofia_detections*, Version 0.0.4 [Computer
   software], https://github.com/PeterKamphuis/deblend_sofia_detections.
3. **The scientific method paper** when the watershed-deblending workflow supports
   the analysis: Huang, Q., Wang, J., Lin, X., et al. (2025), “WALLABY Pilot
   Survey: Star Formation Enhancement and Suppression in Gas-rich Galaxy Pairs,”
   *The Astrophysical Journal*, 980, 157,
   https://doi.org/10.3847/1538-4357/ad9579.

Version `v1.1.0` includes the automatic DR10 counterpart path, optional moment-0
support filter, targeted optical-region deblending, and per-source Gaia-mask
diagnostics. Cite `v1.1.0` when any of these paths contributes to the analysis.

Suggested methods wording:

> H I detections were deblended with `deblend-sofia-detections` v1.1.0
> (Maina 2026), a fork of Peter Kamphuis's `deblend_sofia_detections` v0.0.4
> (Kamphuis 2026), using the watershed-deblending approach described by Huang et
> al. (2025).

Replace `v1.1.0` with the exact later release or full Git commit if another
version is used.

Suggested acknowledgement:

> We thank Peter Kamphuis for creating and openly sharing the original
> `deblend_sofia_detections` package on which this fork is based.

The repository includes machine-readable
[`CITATION.cff`](https://github.com/3rico/deblend_sofia_detections/blob/need_based_modifications/CITATION.cff)
metadata.
Full paper-ready references and BibTeX are provided in the
[citation guide](https://github.com/3rico/deblend_sofia_detections/blob/need_based_modifications/docs/source/Citing.rst).
Adapt the formatting to the target
journal, but retain the software author, title, version, year, and URL. When using
an unreleased checkout, also record the full Git commit hash.

## Reproducibility checklist

Before a production run, record:

- package version from `deblend -v`;
- SoFiA version and executable path;
- the original and generated SoFiA parameter files;
- the deblender YAML and exact command line;
- source IDs attempted;
- manual table and optical-image versions;
- automatic-catalogue cache metadata and DR10 filter settings, when used;
- the per-source Gaia and DR10 debug/audit products, when enabled;
- baseline catalogue and mask checksums;
- per-source QA classification and notes.

## Further documentation and development status

- [Advanced configuration reference](https://github.com/3rico/deblend_sofia_detections/blob/need_based_modifications/docs/source/Advanced.rst)
- [Citation and acknowledgement guide](https://github.com/3rico/deblend_sofia_detections/blob/need_based_modifications/docs/source/Citing.rst)
- [Original project by Peter Kamphuis](https://github.com/PeterKamphuis/deblend_sofia_detections)
- [SoFiA-2 project](https://gitlab.com/SoFiA-Admin/SoFiA-2)
- [Method paper](https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract)

The package is beta research software. The scientific-review warning near the top
of this README and the quality-assurance procedure above apply to every proposed
split. Its software and scientific provenance is recorded in the lineage and
citation sections above.

Bug reports and focused improvements are welcome through the
[fork repository](https://github.com/3rico/deblend_sofia_detections).
