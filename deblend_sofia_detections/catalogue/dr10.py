"""Legacy Surveys DR10 catalogue download and optical-counterpart selection."""

from __future__ import annotations

from functools import lru_cache
import csv
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from urllib.parse import urlencode
from urllib.request import urlopen

import astropy.units as u
from astropy.io import fits
from astropy.table import QTable, Table
from astropy.wcs import WCS
import numpy as np
from scipy.ndimage import label as label_connected
from scipy.ndimage import maximum_filter

from deblend_sofia_detections.support.errors import DownloadError, InputError


DR10_TAP_URL = "https://datalab.noirlab.edu/tap/sync"
DR10_TABLE = "ls_dr10.tractor"
DR10_COLUMNS = (
    "release",
    "brickid",
    "brickname",
    "objid",
    "brick_primary",
    "ls_id",
    "type",
    "ra",
    "dec",
    "flux_g",
)
_CACHE_VERSION = 1
_TYPE_PATTERN = re.compile(r"^[A-Z0-9_]+$")


def has_manual_catalogue(cfg) -> bool:
    """Return whether at least one manual catalogue path was supplied."""
    tables = getattr(cfg.input, "manual_input_tables", None) or []
    return any(table not in (None, "") for table in tables)


def normalise_galaxy_types(galaxy_types) -> tuple[str, ...]:
    """Validate and normalise the configured DR10 morphology allowlist."""
    normalised = []
    for source_type in galaxy_types or []:
        value = str(source_type).strip().upper()
        if not value:
            continue
        if not _TYPE_PATTERN.fullmatch(value):
            raise InputError(
                "input.galaxy_types entries may contain only letters, numbers, "
                f"and underscores; received {source_type!r}."
            )
        if value not in normalised:
            normalised.append(value)
    if not normalised:
        raise InputError(
            "input.galaxy_types must contain at least one DR10 source type when "
            "input.auto_query_catalogue=true."
        )
    return tuple(normalised)


def _perimeter_pixels(width: int, height: int, samples: int = 64):
    """Return pixel coordinates sampled around an image perimeter."""
    horizontal = np.linspace(-0.5, width - 0.5, samples)
    vertical = np.linspace(-0.5, height - 0.5, samples)
    x_values = np.concatenate(
        (
            horizontal,
            horizontal,
            np.full(samples, -0.5),
            np.full(samples, width - 0.5),
        )
    )
    y_values = np.concatenate(
        (
            np.full(samples, -0.5),
            np.full(samples, height - 0.5),
            vertical,
            vertical,
        )
    )
    return x_values, y_values


def _smallest_ra_interval(ra_values):
    """Return the smallest circular RA interval containing all values."""
    values = np.sort(np.mod(np.asarray(ra_values, dtype=float), 360.0))
    if values.size == 0:
        raise InputError("Could not determine an RA range from the image WCS.")
    gaps = np.diff(np.concatenate((values, [values[0] + 360.0])))
    largest_gap = int(np.argmax(gaps))
    start = float(values[(largest_gap + 1) % values.size])
    end = float(values[largest_gap])
    return start, end, start > end


def sky_bounds_from_wcs(wcs, width: int, height: int) -> dict:
    """Calculate an RA/Dec bounding box from a celestial image WCS."""
    if width <= 0 or height <= 0:
        raise InputError("Cannot query DR10 for an image with an empty shape.")
    celestial_wcs = wcs.celestial
    x_values, y_values = _perimeter_pixels(width, height)
    ra_values, dec_values = celestial_wcs.pixel_to_world_values(x_values, y_values)
    finite = np.isfinite(ra_values) & np.isfinite(dec_values)
    if not np.any(finite):
        raise InputError("The image WCS produced no finite sky coordinates.")
    ra_start, ra_end, ra_wraps = _smallest_ra_interval(ra_values[finite])
    return {
        "ra_start": ra_start,
        "ra_end": ra_end,
        "ra_wraps": bool(ra_wraps),
        "dec_min": float(np.min(np.asarray(dec_values)[finite])),
        "dec_max": float(np.max(np.asarray(dec_values)[finite])),
    }


def sky_bounds_from_header(header) -> dict:
    """Calculate catalogue-query bounds from a FITS image header."""
    width = int(header.get("NAXIS1", 0))
    height = int(header.get("NAXIS2", 0))
    wcs = WCS(header)
    if not wcs.has_celestial:
        raise InputError("The field moment-0 FITS header has no celestial WCS.")
    return sky_bounds_from_wcs(wcs, width, height)


def build_dr10_query(bounds: dict, galaxy_types) -> str:
    """Build the anonymous Data Lab TAP query for the current field."""
    source_types = normalise_galaxy_types(galaxy_types)
    type_values = ", ".join(f"'{source_type}'" for source_type in source_types)
    if bounds["ra_wraps"]:
        ra_clause = (
            f"(ra >= {bounds['ra_start']:.12f} OR "
            f"ra <= {bounds['ra_end']:.12f})"
        )
    else:
        ra_clause = (
            f"ra BETWEEN {bounds['ra_start']:.12f} "
            f"AND {bounds['ra_end']:.12f}"
        )
    return (
        f"SELECT {', '.join(DR10_COLUMNS)} FROM {DR10_TABLE} "
        "WHERE brick_primary = 1 "
        f"AND type IN ({type_values}) "
        f"AND {ra_clause} "
        f"AND dec BETWEEN {bounds['dec_min']:.12f} "
        f"AND {bounds['dec_max']:.12f}"
    )


def _cache_paths(cfg) -> tuple[Path, Path]:
    cache_directory = Path(cfg.directories.ancillary_directory) / "catalogues"
    cache_directory.mkdir(parents=True, exist_ok=True)
    safe_basename = re.sub(r"[^A-Za-z0-9_.-]+", "_", cfg.sofia.basename)
    catalogue_path = cache_directory / f"{safe_basename}_ls_dr10_tractor.csv"
    metadata_path = catalogue_path.with_suffix(".json")
    return catalogue_path, metadata_path


def _expected_metadata(bounds, galaxy_types, query):
    return {
        "cache_version": _CACHE_VERSION,
        "tap_url": DR10_TAP_URL,
        "table": DR10_TABLE,
        "columns": list(DR10_COLUMNS),
        "bounds": bounds,
        "galaxy_types": list(normalise_galaxy_types(galaxy_types)),
        "query": query,
    }


def _csv_has_expected_columns(catalogue_path: Path) -> bool:
    try:
        with catalogue_path.open(newline="", encoding="utf-8-sig") as handle:
            columns = next(csv.reader(handle))
    except (OSError, StopIteration, UnicodeError):
        return False
    return set(DR10_COLUMNS).issubset(columns)


def _cache_is_current(catalogue_path, metadata_path, expected_metadata):
    if not catalogue_path.is_file() or not metadata_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    return metadata == expected_metadata and _csv_has_expected_columns(
        catalogue_path
    )


def _tap_url(query: str) -> str:
    parameters = urlencode(
        {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": query,
        }
    )
    return f"{DR10_TAP_URL}?{parameters}"


def _tap_error(response_start: bytes) -> str | None:
    text = response_start.decode("utf-8", errors="replace")
    if "<VOTABLE" not in text.upper() and "QUERY_STATUS" not in text.upper():
        return None
    match = re.search(
        r'<INFO[^>]*name=["\']QUERY_STATUS["\'][^>]*value=["\']ERROR["\'][^>]*>'
        r"(.*?)</INFO>",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if match:
        return re.sub(r"\s+", " ", match.group(1)).strip()
    return "The TAP service returned a VOTable response instead of CSV."


def download_dr10_catalogue(cfg, urlopen_func=urlopen) -> str:
    """Download or reuse the configured field's DR10 Tractor subset."""
    moment0_path = Path(cfg.sofia.directory) / f"{cfg.sofia.basename}_mom0.fits"
    if not moment0_path.is_file():
        raise InputError(
            "Cannot determine the DR10 query footprint because the field "
            f"moment-0 image does not exist: {moment0_path}"
        )
    bounds = sky_bounds_from_header(fits.getheader(moment0_path))
    galaxy_types = normalise_galaxy_types(cfg.input.galaxy_types)
    query = build_dr10_query(bounds, galaxy_types)
    expected_metadata = _expected_metadata(bounds, galaxy_types, query)
    catalogue_path, metadata_path = _cache_paths(cfg)

    if not cfg.input.original_tables and _cache_is_current(
        catalogue_path, metadata_path, expected_metadata
    ):
        if cfg.general.verbose:
            print(f"Reusing cached Legacy Surveys DR10 catalogue {catalogue_path}")
        return str(catalogue_path)

    if cfg.general.verbose:
        print(
            "Downloading the Legacy Surveys DR10 Tractor catalogue for the "
            f"field to {catalogue_path}"
        )

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{catalogue_path.name}.",
            suffix=".tmp",
            dir=catalogue_path.parent,
            delete=False,
        ) as temporary_file:
            temporary_path = Path(temporary_file.name)
            with urlopen_func(_tap_url(query), timeout=3600) as response:
                response_start = response.read(8192)
                error_message = _tap_error(response_start)
                if error_message is not None:
                    raise DownloadError(
                        f"The Legacy Surveys DR10 TAP query failed: {error_message}"
                    )
                temporary_file.write(response_start)
                shutil.copyfileobj(response, temporary_file)

        if not _csv_has_expected_columns(temporary_path):
            raise DownloadError(
                "The Legacy Surveys DR10 TAP response is not a CSV containing "
                f"the required columns: {', '.join(DR10_COLUMNS)}"
            )
        os.replace(temporary_path, catalogue_path)
        temporary_path = None

        with tempfile.NamedTemporaryFile(
            mode="w",
            prefix=f".{metadata_path.name}.",
            suffix=".tmp",
            dir=metadata_path.parent,
            encoding="utf-8",
            delete=False,
        ) as metadata_file:
            json.dump(expected_metadata, metadata_file, indent=2, sort_keys=True)
            metadata_file.write("\n")
            temporary_metadata_path = Path(metadata_file.name)
        os.replace(temporary_metadata_path, metadata_path)
    except DownloadError:
        raise
    except Exception as error:
        raise DownloadError(
            "Could not download the Legacy Surveys DR10 catalogue from "
            f"{DR10_TAP_URL}: {error}"
        ) from error
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    if cfg.general.verbose:
        print(f"Successfully downloaded the DR10 catalogue to {catalogue_path}")
    return str(catalogue_path)


def prepare_dr10_catalogue(cfg):
    """Resolve manual precedence and prepare an automatic DR10 catalogue."""
    cfg.internal.auto_catalogue_path = None
    if has_manual_catalogue(cfg):
        if cfg.input.auto_query_catalogue and cfg.general.verbose:
            print(
                "A manual input catalogue was supplied; it takes precedence, "
                "so the automatic DR10 query will not run."
            )
        return None
    if not cfg.input.auto_query_catalogue:
        return None
    cfg.internal.auto_catalogue_path = download_dr10_catalogue(cfg)
    return cfg.internal.auto_catalogue_path


@lru_cache(maxsize=2)
def load_dr10_catalogue(catalogue_path: str) -> QTable:
    """Load a downloaded DR10 CSV, caching it for per-source selection."""
    try:
        table = QTable(Table.read(catalogue_path, format="ascii.csv"))
    except Exception as error:
        raise InputError(
            f"Could not read the downloaded DR10 catalogue {catalogue_path}: {error}"
        ) from error
    missing = [column for column in DR10_COLUMNS if column not in table.colnames]
    if missing:
        raise InputError(
            f"Downloaded DR10 catalogue {catalogue_path} is missing columns: "
            f"{', '.join(missing)}"
        )
    table["ra"].unit = u.deg
    table["dec"].unit = u.deg
    return table


def _column_values(table, column_name, unit=None):
    column = table[column_name]
    if unit is not None and getattr(column, "unit", None) is not None:
        values = column.to_value(unit)
    else:
        values = getattr(column, "value", column)
    return np.asarray(np.ma.filled(values, np.nan))


def _rows_inside_bounds(ra_values, dec_values, bounds):
    if bounds["ra_wraps"]:
        in_ra = (ra_values >= bounds["ra_start"]) | (
            ra_values <= bounds["ra_end"]
        )
    else:
        in_ra = (ra_values >= bounds["ra_start"]) & (
            ra_values <= bounds["ra_end"]
        )
    return (
        in_ra
        & (dec_values >= bounds["dec_min"])
        & (dec_values <= bounds["dec_max"])
        & np.isfinite(ra_values)
        & np.isfinite(dec_values)
    )


def _empty_selected_table():
    return QTable(
        names=("Name", "RA", "DEC", "type", "flux_g", "optical_label"),
        dtype=("U80", float, float, "U8", float, int),
        units=(None, u.deg, u.deg, None, None, None),
    )


def select_dr10_counterparts(
    catalogue,
    detected_markers,
    optical_header,
    hi_mask,
    galaxy_types,
) -> QTable:
    """Select the brightest allowed DR10 row in each H I-overlapping marker.

    A candidate must lie both inside the projected parent H I footprint and
    inside a positive cyan optical-segmentation label.  At most one row is
    returned per label; finite ``flux_g`` values outrank missing values.
    """
    if catalogue is None or len(catalogue) == 0:
        return _empty_selected_table()
    required = ("type", "ra", "dec", "flux_g")
    missing = [column for column in required if column not in catalogue.colnames]
    if missing:
        raise InputError(
            "The DR10 catalogue cannot be matched to optical detections because "
            f"it is missing columns: {', '.join(missing)}"
        )

    marker_data = np.asarray(np.ma.getdata(detected_markers))
    hi_data = np.asarray(np.ma.getdata(hi_mask))
    if marker_data.shape != hi_data.shape:
        raise InputError(
            "The optical segmentation and projected H I mask must have the same "
            "shape before DR10 counterparts are selected."
        )
    height, width = marker_data.shape
    optical_wcs = WCS(optical_header)
    bounds = sky_bounds_from_wcs(optical_wcs, width, height)

    ra_values = _column_values(catalogue, "ra", u.deg)
    dec_values = _column_values(catalogue, "dec", u.deg)
    allowed_types = set(normalise_galaxy_types(galaxy_types))
    spatial_rows = np.flatnonzero(
        _rows_inside_bounds(ra_values, dec_values, bounds)
    )
    spatial_types = np.asarray(
        [
            str(value).strip().upper()
            for value in catalogue["type"][spatial_rows]
        ]
    )
    candidate_rows = spatial_rows[
        np.isin(spatial_types, list(allowed_types))
    ]
    if candidate_rows.size == 0:
        return _empty_selected_table()

    x_values, y_values = optical_wcs.celestial.world_to_pixel_values(
        ra_values[candidate_rows], dec_values[candidate_rows]
    )
    finite_pixels = np.isfinite(x_values) & np.isfinite(y_values)
    x_pixels = np.full(x_values.shape, -1, dtype=int)
    y_pixels = np.full(y_values.shape, -1, dtype=int)
    x_pixels[finite_pixels] = np.floor(
        x_values[finite_pixels] + 0.5
    ).astype(int)
    y_pixels[finite_pixels] = np.floor(
        y_values[finite_pixels] + 0.5
    ).astype(int)
    inside_image = (
        finite_pixels
        & (x_pixels >= 0)
        & (x_pixels < width)
        & (y_pixels >= 0)
        & (y_pixels < height)
    )
    candidate_rows = candidate_rows[inside_image]
    x_pixels = x_pixels[inside_image]
    y_pixels = y_pixels[inside_image]
    if candidate_rows.size == 0:
        return _empty_selected_table()

    flux_values = _column_values(catalogue[candidate_rows], "flux_g")
    best_by_label = {}
    for row_index, x_pixel, y_pixel, flux in zip(
        candidate_rows, x_pixels, y_pixels, flux_values
    ):
        label = int(marker_data[y_pixel, x_pixel])
        if label <= 0 or not bool(hi_data[y_pixel, x_pixel]):
            continue
        flux = float(flux)
        rank = flux if np.isfinite(flux) else -np.inf
        previous = best_by_label.get(label)
        if previous is None or rank > previous[0]:
            best_by_label[label] = (rank, int(row_index))

    if not best_by_label:
        return _empty_selected_table()

    selected_indices = [
        best_by_label[label][1] for label in sorted(best_by_label)
    ]
    selected = QTable(catalogue[selected_indices], copy=True)
    selected.rename_column("ra", "RA")
    selected.rename_column("dec", "DEC")
    if getattr(selected["RA"], "unit", None) is None:
        selected["RA"] = np.asarray(selected["RA"], dtype=float) * u.deg
    else:
        selected["RA"] = selected["RA"].to(u.deg)
    if getattr(selected["DEC"], "unit", None) is None:
        selected["DEC"] = np.asarray(selected["DEC"], dtype=float) * u.deg
    else:
        selected["DEC"] = selected["DEC"].to(u.deg)
    names = []
    for row in selected:
        if "brickname" in selected.colnames and "objid" in selected.colnames:
            names.append(f"LS_{row['brickname']}_{row['objid']}")
        elif "ls_id" in selected.colnames:
            names.append(f"LS_{row['ls_id']}")
        else:
            names.append("Legacy_Surveys_DR10_source")
    selected.add_column(names, name="Name", index=0)
    selected["optical_label"] = sorted(best_by_label)
    return selected


def _project_parent_mask(parent_mask, expected_shape):
    """Project a parent H I mask to the moment-0 plane."""
    mask = np.asarray(np.ma.getdata(parent_mask))
    if mask.ndim == 3:
        mask = np.any(mask > 0, axis=0)
    elif mask.ndim == 2:
        mask = mask > 0
    else:
        raise InputError(
            "The parent H I mask used by the DR10 moment-0 peak filter must "
            f"be 2-D or 3-D; found shape {mask.shape}."
        )
    if mask.shape != expected_shape:
        raise InputError(
            "The parent H I mask and moment-0 image used by the DR10 peak "
            f"filter have different shapes: {mask.shape} and {expected_shape}."
        )
    return mask


def _beam_maximum_filter_footprint(moment0_header, moment0_wcs):
    """Build a beam-FWHM elliptical footprint in moment-0 pixel space."""
    missing = [
        keyword
        for keyword in ("BMAJ", "BMIN", "BPA")
        if keyword not in moment0_header
    ]
    if missing:
        raise InputError(
            "input.filter_dr10_markers_by_moment0_peaks=true requires "
            "synthesized-beam metadata in the parent moment-0 header; "
            f"missing {', '.join(missing)}."
        )
    try:
        bmaj = float(moment0_header["BMAJ"])
        bmin = float(moment0_header["BMIN"])
        bpa = float(moment0_header["BPA"])
    except (TypeError, ValueError) as error:
        raise InputError(
            "The parent moment-0 BMAJ, BMIN, and BPA values must be finite "
            "numbers for DR10 moment-0 peak filtering."
        ) from error
    if (
        not np.isfinite(bmaj)
        or not np.isfinite(bmin)
        or not np.isfinite(bpa)
        or bmaj <= 0
        or bmin <= 0
    ):
        raise InputError(
            "The parent moment-0 synthesized beam is invalid for DR10 peak "
            f"filtering: BMAJ={bmaj!r}, BMIN={bmin!r}, BPA={bpa!r}."
        )

    try:
        pixel_matrix = np.asarray(moment0_wcs.pixel_scale_matrix, dtype=float)
        singular_values = np.linalg.svd(pixel_matrix, compute_uv=False)
    except Exception as error:
        raise InputError(
            "Could not derive the celestial pixel scale from the parent "
            "moment-0 WCS for DR10 peak filtering."
        ) from error
    if (
        pixel_matrix.shape != (2, 2)
        or not np.all(np.isfinite(pixel_matrix))
        or not np.all(np.isfinite(singular_values))
        or np.min(singular_values) <= 0
    ):
        raise InputError(
            "The parent moment-0 celestial WCS has an invalid or singular "
            "pixel-scale matrix for DR10 peak filtering."
        )

    # A maximum-filter window one full beam wide has half-width BMAJ/2.
    # BPA is measured north through east.  The WCS pixel-scale matrix maps
    # pixel offsets to the projected east/north celestial axes, including
    # image rotation and the usual reversed right-ascension axis.
    maximum_radius = int(
        np.ceil(max(bmaj, bmin) / (2.0 * np.min(singular_values)))
    )
    maximum_radius = max(maximum_radius, 1)
    y_offsets, x_offsets = np.mgrid[
        -maximum_radius : maximum_radius + 1,
        -maximum_radius : maximum_radius + 1,
    ]
    east_offsets = (
        pixel_matrix[0, 0] * x_offsets
        + pixel_matrix[0, 1] * y_offsets
    )
    north_offsets = (
        pixel_matrix[1, 0] * x_offsets
        + pixel_matrix[1, 1] * y_offsets
    )
    position_angle = np.deg2rad(bpa)
    major_offsets = (
        east_offsets * np.sin(position_angle)
        + north_offsets * np.cos(position_angle)
    )
    minor_offsets = (
        east_offsets * np.cos(position_angle)
        - north_offsets * np.sin(position_angle)
    )
    footprint = (
        (major_offsets / (bmaj / 2.0)) ** 2
        + (minor_offsets / (bmin / 2.0)) ** 2
        <= 1.0 + 1e-12
    )
    footprint[maximum_radius, maximum_radius] = True
    return footprint


def detect_positive_beam_scale_moment0_peaks(
    moment0_data,
    moment0_header,
    parent_mask,
):
    """Return one representative pixel per positive beam-scale local maximum.

    No global flux or signal-to-noise threshold is applied.  Connected pixels
    belonging to a flat-topped maximum are collapsed to one deterministic
    representative pixel.
    """
    data = np.squeeze(np.asarray(np.ma.getdata(moment0_data), dtype=float))
    if data.ndim != 2:
        raise InputError(
            "The parent moment-0 image used by the DR10 peak filter must be "
            f"2-D after singleton axes are removed; found shape {data.shape}."
        )
    try:
        celestial_wcs = WCS(moment0_header).celestial
    except Exception as error:
        raise InputError(
            "Could not read the parent moment-0 celestial WCS required by "
            "input.filter_dr10_markers_by_moment0_peaks=true."
        ) from error
    if not celestial_wcs.has_celestial:
        raise InputError(
            "input.filter_dr10_markers_by_moment0_peaks=true requires a valid "
            "celestial WCS in the parent moment-0 header."
        )

    inside_parent = _project_parent_mask(parent_mask, data.shape)
    footprint = _beam_maximum_filter_footprint(moment0_header, celestial_wcs)
    eligible = inside_parent & np.isfinite(data) & (data > 0)
    filtered_data = np.where(eligible, data, -np.inf)
    neighbourhood_maximum = maximum_filter(
        filtered_data,
        footprint=footprint,
        mode="constant",
        cval=-np.inf,
    )
    candidate_pixels = eligible & (filtered_data == neighbourhood_maximum)

    peak_map = np.zeros(data.shape, dtype=np.uint8)
    for plateau_value in np.unique(data[candidate_pixels]):
        plateau_labels, plateau_count = label_connected(
            candidate_pixels & (data == plateau_value),
            structure=np.ones((3, 3), dtype=np.uint8),
        )
        for plateau_label in range(1, plateau_count + 1):
            y_pixels, x_pixels = np.nonzero(
                plateau_labels == plateau_label
            )
            if x_pixels.size == 0:
                continue
            centre_x = float(np.mean(x_pixels))
            centre_y = float(np.mean(y_pixels))
            distances = (
                (x_pixels - centre_x) ** 2
                + (y_pixels - centre_y) ** 2
            )
            # np.argmin supplies a stable row-major tie break after nonzero().
            representative = int(np.argmin(distances))
            peak_map[y_pixels[representative], x_pixels[representative]] = 1
    return peak_map


def _empty_peak_filter_audit(selected):
    audit = QTable(selected, copy=True)
    audit["moment0_peak_supported"] = np.zeros(len(audit), dtype=bool)
    audit["status"] = np.full(len(audit), "rejected", dtype="U8")
    audit["peak_x"] = np.full(len(audit), np.nan, dtype=float)
    audit["peak_y"] = np.full(len(audit), np.nan, dtype=float)
    audit["peak_ra"] = np.full(len(audit), np.nan, dtype=float) * u.deg
    audit["peak_dec"] = np.full(len(audit), np.nan, dtype=float) * u.deg
    audit["peak_value"] = np.full(len(audit), np.nan, dtype=float)
    audit["rejection_reason"] = np.full(
        len(audit),
        "no_positive_beam_scale_moment0_peak_in_optical_region",
        dtype="U80",
    )
    return audit


def _map_moment0_peaks_to_optical_labels(
    peak_map,
    moment0_header,
    detected_markers,
    optical_header,
):
    """Map representative moment-0 peak pixels into optical labels."""
    peaks = np.squeeze(np.asarray(np.ma.getdata(peak_map))) > 0
    if peaks.ndim != 2:
        raise InputError(
            "The moment-0 peak map must be 2-D before it can be mapped into "
            f"cyan optical regions; found shape {peaks.shape}."
        )
    marker_data = np.squeeze(
        np.asarray(np.ma.getdata(detected_markers), dtype=int)
    )
    if marker_data.ndim != 2:
        raise InputError(
            "The optical segmentation used by the DR10 peak filter must be "
            f"2-D; found shape {marker_data.shape}."
        )
    try:
        moment0_wcs = WCS(moment0_header).celestial
        optical_wcs = WCS(optical_header).celestial
    except Exception as error:
        raise InputError(
            "Could not read the moment-0 or optical celestial WCS needed to "
            "associate H I peaks with cyan regions."
        ) from error
    if not moment0_wcs.has_celestial or not optical_wcs.has_celestial:
        raise InputError(
            "Valid celestial WCS metadata is required in both the moment-0 "
            "and optical headers to associate H I peaks with cyan regions."
        )

    peak_y, peak_x = np.nonzero(peaks)
    peak_ra, peak_dec = moment0_wcs.pixel_to_world_values(peak_x, peak_y)
    optical_x, optical_y = optical_wcs.world_to_pixel_values(peak_ra, peak_dec)
    finite_optical = np.isfinite(optical_x) & np.isfinite(optical_y)
    rounded_x = np.full(optical_x.shape, -1, dtype=int)
    rounded_y = np.full(optical_y.shape, -1, dtype=int)
    rounded_x[finite_optical] = np.rint(optical_x[finite_optical]).astype(int)
    rounded_y[finite_optical] = np.rint(optical_y[finite_optical]).astype(int)
    inside_optical = (
        finite_optical
        & (rounded_x >= 0)
        & (rounded_x < marker_data.shape[1])
        & (rounded_y >= 0)
        & (rounded_y < marker_data.shape[0])
    )
    peak_labels = np.zeros(peak_x.size, dtype=int)
    peak_labels[inside_optical] = marker_data[
        rounded_y[inside_optical], rounded_x[inside_optical]
    ]
    return {
        "peak_x": peak_x,
        "peak_y": peak_y,
        "peak_ra": np.asarray(peak_ra, dtype=float),
        "peak_dec": np.asarray(peak_dec, dtype=float),
        "optical_x": np.asarray(optical_x, dtype=float),
        "optical_y": np.asarray(optical_y, dtype=float),
        "optical_label": peak_labels,
    }


def optical_labels_with_multiple_moment0_peaks(
    peak_map,
    moment0_header,
    detected_markers,
    optical_header,
    minimum_peaks=2,
):
    """Return cyan labels containing at least ``minimum_peaks`` H I peaks."""
    try:
        minimum_peaks = int(minimum_peaks)
    except (TypeError, ValueError) as error:
        raise InputError("minimum_peaks must be a positive integer.") from error
    if minimum_peaks < 1:
        raise InputError("minimum_peaks must be a positive integer.")
    mapped = _map_moment0_peaks_to_optical_labels(
        peak_map,
        moment0_header,
        detected_markers,
        optical_header,
    )
    labels, counts = np.unique(
        mapped["optical_label"][mapped["optical_label"] > 0],
        return_counts=True,
    )
    return [
        int(label)
        for label, count in zip(labels, counts)
        if int(count) >= minimum_peaks
    ]


def filter_dr10_counterparts_by_moment0_peaks(
    selected,
    detected_markers,
    optical_header,
    moment0_data,
    moment0_header,
    parent_mask,
    peak_map=None,
):
    """Keep selected DR10 rows whose exact cyan region contains an H I peak.

    Returns ``(accepted_rows, audit_table, peak_map)``.  Peak pixels are mapped
    from the moment-0 WCS into the optical-segmentation WCS, so the two images
    may have different pixel sizes, rotations, and array shapes.
    """
    if selected is None or not hasattr(selected, "colnames"):
        raise InputError(
            "The DR10 moment-0 peak filter received no selected catalogue table."
        )
    if "optical_label" not in selected.colnames:
        raise InputError(
            "The selected DR10 table has no optical_label column required by "
            "the moment-0 peak filter."
        )
    if peak_map is None:
        peak_map = detect_positive_beam_scale_moment0_peaks(
            moment0_data,
            moment0_header,
            parent_mask,
        )
    else:
        peak_map = np.squeeze(np.asarray(np.ma.getdata(peak_map))) > 0
    moment0_array = np.squeeze(
        np.asarray(np.ma.getdata(moment0_data), dtype=float)
    )
    if peak_map.shape != moment0_array.shape:
        raise InputError(
            "The supplied moment-0 peak map and parent moment-0 image have "
            f"different shapes: {peak_map.shape} and {moment0_array.shape}."
        )
    audit = _empty_peak_filter_audit(selected)
    if len(selected) == 0 or not np.any(peak_map):
        return QTable(selected[:0], copy=True), audit, peak_map

    mapped_peaks = _map_moment0_peaks_to_optical_labels(
        peak_map,
        moment0_header,
        detected_markers,
        optical_header,
    )
    peak_x = mapped_peaks["peak_x"]
    peak_y = mapped_peaks["peak_y"]
    peak_ra = mapped_peaks["peak_ra"]
    peak_dec = mapped_peaks["peak_dec"]
    peak_labels = mapped_peaks["optical_label"]
    peak_values = moment0_array[peak_y, peak_x]

    accepted_mask = np.zeros(len(selected), dtype=bool)
    for row_index, optical_label in enumerate(selected["optical_label"]):
        matching_peaks = np.flatnonzero(peak_labels == int(optical_label))
        if matching_peaks.size == 0:
            continue
        best_peak = matching_peaks[
            int(np.argmax(peak_values[matching_peaks]))
        ]
        accepted_mask[row_index] = True
        audit["moment0_peak_supported"][row_index] = True
        audit["status"][row_index] = "accepted"
        audit["peak_x"][row_index] = float(peak_x[best_peak])
        audit["peak_y"][row_index] = float(peak_y[best_peak])
        audit["peak_ra"][row_index] = float(peak_ra[best_peak]) * u.deg
        audit["peak_dec"][row_index] = float(peak_dec[best_peak]) * u.deg
        audit["peak_value"][row_index] = float(peak_values[best_peak])
        audit["rejection_reason"][row_index] = ""

    return QTable(selected[accepted_mask], copy=True), audit, peak_map


def remove_rejected_optical_regions(detected_markers, audit_table):
    """Return a marker copy with rejected DR10-associated labels zeroed."""
    marker_data = np.ma.asarray(detected_markers).copy()
    if audit_table is None or len(audit_table) == 0:
        return marker_data
    for row in audit_table:
        if not bool(row["moment0_peak_supported"]):
            marker_data[marker_data == int(row["optical_label"])] = 0
    return marker_data


def write_moment0_peak_filter_debug_products(
    debug_directory,
    source_id,
    moment0_data,
    moment0_header,
    peak_map,
    audit_table,
):
    """Write the exact parent map, peak map, and DR10 decision audit."""
    output_directory = Path(debug_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    moment0_array = np.squeeze(
        np.asarray(np.ma.getdata(moment0_data), dtype=float)
    )
    parent_name = output_directory / (
        f"parent_moment0_used_for_dr10_peak_filter_source_{source_id}.fits"
    )
    peak_name = output_directory / (
        f"moment0_peak_map_source_{source_id}.fits"
    )
    audit_name = output_directory / (
        f"dr10_moment0_peak_filter_audit_source_{source_id}.ecsv"
    )
    fits.writeto(
        parent_name,
        moment0_array,
        header=moment0_header,
        overwrite=True,
        output_verify="silentfix",
    )
    peak_header = moment0_header.copy()
    peak_header["BUNIT"] = ("1", "Binary positive beam-scale peak map")
    fits.writeto(
        peak_name,
        np.asarray(peak_map, dtype=np.uint8),
        header=peak_header,
        overwrite=True,
        output_verify="silentfix",
    )
    QTable(audit_table, copy=True).write(
        audit_name, format="ascii.ecsv", overwrite=True
    )
    return parent_name, peak_name, audit_name
