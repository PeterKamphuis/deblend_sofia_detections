"""Per-source visual diagnostics for optical and H I deblending."""

from pathlib import Path

import numpy as np


PURPLE = "#8b2a8d"
CONTOUR_PURPLE = "#d8b4fe"
OPTICAL_MARKER_COLOUR = "#35e6e6"
CATALOGUE_COLOUR = "#ffd400"


def marker_centroids(marker_data):
    """Return ``(x, y, label)`` centroids for positive marker labels."""
    markers = np.ma.asarray(marker_data).filled(0)
    centroids = []
    for label in np.unique(markers):
        if label <= 0:
            continue
        y_positions, x_positions = np.nonzero(markers == label)
        if x_positions.size:
            centroids.append(
                (
                    float(np.mean(x_positions)),
                    float(np.mean(y_positions)),
                    int(label),
                )
            )
    return centroids


def contour_levels(moment0_data, source_mask, percentiles=(60, 75, 90)):
    """Calculate stable, unique contour levels inside the source footprint."""
    data = np.asarray(moment0_data, dtype=float)
    mask = np.asarray(source_mask, dtype=bool)
    values = data[mask & np.isfinite(data)]
    positive_values = values[values > 0]
    if positive_values.size:
        values = positive_values
    if values.size < 2 or np.nanmin(values) == np.nanmax(values):
        return []
    levels = np.percentile(values, percentiles)
    return [float(level) for level in np.unique(levels) if np.isfinite(level)]


def _column_name(table, requested_name):
    for column in table.colnames:
        if column.lower() == requested_name.lower():
            return column
    return None


def _coordinate_in_degrees(value):
    if np.ma.is_masked(value):
        return float("nan")
    if hasattr(value, "to_value"):
        try:
            return float(value.to_value("deg"))
        except Exception:
            pass
    if hasattr(value, "value"):
        value = value.value
    return float(value)


def catalogue_positions_from_table(table, catalogue="Manual input"):
    """Extract plottable positions from a manual or query-result table."""
    if table is None or not hasattr(table, "colnames"):
        return []
    ra_column = _column_name(table, "RA")
    dec_column = _column_name(table, "DEC")
    if ra_column is None or dec_column is None:
        return []
    name_column = (
        _column_name(table, "Name")
        or _column_name(table, "Object Name")
        or _column_name(table, "ID")
    )
    positions = []
    for index, row in enumerate(table):
        try:
            ra = _coordinate_in_degrees(row[ra_column])
            dec = _coordinate_in_degrees(row[dec_column])
        except (TypeError, ValueError):
            continue
        if not np.isfinite(ra) or not np.isfinite(dec):
            continue
        name = (
            str(row[name_column])
            if name_column is not None
            else f"{catalogue} source {index + 1}"
        )
        positions.append(
            {
                "ra_deg": ra,
                "dec_deg": dec,
                "name": name,
                "catalogue": catalogue,
            }
        )
    return positions


def matched_counterpart_positions(match_table):
    """Extract the manual/NED counterpart selected for one trial SoFiA source."""
    if match_table is None or not hasattr(match_table, "colnames"):
        return []
    for prefix, catalogue, name_column in (
        ("Manual", "Manual match", "Manual_Name"),
        ("NED", "NED match", "NED_Object Name"),
    ):
        confirmed_column = _column_name(
            match_table, f"{prefix}_spectroscopic"
        )
        if confirmed_column is None:
            continue
        confirmed = np.asarray(match_table[confirmed_column]).reshape(-1)
        if confirmed.size == 0 or not bool(confirmed[0]):
            continue
        ra_column = _column_name(match_table, f"{prefix}_RA")
        dec_column = _column_name(match_table, f"{prefix}_DEC")
        if ra_column is None or dec_column is None:
            continue
        ra = _coordinate_in_degrees(match_table[ra_column][0])
        dec = _coordinate_in_degrees(match_table[dec_column][0])
        if not np.isfinite(ra) or not np.isfinite(dec):
            continue
        actual_name_column = _column_name(match_table, name_column)
        name = (
            str(match_table[actual_name_column][0])
            if actual_name_column is not None
            else catalogue
        )
        return [
            {
                "ra_deg": ra,
                "dec_deg": dec,
                "name": name,
                "catalogue": catalogue,
            }
        ]
    return []


def deduplicate_positions(positions, tolerance_degrees=1e-7):
    """Remove repeated positions while preferring later match information."""
    unique = []
    for position in positions or []:
        try:
            candidate = {
                "ra_deg": float(position["ra_deg"]),
                "dec_deg": float(position["dec_deg"]),
                "name": str(position.get("name", "Catalogue source")),
                "catalogue": str(position.get("catalogue", "Catalogue")),
            }
        except (KeyError, TypeError, ValueError):
            continue
        if not np.isfinite(candidate["ra_deg"]) or not np.isfinite(
            candidate["dec_deg"]
        ):
            continue
        duplicate_index = None
        for index, existing in enumerate(unique):
            if (
                abs(candidate["ra_deg"] - existing["ra_deg"])
                <= tolerance_degrees
                and abs(candidate["dec_deg"] - existing["dec_deg"])
                <= tolerance_degrees
            ):
                duplicate_index = index
                break
        if duplicate_index is None:
            unique.append(candidate)
        else:
            unique[duplicate_index] = candidate
    return unique


def _positions_inside_image(positions, optical_wcs, shape):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    selected = []
    height, width = shape
    for position in deduplicate_positions(positions):
        coordinate = SkyCoord(
            position["ra_deg"] * u.deg,
            position["dec_deg"] * u.deg,
            frame="icrs",
        )
        x_position, y_position = optical_wcs.world_to_pixel(coordinate)
        if (
            np.isfinite(x_position)
            and np.isfinite(y_position)
            and -0.5 <= x_position < width - 0.5
            and -0.5 <= y_position < height - 0.5
        ):
            selected.append(position)
    return selected


def _display_limits(data):
    values = np.asarray(data, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    lower, upper = np.percentile(values, (1, 99.5))
    if lower == upper:
        upper = lower + 1.0
    return float(lower), float(upper)


def write_source_debug_overlay(
    *,
    optical_image_name,
    moment0_data,
    moment0_header,
    source_mask,
    marker_data,
    catalogue_positions,
    output_name,
    catalogue_output_name,
    source_id,
):
    """Write the optical/H I/catalogue QA PNG and its plotted-position table."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.table import Table
    from astropy.wcs import WCS
    import astropy.units as u

    optical_data = np.squeeze(fits.getdata(optical_image_name))
    optical_header = fits.getheader(optical_image_name)
    if optical_data.ndim != 2:
        raise ValueError(
            f"Expected a 2-D optical image, found shape {optical_data.shape}."
        )
    optical_wcs = WCS(optical_header).celestial

    hi_data = np.squeeze(np.asarray(moment0_data, dtype=float))
    if hi_data.ndim != 2:
        raise ValueError(
            f"Expected a 2-D H I moment-0 image, found shape {hi_data.shape}."
        )
    hi_wcs = WCS(moment0_header).celestial

    hi_mask = np.asarray(source_mask)
    if hi_mask.ndim == 3:
        hi_mask = np.any(hi_mask > 0, axis=0)
    else:
        hi_mask = hi_mask > 0
    if hi_mask.shape != hi_data.shape:
        raise ValueError(
            "The projected source mask and moment-0 image have different "
            f"shapes: {hi_mask.shape} and {hi_data.shape}."
        )

    plotted_positions = _positions_inside_image(
        catalogue_positions, optical_wcs, optical_data.shape
    )
    position_table = Table(
        rows=[
            (
                position["catalogue"],
                position["name"],
                position["ra_deg"],
                position["dec_deg"],
            )
            for position in plotted_positions
        ],
        names=("catalogue", "name", "ra_deg", "dec_deg"),
        dtype=("U32", "U128", float, float),
    )
    position_table["ra_deg"].unit = u.deg
    position_table["dec_deg"].unit = u.deg
    catalogue_path = Path(catalogue_output_name)
    catalogue_path.parent.mkdir(parents=True, exist_ok=True)
    position_table.write(catalogue_path, format="ascii.ecsv", overwrite=True)

    figure = plt.figure(figsize=(8, 8))
    axes = figure.add_subplot(111, projection=optical_wcs)
    lower, upper = _display_limits(optical_data)
    axes.imshow(
        optical_data,
        origin="lower",
        cmap="gray",
        vmin=lower,
        vmax=upper,
        interpolation="nearest",
    )

    hi_transform = axes.get_transform(hi_wcs)
    footprint = np.ma.masked_where(
        ~hi_mask, np.ones(hi_mask.shape, dtype=float)
    )
    axes.imshow(
        footprint,
        origin="lower",
        cmap=ListedColormap([PURPLE]),
        alpha=0.40,
        interpolation="nearest",
        transform=hi_transform,
    )
    levels = contour_levels(hi_data, hi_mask)
    if levels:
        axes.contour(
            np.ma.masked_where(~hi_mask, hi_data),
            levels=levels,
            colors=CONTOUR_PURPLE,
            linewidths=0.8,
            alpha=0.65,
            origin="lower",
            transform=hi_transform,
        )

    detected_markers = np.ma.asarray(marker_data).filled(0)
    centroids = marker_centroids(detected_markers)
    if centroids:
        for _, _, marker_label in centroids:
            axes.contour(
                detected_markers == marker_label,
                levels=[0.5],
                colors=OPTICAL_MARKER_COLOUR,
                linewidths=1.0,
                alpha=0.75,
                origin="lower",
            )
        axes.scatter(
            [centroid[0] for centroid in centroids],
            [centroid[1] for centroid in centroids],
            marker="x",
            s=70,
            linewidths=1.8,
            color=OPTICAL_MARKER_COLOUR,
            zorder=5,
        )

    world_transform = axes.get_transform("world")
    for position in plotted_positions:
        coordinate = SkyCoord(
            position["ra_deg"] * u.deg,
            position["dec_deg"] * u.deg,
            frame="icrs",
        )
        axes.scatter(
            coordinate.ra.deg,
            coordinate.dec.deg,
            transform=world_transform,
            marker="o",
            s=70,
            color=CATALOGUE_COLOUR,
            edgecolor="#191919",
            linewidth=0.8,
            zorder=6,
        )
        axes.annotate(
            position["name"],
            (coordinate.ra.deg, coordinate.dec.deg),
            xycoords=world_transform,
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=7,
            color=CATALOGUE_COLOUR,
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "#191919",
                "edgecolor": "none",
                "alpha": 0.65,
            },
            zorder=7,
        )

    axes.coords[0].set_axislabel("Right Ascension")
    axes.coords[1].set_axislabel("Declination")
    axes.coords.grid(color="white", alpha=0.18, linestyle=":")
    axes.set_title(
        f"SoFiA source {source_id}: optical / H I / catalogue QA\n"
        f"{len(centroids)} optical detection(s), "
        f"{len(plotted_positions)} catalogue position(s)"
    )
    legend_items = [
        Patch(facecolor=PURPLE, alpha=0.40, label="H I source footprint"),
        Line2D(
            [0],
            [0],
            color=CONTOUR_PURPLE,
            linewidth=1.0,
            label="H I moment-0 contours",
        ),
    ]
    if centroids:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="x",
                color="none",
                markeredgecolor=OPTICAL_MARKER_COLOUR,
                markeredgewidth=1.8,
                markersize=7,
                label="Detected optical source",
            )
        )
    if plotted_positions:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=CATALOGUE_COLOUR,
                markeredgecolor="#191919",
                markersize=7,
                label="Catalogue position",
            )
        )
    else:
        axes.text(
            0.02,
            0.02,
            "No catalogue position found in this cutout",
            transform=axes.transAxes,
            color=CATALOGUE_COLOUR,
            fontsize=8,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "#191919",
                "edgecolor": "none",
                "alpha": 0.70,
            },
            zorder=7,
        )
    axes.set_xlim(-0.5, optical_data.shape[1] - 0.5)
    axes.set_ylim(-0.5, optical_data.shape[0] - 0.5)
    axes.legend(handles=legend_items, loc="upper right", fontsize=8)

    output_path = Path(output_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return output_path


def write_source_debug_overlay_safely(**kwargs):
    """Write the QA products without allowing a plotting error to stop a run."""
    try:
        output_path = write_source_debug_overlay(**kwargs)
        print(
            "Wrote optical/H I/catalogue debug products to "
            f"{output_path} and {kwargs['catalogue_output_name']}."
        )
        return output_path
    except Exception as error:
        try:
            from matplotlib import pyplot as plt

            plt.close("all")
        except Exception:
            pass
        source_id = kwargs.get("source_id", "unknown")
        print(
            "WARNING: Could not create the optical/H I/catalogue debug "
            f"overlay for SoFiA source {source_id}: "
            f"{type(error).__name__}: {error}"
        )
        return None
