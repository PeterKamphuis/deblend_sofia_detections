"""Per-source visual diagnostics for optical and H I deblending."""

from pathlib import Path

import numpy as np


PURPLE = "#8b2a8d"
CONTOUR_PURPLE = "#d8b4fe"
OPTICAL_MARKER_COLOUR = "#35e6e6"
CATALOGUE_COLOUR = "#ffd400"


def _normalise_component_id(value):
    """Return the file/legend representation used for a SoFiA source ID."""
    if np.ma.is_masked(value):
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    try:
        numeric_value = float(value)
        if np.isfinite(numeric_value) and numeric_value.is_integer():
            return str(int(numeric_value))
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _component_sort_key(component_id):
    try:
        return (0, int(component_id))
    except (TypeError, ValueError):
        return (1, str(component_id))


def component_colour_mapping(component_ids):
    """Assign stable, distinct colours to sorted trial-component IDs."""
    import matplotlib
    from matplotlib.colors import to_hex

    normalised_ids = sorted(
        {
            component_id
            for component_id in (
                _normalise_component_id(value) for value in component_ids
            )
            if component_id
        },
        key=_component_sort_key,
    )
    if not normalised_ids:
        return {}

    colour_map = matplotlib.colormaps[
        "tab10" if len(normalised_ids) <= 10 else "hsv"
    ]
    colours = [
        to_hex(
            colour_map(index)
            if len(normalised_ids) <= 10
            else colour_map(index / len(normalised_ids))
        )
        for index in range(len(normalised_ids))
    ]
    return dict(zip(normalised_ids, colours))


def trial_hi_components_from_table(table, cubelet_directory, basename):
    """Describe trial SoFiA children and their moment-0 products."""
    if table is None or not hasattr(table, "colnames"):
        return []
    id_column = _column_name(table, "id")
    if id_column is None:
        print(
            "WARNING: Trial SoFiA catalogue has no ID column; "
            "cannot create the H I component overlay."
        )
        return []

    ra_column = _column_name(table, "ra")
    dec_column = _column_name(table, "dec")
    components = []
    for row in table:
        component_id = _normalise_component_id(row[id_column])
        if not component_id:
            continue
        moment0_name = (
            Path(cubelet_directory)
            / f"{basename}_{component_id}_mom0.fits"
        )
        if not moment0_name.is_file():
            print(
                "WARNING: Missing trial H I moment-0 map for component "
                f"{component_id}: {moment0_name}"
            )
            continue

        ra_deg = float("nan")
        dec_deg = float("nan")
        if ra_column is not None and dec_column is not None:
            try:
                ra_deg = _coordinate_in_degrees(row[ra_column])
                dec_deg = _coordinate_in_degrees(row[dec_column])
            except (TypeError, ValueError):
                pass
        components.append(
            {
                "id": component_id,
                "moment0_name": str(moment0_name),
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
            }
        )

    return sorted(
        components,
        key=lambda component: _component_sort_key(component["id"]),
    )


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
    marker_mode="automatic + manual catalogue",
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
        f"{len(plotted_positions)} catalogue position(s)\n"
        f"Watershed markers: {marker_mode}"
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


def write_hi_component_debug_overlay(
    *,
    optical_image_name,
    moment0_data,
    moment0_header,
    source_mask,
    marker_data,
    catalogue_positions,
    components,
    output_name,
    source_id,
    marker_mode="automatic + manual catalogue",
):
    """Write raw trial-child H I contours and centres on the optical QA view."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS
    import astropy.units as u

    optical_data = np.squeeze(fits.getdata(optical_image_name))
    optical_header = fits.getheader(optical_image_name)
    if optical_data.ndim != 2:
        raise ValueError(
            f"Expected a 2-D optical image, found shape {optical_data.shape}."
        )
    optical_wcs = WCS(optical_header).celestial

    parent_data = np.squeeze(np.asarray(moment0_data, dtype=float))
    if parent_data.ndim != 2:
        raise ValueError(
            "Expected a 2-D parent H I moment-0 image, "
            f"found shape {parent_data.shape}."
        )
    parent_wcs = WCS(moment0_header).celestial

    parent_mask = np.asarray(source_mask)
    if parent_mask.ndim == 3:
        parent_mask = np.any(parent_mask > 0, axis=0)
    else:
        parent_mask = parent_mask > 0
    if parent_mask.shape != parent_data.shape:
        raise ValueError(
            "The projected parent mask and moment-0 image have different "
            f"shapes: {parent_mask.shape} and {parent_data.shape}."
        )

    colours = component_colour_mapping(
        [component.get("id", "") for component in components]
    )
    usable_components = []
    for component in components:
        component_id = _normalise_component_id(component.get("id", ""))
        try:
            child_data = np.squeeze(
                np.asarray(
                    fits.getdata(component["moment0_name"]),
                    dtype=float,
                )
            )
            child_header = fits.getheader(component["moment0_name"])
            if child_data.ndim != 2:
                raise ValueError(
                    f"expected two dimensions, found {child_data.shape}"
                )
            child_wcs = WCS(child_header).celestial
            if not child_wcs.has_celestial:
                raise ValueError("no valid celestial WCS")
            child_mask = np.isfinite(child_data) & (child_data != 0)
            levels = contour_levels(child_data, child_mask)
            if not levels:
                raise ValueError("no finite, varying contour values")
        except Exception as error:
            print(
                "WARNING: Skipping trial H I component "
                f"{component_id or 'unknown'} in the debug overlay: "
                f"{type(error).__name__}: {error}"
            )
            continue

        usable_component = dict(component)
        usable_component.update(
            {
                "id": component_id,
                "data": child_data,
                "mask": child_mask,
                "levels": levels,
                "wcs": child_wcs,
                "colour": colours[component_id],
            }
        )
        usable_components.append(usable_component)

    if not usable_components:
        raise ValueError("No usable trial H I child moment-0 maps were found.")

    plotted_positions = _positions_inside_image(
        catalogue_positions, optical_wcs, optical_data.shape
    )
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

    parent_transform = axes.get_transform(parent_wcs)
    footprint = np.ma.masked_where(
        ~parent_mask, np.ones(parent_mask.shape, dtype=float)
    )
    axes.imshow(
        footprint,
        origin="lower",
        cmap=ListedColormap([PURPLE]),
        alpha=0.40,
        interpolation="nearest",
        transform=parent_transform,
    )
    parent_levels = contour_levels(parent_data, parent_mask)
    if parent_levels:
        axes.contour(
            np.ma.masked_where(~parent_mask, parent_data),
            levels=parent_levels,
            colors=CONTOUR_PURPLE,
            linewidths=0.8,
            alpha=0.55,
            origin="lower",
            transform=parent_transform,
        )

    detected_markers = np.ma.asarray(marker_data).filled(0)
    optical_centroids = marker_centroids(detected_markers)
    if optical_centroids:
        for _, _, marker_label in optical_centroids:
            axes.contour(
                detected_markers == marker_label,
                levels=[0.5],
                colors=OPTICAL_MARKER_COLOUR,
                linewidths=1.0,
                alpha=0.65,
                origin="lower",
            )
        axes.scatter(
            [centroid[0] for centroid in optical_centroids],
            [centroid[1] for centroid in optical_centroids],
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
            zorder=7,
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
            zorder=8,
        )

    component_legend_items = []
    plotted_centres = 0
    for component in usable_components:
        axes.contour(
            np.ma.masked_where(~component["mask"], component["data"]),
            levels=component["levels"],
            colors=component["colour"],
            linewidths=1.35,
            alpha=0.95,
            origin="lower",
            transform=axes.get_transform(component["wcs"]),
            zorder=6,
        )
        component_legend_items.append(
            Line2D(
                [0],
                [0],
                color=component["colour"],
                linewidth=1.5,
                marker="o",
                markerfacecolor=component["colour"],
                markeredgecolor="#191919",
                markersize=6,
                label=f"H I child {component['id']}",
            )
        )

        try:
            centre = SkyCoord(
                float(component["ra_deg"]) * u.deg,
                float(component["dec_deg"]) * u.deg,
                frame="icrs",
            )
            x_centre, y_centre = optical_wcs.world_to_pixel(centre)
            centre_is_inside = (
                np.isfinite(x_centre)
                and np.isfinite(y_centre)
                and -0.5 <= x_centre < optical_data.shape[1] - 0.5
                and -0.5 <= y_centre < optical_data.shape[0] - 0.5
            )
        except (TypeError, ValueError):
            centre_is_inside = False

        if not centre_is_inside:
            print(
                "WARNING: Trial H I component "
                f"{component['id']} has no valid centre inside the optical "
                "cutout; its contours will still be plotted."
            )
            continue

        axes.scatter(
            centre.ra.deg,
            centre.dec.deg,
            transform=world_transform,
            marker="o",
            s=62,
            color=component["colour"],
            edgecolor="#191919",
            linewidth=1.0,
            zorder=9,
        )
        axes.annotate(
            f"H I {component['id']}",
            (centre.ra.deg, centre.dec.deg),
            xycoords=world_transform,
            xytext=(6, -10),
            textcoords="offset points",
            fontsize=7,
            color=component["colour"],
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "#191919",
                "edgecolor": "none",
                "alpha": 0.68,
            },
            zorder=10,
        )
        plotted_centres += 1

    axes.coords[0].set_axislabel("Right Ascension")
    axes.coords[1].set_axislabel("Declination")
    axes.coords.grid(color="white", alpha=0.18, linestyle=":")
    axes.set_title(
        f"SoFiA source {source_id}: raw candidate H I component QA\n"
        f"{len(usable_components)} candidate component(s), "
        f"{plotted_centres} measured centre(s)\n"
        f"Watershed markers: {marker_mode}"
    )
    legend_items = [
        Patch(facecolor=PURPLE, alpha=0.40, label="Parent H I footprint"),
        Line2D(
            [0],
            [0],
            color=CONTOUR_PURPLE,
            linewidth=1.0,
            label="Parent H I moment-0",
        ),
    ]
    if optical_centroids:
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
    legend_items.extend(component_legend_items)
    axes.set_xlim(-0.5, optical_data.shape[1] - 0.5)
    axes.set_ylim(-0.5, optical_data.shape[0] - 0.5)
    axes.legend(
        handles=legend_items,
        loc="upper right",
        fontsize=7,
        ncol=2 if len(legend_items) > 7 else 1,
    )

    output_path = Path(output_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return output_path


def write_hi_component_debug_overlay_safely(**kwargs):
    """Write raw trial-child QA without allowing plotting to stop a run."""
    components = kwargs.get("components") or []
    if not components:
        return None
    try:
        output_path = write_hi_component_debug_overlay(**kwargs)
        print(f"Wrote raw trial H I component debug overlay to {output_path}.")
        return output_path
    except Exception as error:
        try:
            from matplotlib import pyplot as plt

            plt.close("all")
        except Exception:
            pass
        source_id = kwargs.get("source_id", "unknown")
        print(
            "WARNING: Could not create the raw trial H I component debug "
            f"overlay for SoFiA source {source_id}: "
            f"{type(error).__name__}: {error}"
        )
        return None
