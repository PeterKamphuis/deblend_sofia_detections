"""Per-source visual diagnostics for optical and H I deblending."""

from pathlib import Path

import numpy as np


PURPLE = "#8b2a8d"
CONTOUR_PURPLE = "#d8b4fe"
OPTICAL_MARKER_COLOUR = "#35e6e6"
CATALOGUE_COLOUR = "#ffd400"
REJECTED_CATALOGUE_COLOUR = "#ff5a5f"


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
    """Describe trial SoFiA children and their moment-0/cube products."""
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
    velocity_column = next(
        (
            _column_name(table, name)
            for name in ("v_sofia", "v_rad", "v_opt", "v_app")
            if _column_name(table, name) is not None
        ),
        None,
    )
    velocity_unit = (
        getattr(table[velocity_column], "unit", None)
        if velocity_column is not None
        else None
    )
    components = []
    for row in table:
        component_id = _normalise_component_id(row[id_column])
        if not component_id:
            continue
        moment0_name = (
            Path(cubelet_directory)
            / f"{basename}_{component_id}_mom0.fits"
        )
        cube_name = (
            Path(cubelet_directory)
            / f"{basename}_{component_id}_cube.fits"
        )
        mask_name = (
            Path(cubelet_directory)
            / f"{basename}_{component_id}_mask.fits"
        )
        if not moment0_name.is_file():
            print(
                "WARNING: Missing trial H I moment-0 map for component "
                f"{component_id}: {moment0_name}"
            )
            continue

        ra_deg = float("nan")
        dec_deg = float("nan")
        velocity_kms = float("nan")
        if ra_column is not None and dec_column is not None:
            try:
                ra_deg = _coordinate_in_degrees(row[ra_column])
                dec_deg = _coordinate_in_degrees(row[dec_column])
            except (TypeError, ValueError):
                pass
        if velocity_column is not None:
            try:
                velocity_kms = _velocity_in_kms(
                    row[velocity_column], velocity_unit
                )
            except (TypeError, ValueError):
                pass
        components.append(
            {
                "id": component_id,
                "moment0_name": str(moment0_name),
                "cube_name": str(cube_name) if cube_name.is_file() else None,
                "mask_name": str(mask_name) if mask_name.is_file() else None,
                "ra_deg": ra_deg,
                "dec_deg": dec_deg,
                "velocity_kms": velocity_kms,
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
    peak_supported_column = _column_name(
        table, "moment0_peak_supported"
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
        position = {
            "ra_deg": ra,
            "dec_deg": dec,
            "name": name,
            "catalogue": catalogue,
        }
        if peak_supported_column is not None:
            position["marker_status"] = (
                "accepted"
                if bool(row[peak_supported_column])
                else "rejected"
            )
        positions.append(position)
    return positions


def matched_counterpart_positions(match_table):
    """Extract the manual, DR10, or NED match for one trial SoFiA source."""
    if match_table is None or not hasattr(match_table, "colnames"):
        return []
    for prefix, catalogue, name_column, confirmation_name in (
        ("Manual", "Manual match", "Manual_Name", "Manual_spectroscopic"),
        ("DR10", "Legacy Surveys DR10 match", "DR10_Name", "DR10_counterpart"),
        ("NED", "NED match", "NED_Object Name", "NED_spectroscopic"),
    ):
        confirmed_column = _column_name(
            match_table, confirmation_name
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
            if position.get("marker_status") in ("accepted", "rejected"):
                candidate["marker_status"] = position["marker_status"]
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


def _velocity_in_kms(value, unit=None):
    """Return a finite velocity value in km/s without guessing its unit."""
    import astropy.units as u

    if np.ma.is_masked(value):
        return float("nan")
    if hasattr(value, "to_value"):
        converted = float(value.to_value(u.km / u.s))
    elif unit is not None:
        converted = float((float(value) * u.Unit(unit)).to_value(u.km / u.s))
    else:
        raise ValueError("A velocity unit is required.")
    return converted if np.isfinite(converted) else float("nan")


def _quiet_wcs(header):
    """Construct a WCS without repeating harmless FITS-fix warnings per plot."""
    import warnings
    from astropy.wcs import FITSFixedWarning, WCS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header)


def _spectral_velocity_axis(cube_header, channel_count):
    """Convert the cube spectral WCS to a radio-velocity axis in km/s."""
    import astropy.units as u

    try:
        spectral_wcs = _quiet_wcs(cube_header).spectral
        spectral_values = np.asarray(
            spectral_wcs.pixel_to_world_values(
                np.arange(channel_count, dtype=float)
            ),
            dtype=float,
        )
        unit_name = spectral_wcs.world_axis_units[0]
        spectral_unit = u.Unit(unit_name)
    except Exception as error:
        raise ValueError(
            "Could not derive the spectral coordinate from the cube WCS."
        ) from error
    if spectral_values.shape != (channel_count,) or not np.all(
        np.isfinite(spectral_values)
    ):
        raise ValueError("The cube spectral WCS produced invalid coordinates.")

    spectral_quantity = spectral_values * spectral_unit
    if spectral_unit.is_equivalent(u.m / u.s):
        return spectral_quantity.to_value(u.km / u.s)

    rest_frequency = cube_header.get(
        "RESTFRQ", cube_header.get("RESTFREQ")
    )
    if rest_frequency is None:
        raise ValueError(
            "A frequency/wavelength spectral cube requires RESTFRQ or RESTFREQ "
            "to create velocity PV plots."
        )
    try:
        rest = float(rest_frequency) * u.Hz
        return spectral_quantity.to_value(
            u.km / u.s,
            equivalencies=u.doppler_radio(rest),
        )
    except Exception as error:
        raise ValueError(
            "Could not convert the cube spectral WCS to radio velocity."
        ) from error


def _spatial_coordinate_axis(cube_header, shape, spatial_axis):
    """Return RA or Dec pixel-centre coordinates for a 3-D cube."""
    _, height, width = shape
    try:
        celestial_wcs = _quiet_wcs(cube_header).celestial
        if not celestial_wcs.has_celestial:
            raise ValueError("no celestial WCS")
        if spatial_axis == "ra":
            pixels = np.arange(width, dtype=float)
            coordinates, _ = celestial_wcs.pixel_to_world_values(
                pixels, np.full(width, (height - 1.0) / 2.0)
            )
            # Keep a continuous longitude axis for fields crossing RA=0.
            coordinates = np.rad2deg(
                np.unwrap(np.deg2rad(np.asarray(coordinates, dtype=float)))
            )
        elif spatial_axis == "dec":
            pixels = np.arange(height, dtype=float)
            _, coordinates = celestial_wcs.pixel_to_world_values(
                np.full(height, (width - 1.0) / 2.0), pixels
            )
            coordinates = np.asarray(coordinates, dtype=float)
        else:
            raise ValueError("spatial_axis must be 'ra' or 'dec'.")
    except Exception as error:
        if isinstance(error, ValueError) and "spatial_axis" in str(error):
            raise
        raise ValueError(
            f"Could not derive the {spatial_axis.upper()} coordinate from "
            "the cube celestial WCS."
        ) from error
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(
            f"The cube celestial WCS produced invalid {spatial_axis.upper()} "
            "coordinates."
        )
    return coordinates


def pv_projection_data(cube_data, cube_header, source_mask, spatial_axis):
    """Build a parent-mask-aware RA/Dec-versus-velocity PV projection.

    The cube is summed over the orthogonal celestial axis. ``background`` uses
    all finite cubelet voxels, while ``source`` and ``support`` use only voxels
    inside the supplied 3-D parent/source mask.
    """
    cube = np.squeeze(
        np.asarray(np.ma.getdata(cube_data), dtype=float)
    )
    if cube.ndim != 3:
        raise ValueError(
            f"Expected a 3-D cube for PV plotting, found shape {cube.shape}."
        )
    mask = np.squeeze(np.asarray(np.ma.getdata(source_mask)))
    if mask.ndim == 2:
        mask = np.broadcast_to(mask > 0, cube.shape)
    elif mask.ndim == 3:
        mask = mask > 0
    else:
        raise ValueError(
            "The source mask for PV plotting must be 2-D or 3-D; "
            f"found shape {mask.shape}."
        )
    if mask.shape != cube.shape:
        raise ValueError(
            "The source mask and cube used for PV plotting have different "
            f"shapes: {mask.shape} and {cube.shape}."
        )

    if spatial_axis == "ra":
        collapse_axis = 1
    elif spatial_axis == "dec":
        collapse_axis = 2
    else:
        raise ValueError("spatial_axis must be 'ra' or 'dec'.")

    finite = np.isfinite(cube)
    background = np.sum(
        np.where(finite, cube, 0.0), axis=collapse_axis
    )
    source = np.sum(
        np.where(finite & mask, cube, 0.0), axis=collapse_axis
    )
    support = np.any(mask, axis=collapse_axis)
    spatial_coordinates = _spatial_coordinate_axis(
        cube_header, cube.shape, spatial_axis
    )
    velocities = _spectral_velocity_axis(cube_header, cube.shape[0])
    if velocities.size > 1 and velocities[0] > velocities[-1]:
        velocities = velocities[::-1]
        background = background[::-1]
        source = source[::-1]
        support = support[::-1]
    return {
        "spatial_axis": spatial_axis,
        "spatial_coordinates": spatial_coordinates,
        "velocity_kms": velocities,
        "background": background,
        "source": source,
        "support": support,
    }


def _coordinate_edges(coordinates):
    """Convert monotonic pixel-centre coordinates to plotting edges."""
    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("PV plot coordinates must be a non-empty 1-D array.")
    if values.size == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    differences = np.diff(values)
    if not np.all(np.isfinite(differences)) or not (
        np.all(differences > 0) or np.all(differences < 0)
    ):
        raise ValueError("PV plot coordinates must be strictly monotonic.")
    edges = np.empty(values.size + 1, dtype=float)
    edges[1:-1] = (values[:-1] + values[1:]) / 2.0
    edges[0] = values[0] - differences[0] / 2.0
    edges[-1] = values[-1] + differences[-1] / 2.0
    return edges


def _align_ra_to_axis(ra_deg, spatial_coordinates):
    reference = float(np.mean(spatial_coordinates))
    return reference + ((float(ra_deg) - reference + 180.0) % 360.0) - 180.0


def _pv_optical_centres(optical_image_name, marker_data):
    """Return world-coordinate centroids of positive cyan marker labels."""
    from astropy.io import fits
    optical_wcs = _quiet_wcs(fits.getheader(optical_image_name)).celestial
    if not optical_wcs.has_celestial:
        raise ValueError("The optical marker image has no celestial WCS.")
    centres = []
    for x_position, y_position, label in marker_centroids(marker_data):
        ra_deg, dec_deg = optical_wcs.pixel_to_world_values(
            x_position, y_position
        )
        if np.isfinite(ra_deg) and np.isfinite(dec_deg):
            centres.append((float(ra_deg), float(dec_deg), label))
    return centres


def _add_pv_position_guides(
    axes,
    projection,
    optical_centres,
    catalogue_positions,
):
    """Draw spatial-only optical/catalogue locations as vertical PV guides."""
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory

    spatial_axis = projection["spatial_axis"]
    coordinates = projection["spatial_coordinates"]
    minimum, maximum = sorted((float(coordinates[0]), float(coordinates[-1])))

    def spatial_value(ra_deg, dec_deg):
        if spatial_axis == "ra":
            return _align_ra_to_axis(ra_deg, coordinates)
        return float(dec_deg)

    visible_optical = []
    for ra_deg, dec_deg, label in optical_centres:
        value = spatial_value(ra_deg, dec_deg)
        if minimum <= value <= maximum:
            visible_optical.append((value, label))
            axes.axvline(
                value,
                color=OPTICAL_MARKER_COLOUR,
                linewidth=1.0,
                linestyle=":",
                alpha=0.75,
                zorder=5,
            )

    visible_positions = []
    text_transform = blended_transform_factory(
        axes.transData, axes.transAxes
    )
    for index, position in enumerate(deduplicate_positions(catalogue_positions)):
        value = spatial_value(position["ra_deg"], position["dec_deg"])
        if not minimum <= value <= maximum:
            continue
        rejected = position.get("marker_status") == "rejected"
        colour = (
            REJECTED_CATALOGUE_COLOUR if rejected else CATALOGUE_COLOUR
        )
        axes.axvline(
            value,
            color=colour,
            linewidth=1.5,
            linestyle="--" if rejected else "-",
            alpha=0.90,
            zorder=7,
        )
        axes.text(
            value,
            0.02 + 0.05 * (index % 2),
            position["name"],
            transform=text_transform,
            rotation=90,
            rotation_mode="anchor",
            ha="left",
            va="bottom",
            fontsize=6.5,
            color=colour,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "#191919",
                "edgecolor": "none",
                "alpha": 0.60,
            },
            zorder=8,
        )
        visible_positions.append(position)

    legend_items = []
    if visible_optical:
        legend_items.append(
            Line2D(
                [0], [0], color=OPTICAL_MARKER_COLOUR, linestyle=":",
                linewidth=1.5, label="Detected optical source position"
            )
        )
    if any(
        position.get("marker_status") != "rejected"
        for position in visible_positions
    ):
        legend_items.append(
            Line2D(
                [0], [0], color=CATALOGUE_COLOUR, linewidth=1.7,
                label="Accepted catalogue position (no velocity)"
            )
        )
    if any(
        position.get("marker_status") == "rejected"
        for position in visible_positions
    ):
        legend_items.append(
            Line2D(
                [0], [0], color=REJECTED_CATALOGUE_COLOUR,
                linestyle="--", linewidth=1.7,
                label="Rejected catalogue position (no velocity)"
            )
        )
    return visible_optical, visible_positions, legend_items


def _load_component_pv(component, spatial_axis):
    from astropy.io import fits

    cube_name = component.get("cube_name")
    mask_name = component.get("mask_name")
    if not cube_name or not mask_name:
        raise ValueError("missing child cube or mask")
    child_cube, child_header = fits.getdata(cube_name, header=True)
    child_mask = fits.getdata(mask_name)
    return pv_projection_data(
        child_cube, child_header, child_mask, spatial_axis
    )


def write_pv_debug_overlay(
    *,
    cube_data,
    cube_header,
    source_mask,
    optical_image_name,
    marker_data,
    catalogue_positions,
    output_name,
    source_id,
    spatial_axis,
    marker_mode="automatic + manual catalogue",
    components=None,
):
    """Write one RA/Dec-versus-velocity QA projection."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter
    from astropy.coordinates import Angle
    import astropy.units as u

    projection = pv_projection_data(
        cube_data, cube_header, source_mask, spatial_axis
    )
    spatial_coordinates = projection["spatial_coordinates"]
    velocities = projection["velocity_kms"]
    x_edges = _coordinate_edges(spatial_coordinates)
    y_edges = _coordinate_edges(velocities)

    figure, axes = plt.subplots(figsize=(9, 7))
    lower, upper = _display_limits(projection["background"])
    axes.pcolormesh(
        x_edges,
        y_edges,
        projection["background"],
        cmap="gray",
        vmin=lower,
        vmax=upper,
        shading="flat",
        rasterized=True,
    )
    footprint = np.ma.masked_where(
        ~projection["support"],
        np.ones(projection["support"].shape, dtype=float),
    )
    axes.pcolormesh(
        x_edges,
        y_edges,
        footprint,
        cmap=ListedColormap([PURPLE]),
        alpha=0.40,
        shading="flat",
        rasterized=True,
    )
    parent_levels = contour_levels(
        projection["source"], projection["support"]
    )
    if parent_levels:
        axes.contour(
            spatial_coordinates,
            velocities,
            np.ma.masked_where(
                ~projection["support"], projection["source"]
            ),
            levels=parent_levels,
            colors=CONTOUR_PURPLE,
            linewidths=0.8,
            alpha=0.70,
            zorder=4,
        )

    optical_centres = _pv_optical_centres(
        optical_image_name, marker_data
    )
    visible_optical, visible_positions, guide_legend = (
        _add_pv_position_guides(
            axes,
            projection,
            optical_centres,
            catalogue_positions,
        )
    )

    usable_components = []
    component_legend = []
    component_centres = 0
    component_list = components or []
    colours = component_colour_mapping(
        [component.get("id", "") for component in component_list]
    )
    for component in component_list:
        component_id = _normalise_component_id(component.get("id", ""))
        try:
            child_projection = _load_component_pv(
                component, spatial_axis
            )
            child_levels = contour_levels(
                child_projection["source"], child_projection["support"]
            )
            if not child_levels:
                raise ValueError("no finite, varying child PV values")
        except Exception as error:
            print(
                "WARNING: Skipping trial H I component "
                f"{component_id or 'unknown'} in the {spatial_axis.upper()} "
                f"PV overlay: {type(error).__name__}: {error}"
            )
            continue
        colour = colours[component_id]
        axes.contour(
            child_projection["spatial_coordinates"],
            child_projection["velocity_kms"],
            np.ma.masked_where(
                ~child_projection["support"], child_projection["source"]
            ),
            levels=child_levels,
            colors=colour,
            linewidths=1.35,
            alpha=0.95,
            zorder=6,
        )
        usable_components.append(component)
        component_legend.append(
            Line2D(
                [0], [0], color=colour, linewidth=1.6, marker="o",
                markerfacecolor=colour, markeredgecolor="#191919",
                markersize=5.5, label=f"H I child {component_id}"
            )
        )
        try:
            spatial_value = (
                _align_ra_to_axis(
                    component["ra_deg"], spatial_coordinates
                )
                if spatial_axis == "ra"
                else float(component["dec_deg"])
            )
            velocity = float(component["velocity_kms"])
            minimum, maximum = sorted(
                (float(spatial_coordinates[0]), float(spatial_coordinates[-1]))
            )
            centre_is_valid = (
                np.isfinite(spatial_value)
                and np.isfinite(velocity)
                and minimum <= spatial_value <= maximum
                and velocities[0] <= velocity <= velocities[-1]
            )
        except (KeyError, TypeError, ValueError):
            centre_is_valid = False
        if centre_is_valid:
            axes.scatter(
                spatial_value,
                velocity,
                marker="o",
                s=58,
                color=colour,
                edgecolor="#191919",
                linewidth=1.0,
                zorder=9,
            )
            axes.annotate(
                f"H I {component_id}",
                (spatial_value, velocity),
                xytext=(6, -10),
                textcoords="offset points",
                fontsize=7,
                color=colour,
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": "#191919",
                    "edgecolor": "none",
                    "alpha": 0.68,
                },
                zorder=10,
            )
            component_centres += 1

    spatial_name = "Right Ascension" if spatial_axis == "ra" else "Declination"
    collapsed_name = "Declination" if spatial_axis == "ra" else "Right Ascension"
    if spatial_axis == "ra":
        axes.xaxis.set_major_formatter(
            FuncFormatter(
                lambda value, _: Angle(value * u.deg).to_string(
                    unit=u.hourangle, sep=":", precision=0, pad=True
                )
            )
        )
    else:
        axes.xaxis.set_major_formatter(
            FuncFormatter(
                lambda value, _: Angle(value * u.deg).to_string(
                    unit=u.deg, sep=":", precision=0, alwayssign=True
                )
            )
        )
    axes.set_xlabel(spatial_name)
    axes.set_ylabel(r"Velocity (km s$^{-1}$)")
    axes.grid(color="white", alpha=0.18, linestyle=":")
    component_mode = bool(component_list)
    if component_mode:
        axes.set_title(
            f"SoFiA source {source_id}: raw candidate H I component QA — "
            f"{spatial_name} / velocity PV\n"
            f"{len(usable_components)} candidate component(s), "
            f"{component_centres} measured centre(s); summed over "
            f"{collapsed_name}\nWatershed markers: {marker_mode}"
        )
    else:
        axes.set_title(
            f"SoFiA source {source_id}: catalogue QA — {spatial_name} / "
            f"velocity PV\n{len(visible_optical)} optical detection(s), "
            f"{len(visible_positions)} catalogue position(s); summed over "
            f"{collapsed_name}\nWatershed markers: {marker_mode}"
        )
    legend_items = [
        Patch(facecolor="0.45", label="Parent H I cube projection"),
        Patch(facecolor=PURPLE, alpha=0.40, label="Parent H I footprint"),
        Line2D(
            [0], [0], color=CONTOUR_PURPLE, linewidth=1.0,
            label="Parent H I PV contours"
        ),
    ]
    legend_items.extend(guide_legend)
    legend_items.extend(component_legend)
    axes.legend(
        handles=legend_items,
        loc="upper right",
        fontsize=7,
        ncol=2 if len(legend_items) > 7 else 1,
    )
    axes.set_xlim(float(x_edges[0]), float(x_edges[-1]))
    axes.set_ylim(float(y_edges[0]), float(y_edges[-1]))

    output_path = Path(output_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return output_path


def write_source_pv_debug_overlays(*, output_ra_name, output_dec_name, **kwargs):
    """Write catalogue QA projections in RA-velocity and Dec-velocity."""
    ra_path = write_pv_debug_overlay(
        output_name=output_ra_name,
        spatial_axis="ra",
        **kwargs,
    )
    dec_path = write_pv_debug_overlay(
        output_name=output_dec_name,
        spatial_axis="dec",
        **kwargs,
    )
    return ra_path, dec_path


def write_source_pv_debug_overlays_safely(**kwargs):
    """Write catalogue PV QA without allowing plotting to stop a run."""
    try:
        output_paths = write_source_pv_debug_overlays(**kwargs)
        print(
            "Wrote catalogue PV debug overlays to "
            f"{output_paths[0]} and {output_paths[1]}."
        )
        return output_paths
    except Exception as error:
        try:
            from matplotlib import pyplot as plt

            plt.close("all")
        except Exception:
            pass
        source_id = kwargs.get("source_id", "unknown")
        print(
            "WARNING: Could not create the catalogue PV debug overlays for "
            f"SoFiA source {source_id}: {type(error).__name__}: {error}"
        )
        return None


def write_hi_component_pv_debug_overlays(
    *, output_ra_name, output_dec_name, **kwargs
):
    """Write raw-child QA projections in RA-velocity and Dec-velocity."""
    ra_path = write_pv_debug_overlay(
        output_name=output_ra_name,
        spatial_axis="ra",
        **kwargs,
    )
    dec_path = write_pv_debug_overlay(
        output_name=output_dec_name,
        spatial_axis="dec",
        **kwargs,
    )
    return ra_path, dec_path


def write_hi_component_pv_debug_overlays_safely(**kwargs):
    """Write raw-child PV QA without allowing plotting to stop a run."""
    components = kwargs.get("components") or []
    if not components:
        return None
    try:
        output_paths = write_hi_component_pv_debug_overlays(**kwargs)
        print(
            "Wrote raw trial H I component PV debug overlays to "
            f"{output_paths[0]} and {output_paths[1]}."
        )
        return output_paths
    except Exception as error:
        try:
            from matplotlib import pyplot as plt

            plt.close("all")
        except Exception:
            pass
        source_id = kwargs.get("source_id", "unknown")
        print(
            "WARNING: Could not create the raw trial H I component PV debug "
            f"overlays for SoFiA source {source_id}: "
            f"{type(error).__name__}: {error}"
        )
        return None


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
        rejected = position.get("marker_status") == "rejected"
        scatter_options = {
            "marker": "x" if rejected else "o",
            "s": 78 if rejected else 70,
            "color": (
                REJECTED_CATALOGUE_COLOUR
                if rejected
                else CATALOGUE_COLOUR
            ),
            "linewidth": 1.6 if rejected else 0.8,
            "zorder": 6,
        }
        if not rejected:
            scatter_options["edgecolor"] = "#191919"
        axes.scatter(
            coordinate.ra.deg,
            coordinate.dec.deg,
            transform=world_transform,
            **scatter_options,
        )
        axes.annotate(
            position["name"],
            (coordinate.ra.deg, coordinate.dec.deg),
            xycoords=world_transform,
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=7,
            color=(
                REJECTED_CATALOGUE_COLOUR
                if rejected
                else CATALOGUE_COLOUR
            ),
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
    accepted_positions = [
        position
        for position in plotted_positions
        if position.get("marker_status") != "rejected"
    ]
    rejected_positions = [
        position
        for position in plotted_positions
        if position.get("marker_status") == "rejected"
    ]
    if accepted_positions:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=CATALOGUE_COLOUR,
                markeredgecolor="#191919",
                markersize=7,
                label="Accepted catalogue position",
            )
        )
    if rejected_positions:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="x",
                color="none",
                markeredgecolor=REJECTED_CATALOGUE_COLOUR,
                markeredgewidth=1.6,
                markersize=7,
                label="Rejected by moment-0 peak filter",
            )
        )
    if not plotted_positions:
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
        rejected = position.get("marker_status") == "rejected"
        scatter_options = {
            "marker": "x" if rejected else "o",
            "s": 78 if rejected else 70,
            "color": (
                REJECTED_CATALOGUE_COLOUR
                if rejected
                else CATALOGUE_COLOUR
            ),
            "linewidth": 1.6 if rejected else 0.8,
            "zorder": 7,
        }
        if not rejected:
            scatter_options["edgecolor"] = "#191919"
        axes.scatter(
            coordinate.ra.deg,
            coordinate.dec.deg,
            transform=world_transform,
            **scatter_options,
        )
        axes.annotate(
            position["name"],
            (coordinate.ra.deg, coordinate.dec.deg),
            xycoords=world_transform,
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=7,
            color=(
                REJECTED_CATALOGUE_COLOUR
                if rejected
                else CATALOGUE_COLOUR
            ),
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
    accepted_positions = [
        position
        for position in plotted_positions
        if position.get("marker_status") != "rejected"
    ]
    rejected_positions = [
        position
        for position in plotted_positions
        if position.get("marker_status") == "rejected"
    ]
    if accepted_positions:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=CATALOGUE_COLOUR,
                markeredgecolor="#191919",
                markersize=7,
                label="Accepted catalogue position",
            )
        )
    if rejected_positions:
        legend_items.append(
            Line2D(
                [0],
                [0],
                marker="x",
                color="none",
                markeredgecolor=REJECTED_CATALOGUE_COLOUR,
                markeredgewidth=1.6,
                markersize=7,
                label="Rejected by moment-0 peak filter",
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
