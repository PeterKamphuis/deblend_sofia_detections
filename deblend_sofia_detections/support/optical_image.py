"""Helpers for normalising multi-plane optical FITS image data."""

import warnings

import numpy as np


def celestial_numpy_axes(
    axis_types,
    axis_correlation_matrix,
    pixel_n_dim,
    data_ndim,
):
    """Map celestial FITS/WCS pixel axes to NumPy array-axis indices."""
    if int(pixel_n_dim) != int(data_ndim):
        return None
    celestial_world_axes = [
        index
        for index, axis_type in enumerate(axis_types)
        if axis_type.get("coordinate_type") == "celestial"
    ]
    if not celestial_world_axes:
        return None
    correlation = np.asarray(axis_correlation_matrix, dtype=bool)
    if correlation.ndim != 2 or correlation.shape[0] != len(axis_types):
        return None
    celestial_pixel_axes = np.flatnonzero(
        np.any(correlation[celestial_world_axes, :], axis=0)
    )
    if len(celestial_pixel_axes) != 2:
        return None
    return tuple(
        sorted(data_ndim - 1 - int(axis) for axis in celestial_pixel_axes)
    )


def collapse_optical_data(data, celestial_numpy_axes=None):
    """Return a 2-D optical image and the array axes collapsed to create it.

    FITS stores a conventional cube as ``(plane, y, x)``, while some colour
    products use a different position for the RGB axis. When the caller can
    identify the two celestial array axes from WCS, every other axis is
    averaged. Without that information, singleton axes and an obvious RGB(A)
    axis are handled before falling back to the conventional leading-plane
    layout.
    """
    image = np.asanyarray(data)
    if image.ndim < 2:
        raise ValueError(
            f"Optical image data must have at least two dimensions; "
            f"found shape {image.shape}."
        )
    if image.ndim == 2:
        return image, ()

    if celestial_numpy_axes is not None:
        celestial_axes = tuple(
            sorted(int(axis) % image.ndim for axis in celestial_numpy_axes)
        )
        if len(celestial_axes) != 2 or len(set(celestial_axes)) != 2:
            raise ValueError(
                "Exactly two distinct celestial array axes are required; "
                f"received {celestial_numpy_axes}."
            )
        collapse_axes = tuple(
            axis for axis in range(image.ndim) if axis not in celestial_axes
        )
        return _mean_over_axes(image, collapse_axes), collapse_axes

    singleton_axes = tuple(
        axis for axis, length in enumerate(image.shape) if length == 1
    )
    squeezed = np.squeeze(image)
    if squeezed.ndim == 2:
        return squeezed, singleton_axes

    if image.ndim == 3:
        colour_axes = [
            axis for axis, length in enumerate(image.shape) if length in (3, 4)
        ]
        if len(colour_axes) == 1:
            colour_axis = colour_axes[0]
            return _mean_over_axes(image, (colour_axis,)), (colour_axis,)

    # Conventional FITS image cubes have non-spatial planes before (y, x).
    collapse_axes = tuple(range(image.ndim - 2))
    reduced = _mean_over_axes(image, collapse_axes)
    if reduced.ndim != 2:
        raise ValueError(
            f"Could not reduce optical image shape {image.shape} to two "
            "spatial dimensions."
        )
    return reduced, collapse_axes


def _mean_over_axes(data, axes):
    if not axes:
        reduced = np.asanyarray(data)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            reduced = np.nanmean(data, axis=axes, dtype=np.float32)
    if reduced.ndim != 2:
        raise ValueError(
            f"Collapsing array axes {axes} produced shape {reduced.shape}, "
            "not a 2-D optical image."
        )
    return reduced
