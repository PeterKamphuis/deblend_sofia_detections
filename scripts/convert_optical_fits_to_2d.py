#!/usr/bin/env python3
"""Convert a multi-dimensional optical FITS image to a 2-D celestial FITS.

This script is independent of ``deblend_sofia_detections`` itself. It requires
only NumPy and Astropy, both of which are already dependencies of the package.

Examples
--------
Convert an RGB or multi-band image by averaging non-celestial planes::

    python scripts/convert_optical_fits_to_2d.py input.fits

This automatically writes ``input_2d.fits`` beside the input. An explicit
output path can still be supplied::

    python scripts/convert_optical_fits_to_2d.py input.fits output_2d.fits

Use the median instead of the mean::

    python scripts/convert_optical_fits_to_2d.py \
        input.fits output_2d.fits --method median

Select the first plane along every non-celestial axis::

    python scripts/convert_optical_fits_to_2d.py \
        input.fits output_2d.fits --method first
"""

import argparse
from pathlib import Path
import sys
import warnings

import numpy as np


PRESERVED_HEADER_KEYWORDS = (
    "OBJECT",
    "TELESCOP",
    "INSTRUME",
    "FILTER",
    "BUNIT",
    "DATE-OBS",
    "ORIGIN",
    "AUTHOR",
)


def default_output_path(input_name):
    """Return ``<input>_2d`` beside the input while preserving FITS suffixes."""
    input_path = Path(input_name).expanduser()
    lower_name = input_path.name.lower()
    for suffix in (".fits.gz", ".fit.gz", ".fts.gz", ".fits", ".fit", ".fts"):
        if lower_name.endswith(suffix):
            base_name = input_path.name[: -len(suffix)]
            return input_path.with_name(f"{base_name}_2d{suffix}")
    return input_path.with_name(f"{input_path.stem}_2d.fits")


def celestial_numpy_axes(
    axis_types,
    axis_correlation_matrix,
    pixel_n_dim,
    data_ndim,
):
    """Map celestial WCS pixel axes to NumPy array-axis indices."""
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


def infer_spatial_axes(data, full_wcs):
    """Identify the two NumPy axes that contain the celestial image."""
    try:
        axes = celestial_numpy_axes(
            full_wcs.get_axis_types(),
            full_wcs.axis_correlation_matrix,
            full_wcs.pixel_n_dim,
            data.ndim,
        )
    except (AttributeError, IndexError, TypeError, ValueError):
        axes = None
    if axes is not None:
        return axes, "celestial WCS"

    # A unique RGB(A) axis is a safe fallback for common colour FITS products.
    if data.ndim == 3:
        colour_axes = [
            axis for axis, length in enumerate(data.shape) if length in (3, 4)
        ]
        if len(colour_axes) == 1:
            colour_axis = colour_axes[0]
            spatial_axes = tuple(
                axis for axis in range(data.ndim) if axis != colour_axis
            )
            return spatial_axes, "RGB(A) axis length"

    # Conventional FITS cubes are ordered (..., y, x).
    return (data.ndim - 2, data.ndim - 1), "conventional (..., y, x) layout"


def collapse_to_2d(data, spatial_axes, method="mean", dtype=np.float32):
    """Collapse every non-spatial axis and return a two-dimensional array."""
    image = np.asanyarray(data)
    if image.ndim < 2:
        raise ValueError(
            f"Input image must have at least two dimensions; "
            f"found shape {image.shape}."
        )
    spatial_axes = tuple(sorted(int(axis) % image.ndim for axis in spatial_axes))
    if len(spatial_axes) != 2 or len(set(spatial_axes)) != 2:
        raise ValueError(
            f"Exactly two distinct spatial axes are required; "
            f"received {spatial_axes}."
        )
    collapsed_axes = tuple(
        axis for axis in range(image.ndim) if axis not in spatial_axes
    )

    if not collapsed_axes:
        reduced = image
    elif method == "first":
        selection = [
            0 if axis in collapsed_axes else slice(None)
            for axis in range(image.ndim)
        ]
        reduced = image[tuple(selection)]
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if method == "mean":
                reduced = np.nanmean(
                    image, axis=collapsed_axes, dtype=dtype
                )
            elif method == "median":
                reduced = np.nanmedian(image, axis=collapsed_axes)
            else:
                raise ValueError(f"Unknown collapse method: {method}")

    reduced = np.asarray(reduced, dtype=dtype)
    if reduced.ndim != 2:
        raise ValueError(
            f"Collapsing array axes {collapsed_axes} produced shape "
            f"{reduced.shape}, not a 2-D image."
        )
    return reduced, collapsed_axes


def resolve_hdu(hdul, requested_hdu=None):
    """Return the requested HDU or the first HDU containing image data."""
    if requested_hdu is not None:
        try:
            key = int(requested_hdu)
        except ValueError:
            key = requested_hdu
        try:
            hdu = hdul[key]
        except (IndexError, KeyError) as error:
            raise ValueError(
                f"FITS HDU {requested_hdu!r} was not found."
            ) from error
        if hdu.data is None or np.ndim(hdu.data) < 2:
            raise ValueError(
                f"FITS HDU {requested_hdu!r} does not contain image data "
                "with at least two dimensions."
            )
        return hdu, requested_hdu

    for index, hdu in enumerate(hdul):
        if hdu.data is not None and np.ndim(hdu.data) >= 2:
            return hdu, index
    raise ValueError(
        "The FITS file contains no HDU with image data of at least two "
        "dimensions."
    )


def convert_fits(
    input_name,
    output_name=None,
    *,
    requested_hdu=None,
    method="mean",
    dtype_name="float32",
    overwrite=False,
):
    """Convert one FITS image and return a summary dictionary."""
    from astropy.io import fits
    from astropy.wcs import WCS

    input_path = Path(input_name).expanduser().resolve()
    if output_name is None:
        output_name = default_output_path(input_path)
    output_path = Path(output_name).expanduser().resolve()
    if input_path == output_path:
        raise ValueError(
            "Input and output paths are identical. Choose a different output "
            "file so the original FITS image is preserved."
        )
    if not input_path.is_file():
        raise FileNotFoundError(f"Input FITS image does not exist: {input_path}")
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {output_path}. Use --overwrite to "
            "replace it."
        )

    dtype = np.dtype(dtype_name)
    with fits.open(input_path, memmap=False) as hdul:
        image_hdu, selected_hdu = resolve_hdu(hdul, requested_hdu)
        data = np.asanyarray(image_hdu.data)
        original_shape = tuple(data.shape)
        original_header = image_hdu.header.copy()
        try:
            full_wcs = WCS(original_header)
        except Exception as error:
            raise ValueError(
                f"Input FITS image has an invalid WCS: {error}"
            ) from error
        if not full_wcs.has_celestial:
            raise ValueError(
                "Input FITS image does not contain a celestial WCS."
            )
        spatial_axes, axis_source = infer_spatial_axes(data, full_wcs)
        reduced, collapsed_axes = collapse_to_2d(
            data, spatial_axes, method=method, dtype=dtype
        )
        celestial_header = full_wcs.celestial.to_header(relax=True)

    for keyword in PRESERVED_HEADER_KEYWORDS:
        if keyword in original_header:
            celestial_header[keyword] = original_header[keyword]
    celestial_header.add_history(
        f"Converted to 2-D from {input_path.name}; HDU {selected_hdu}."
    )
    celestial_header.add_history(
        f"Original array shape: {original_shape}; spatial axes: "
        f"{spatial_axes} ({axis_source})."
    )
    celestial_header.add_history(
        f"Collapsed array axes: {collapsed_axes}; method: {method}; "
        f"output dtype: {dtype.name}."
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(
        data=reduced,
        header=celestial_header,
    ).writeto(
        output_path,
        overwrite=overwrite,
        output_verify="fix",
    )

    with fits.open(output_path, memmap=False) as output_hdul:
        output_data = output_hdul[0].data
        output_wcs = WCS(output_hdul[0].header)
        if output_data is None or output_data.ndim != 2:
            raise RuntimeError(
                f"Output verification failed: found shape "
                f"{None if output_data is None else output_data.shape}."
            )
        if not output_wcs.has_celestial:
            raise RuntimeError(
                "Output verification failed: celestial WCS is missing."
            )

    return {
        "input": input_path,
        "output": output_path,
        "hdu": selected_hdu,
        "original_shape": original_shape,
        "output_shape": tuple(reduced.shape),
        "spatial_axes": spatial_axes,
        "collapsed_axes": collapsed_axes,
        "axis_source": axis_source,
        "method": method,
        "dtype": dtype.name,
    }


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a 2-D or multi-dimensional optical FITS image into a "
            "2-D FITS image with celestial WCS."
        )
    )
    parser.add_argument("input_fits", help="Input optical FITS image.")
    parser.add_argument(
        "output_fits",
        nargs="?",
        help=(
            "Optional new 2-D FITS image. If omitted, <input>_2d.fits is "
            "created beside the input. It cannot be the input path."
        ),
    )
    parser.add_argument(
        "--hdu",
        help=(
            "HDU index or extension name. By default, the first HDU with "
            "image data of at least two dimensions is used."
        ),
    )
    parser.add_argument(
        "--method",
        choices=("mean", "median", "first"),
        default="mean",
        help=(
            "How to collapse non-celestial axes: mean (default), median, "
            "or first."
        ),
    )
    parser.add_argument(
        "--dtype",
        choices=("float32", "float64"),
        default="float32",
        help="Output data type (default: float32).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output file. The input is never replaced.",
    )
    return parser


def main(argv=None):
    arguments = build_argument_parser().parse_args(argv)
    try:
        summary = convert_fits(
            arguments.input_fits,
            arguments.output_fits,
            requested_hdu=arguments.hdu,
            method=arguments.method,
            dtype_name=arguments.dtype,
            overwrite=arguments.overwrite,
        )
    except Exception as error:
        print(
            f"ERROR: {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        return 1

    print(f"Input:            {summary['input']}")
    print(f"Selected HDU:     {summary['hdu']}")
    print(f"Original shape:   {summary['original_shape']}")
    print(
        f"Spatial axes:     {summary['spatial_axes']} "
        f"({summary['axis_source']})"
    )
    print(f"Collapsed axes:   {summary['collapsed_axes']}")
    print(f"Collapse method:  {summary['method']}")
    print(f"Output shape:     {summary['output_shape']}")
    print(f"Output dtype:     {summary['dtype']}")
    print(f"Wrote:            {summary['output']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
