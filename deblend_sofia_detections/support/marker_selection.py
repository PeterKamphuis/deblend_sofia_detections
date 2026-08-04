"""Helpers for choosing which optical markers seed the watershed."""

import numpy as np


def catalogue_marker_base(detected_markers, manual_markers_only=False):
    """Return automatic markers or an empty masked base for catalogue-only mode."""
    if not manual_markers_only:
        return detected_markers

    return np.ma.masked_array(
        np.zeros_like(np.ma.getdata(detected_markers)),
        mask=np.ma.getmaskarray(detected_markers).copy(),
    )
