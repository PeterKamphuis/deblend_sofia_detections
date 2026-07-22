"""Helpers for selecting which SoFiA catalogue detections to process."""

from deblend_sofia_detections.support.errors import InputError


def select_source_ids(catalogue_source_ids, requested_source_ids):
    """Return catalogue IDs selected by an optional user allowlist.

    The returned values retain their catalogue types and catalogue ordering so they
    can be used directly when constructing cubelet filenames. Matching is performed
    on string representations because Astropy tables and OmegaConf may represent
    the same SoFiA ID with different scalar types.
    """

    catalogue_source_ids = list(catalogue_source_ids)
    if not requested_source_ids:
        return catalogue_source_ids

    available_by_text = {
        str(source_id).strip(): source_id for source_id in catalogue_source_ids
    }
    requested_as_text = {
        str(source_id).strip() for source_id in requested_source_ids
    }

    missing_source_ids = requested_as_text - set(available_by_text)
    if missing_source_ids:
        missing = sorted(missing_source_ids)
        raise InputError(
            "The following requested source IDs are not present in the SoFiA "
            f"catalogue: {missing}"
        )

    return [
        source_id
        for source_id in catalogue_source_ids
        if str(source_id).strip() in requested_as_text
    ]
