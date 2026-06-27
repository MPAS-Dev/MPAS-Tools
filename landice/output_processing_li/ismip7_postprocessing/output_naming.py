"""Helpers for building ISMIP7-compliant output filenames."""

import os


REQUIRED_NAMING_KEYS = [
    'domain_id',
    'source_id',
    'ism_id',
    'ism_member_id',
    'esm_id',
    'forcing_member_id',
    'experiment_id',
    'set_counter',
    'time_range',
]


def _sanitize_component(value):
    """Sanitize one filename component into a compact token."""
    text = str(value).strip()
    text = text.replace(' ', '-')
    text = text.replace('/', '-')
    return text


def build_output_filename(output_path, variable_id, metadata):
    """
    Build a filename of the form:
    <variable_id>_<domain_id>_<source_id>_<ism_id>_<ISM_member_id>_<ESM_id>_
    <forcing_member_id>_<experiment_id>_<set_counter>_<time_range>.nc
    """
    missing = [key for key in REQUIRED_NAMING_KEYS if key not in metadata]
    if missing:
        raise ValueError(
            "Metadata is missing required naming fields: "
            f"{missing}"
        )

    pieces = [
        _sanitize_component(variable_id),
        _sanitize_component(metadata['domain_id']),
        _sanitize_component(metadata['source_id']),
        _sanitize_component(metadata['ism_id']),
        _sanitize_component(metadata['ism_member_id']),
        _sanitize_component(metadata['esm_id']),
        _sanitize_component(metadata['forcing_member_id']),
        _sanitize_component(metadata['experiment_id']),
        _sanitize_component(metadata['set_counter']),
        _sanitize_component(metadata['time_range']),
    ]
    return os.path.join(output_path, '_'.join(pieces) + '.nc')
