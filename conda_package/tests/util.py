import os


def get_test_data_file(filename):
    """
    Get the full path to a data file in the tests/data directory.

    Parameters
    ----------
    filename : str
        The name of the data file.

    Returns
    -------
    str
        The full relative path to the data file.
    """

    local_path = os.path.join(
        'mesh_tools', 'mesh_conversion_tools', 'test', filename
    )
    repo_path = os.path.join('..', '..', local_path)
    if os.path.exists(local_path):
        return local_path
    elif os.path.exists(repo_path):
        return repo_path

    raise FileNotFoundError(
        f"Data file '{filename}' not found in expected locations: "
        f'{local_path} or {repo_path}'
    )
