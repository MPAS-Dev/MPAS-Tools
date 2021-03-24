import multiprocessing


def create_pool(process_count=None, method='forkserver'):
    """
    Crate a pool for creating masks with Python multiprocessing.  This should
    be called only once at the beginning of the script performing cell culling.
    ``pool.terminate()`` should be called before exiting the script.

    Parameters
    ----------
    process_count : int, optional
        The number of processors or None to use all available processors

    method : {'fork', 'spawn', 'forkserver'}
        The mutiprocessing method

    Returns
    -------
    pool : multiprocessing.Pool
        A pool to use for python-based mask creation.
    """
    pool = None
    multiprocessing.set_start_method(method)
    if process_count is None:
        process_count = multiprocessing.cpu_count()
    else:
        process_count = min(process_count, multiprocessing.cpu_count())

    if process_count > 1:
        pool = multiprocessing.Pool(process_count)

    return pool

