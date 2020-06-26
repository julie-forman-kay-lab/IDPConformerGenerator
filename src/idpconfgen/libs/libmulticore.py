"""Multi-core related objects."""
from functools import partial
from multiprocessing import Pool

from idpconfgen import log
from idpconfgen.libs.libtimer import ProgressWatcher


# USED OKAY
# here can't use closures because those can't be pickled
def consume_iterable_in_list(func, *args, **kwargs):
    """
    Consumes a generator to a list.

    Used as a wrapper to send generators to Python multiprocessing lib.

    Returns
    -------
    list of length N
        Where N is the number of times `func` yields.
    """
    return list(func(*args, **kwargs))


# USED OKAY
def flat_results_from_chunk(execute, func, *args, **kwargs):
    """
    Flattens a result coming from consume_iterable_in_list execution.

    Parameters
    ----------
    execute : callable
        The prepared multiprocessing function.
        Each item of its iteration should be iterable of iterable.

    func : callable
        The function to process each yielded result from the results
        in the multiprocessing chunk.
    """
    for chunk in execute():
        # each result is a list coming from consume_iterable_in_list
        # where each element is a list yield from
        # extract_secondary_structure
        flatted = (yielded for result in chunk for yielded in result)
        func(flatted, *args, **kwargs)


# USED OKAY
def pool_function(
        func,
        items,
        *args,
        method='imap_unordered',
        ncores=1,
        **kwargs,
        ):
    """
    Multiprocess Pools a function.

    Parameters
    ----------
    func : callable
        The function to execute along the iterable `items`.

    items : interable
        Elements to pass to the `func`.

    *args : any
        Additional positional arguments to pass to `func`.
        In order to workd, positional arguments should preceed `item`
        in function signature. Is preferable to use `kwargs` whenever
        possible.

    method : str
        The :class:`Pool` method to execute.
        Defaults to `imap_unordered`.

    ncores : int
        The number of cores to use. Defaults to `1`.

    **kwargs :
        The named arguments to pass to `func`.
    """
    f = partial(func, *args, **kwargs)

    with Pool(ncores) as pool, ProgressWatcher(items) as pb:
        imap = getattr(pool, method)(f, items)
        # the use of `while` here is needed, instead of for
        # to allo try/catch options
        while True:
            try:
                yield next(imap)
                pb.increment()
            except StopIteration:
                break
            except IndexError:
                log.info('IndexError of multiprocessing, ignoring something')


# USED OKAY
def pool_function_in_chunks(func, items, *args, chunks=5_000, **kwargs):
    """
    Execute ``func`` in ``chunks`` of ``items`` using `Pool`.

    Yields the results after each chunk.

    Parameters
    ----------
    func : callable
        The function to execute over ``items``.

    items : iterable
        The items to process by the ``funct``.

    *args : any
        Additional positional arguments to send to ``func``.

    chunks : int
        The size of each chunk processed multiprocessing before yielding.

    **kwargs : any
        Additional keyword arguments to send to ``func``.

    Yields
    ------
    list
        Containing the results after each chunk.
    """
    for i in range(0, len(items), chunks):
        task = items[i: i + chunks]
        yield list(pool_function(func, task, *args, **kwargs))
