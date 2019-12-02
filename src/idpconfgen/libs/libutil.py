"""General utilities."""
import random


def random_fragment(iterable, fragsize=None):
    """
    Generate a slice object that reflects a random fragment from iterable.
    
    (iterable, int -> slice object)
    
    Parameters
    ----------
    iterable : iterable-type
        An interable: string, list, tuple, etc.

    fragsize : int or NoneType
        The size of the fragment to generate.

    Returns
    -------
    slice object
        A slice object reflecting a random fragment of the `iterable`
        of size `fragsize`.
        If `fragsize` is ``None``, returns a full range slice object.
    """
    
    # this is separate from the try: block to account input of for types
    # that are now iterable and have not to do with the usage
    # of fragsize=None
    len_ = len(iterable)

    try:
        start = random.randint(0, len_ - fragsize)
    except TypeError:  # fragsize is None
        return slice(None, None, None)
    else:
        return slice(start, start + fragsize, None)
