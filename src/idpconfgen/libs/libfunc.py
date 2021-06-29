def flatlist(list_):
    """
    Flat a list recursively.

    This is a generator.
    """
    # escape strings which would yield infinite recursion
    if isinstance(list_, str):
        yield list_
    else:
       try:
           for sublist in list_:
               yield from flatlist(sublist)
       except TypeError:  # sublist is not iterable
           yield sublist
