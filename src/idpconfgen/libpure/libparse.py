"""
Parse data from one format to another.

This module contains pure functions that parse data from one format to
another.

Functions in this module should only depend on Python STD.
"""


def parse_number_ranges(string):
    """
    Parse a number range in string to a tuple of tuples.

    Numbers should be integers or they will be assumed as integers.

    Example
    -------
    >>> parse_number_range('1-59,80-102,150-170')
    ((1, 59), (80, 102), (150, 170))
    """
    err_msg = f'Not a valid input: {string!r}. Please read docstring.'
    # removes empty strings for the cases where single range is given
    input_ranges = (i for i in string.split(',') if i)
    ranges = []
    for range_ in input_ranges:
        try:
            start, end = map(int, range_.split('-'))
        except TypeError as err:
            raise ValueError(err_msg) from err

        ranges.append((start, end))

    if input_ranges:
        return tuple(ranges)
    else:
        raise ValueError(err_msg)
