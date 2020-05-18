"""Manages time."""
import functools
import os
import sys
import time

import numpy as np

from idpconfgen import log
from idpconfgen.logger import S


def record_time(process_name='', *args, **kwargs):
    """
    Record time of function execution.

    Use as decorator.
    """
    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            start = time.time()

            result = func(*args, **kwargs)

            log.info(S(f'elapsed time :{process_name}: {time.time() - start}'))

            return result
        return wrapper
    return decorator


class ProgressBar:
    """
    Contextualizes a Progress Bar.

    Parameters
    ----------
    total : int convertable
        The total number o iteractions expected.

    prefix : str
        Some prefix to enhance readability.

    suffix : str
        Some suffix to enhance readability.

    decimals : int-convertable
        The demicals to show in percentage.
        Defaults to `1`.

    bar_length : int, float, -convertable
        The length of the bar.
        If not provided (``None``), uses half of the terminal window.

    Thanks to for the initial template function:
    https://dev.to/natamacm/progressbar-in-python-a3n

    Examples
    --------

    >>> with ProgressBar(5000, suffix='frames') as PB:
    >>>     for i in trajectory:
    >>>         # do something
    >>>         PB.increment()
    """

    def __init__(
            self,
            total,
            prefix='',
            suffix='',
            decimals=1,
            bar_length=None,
            ):

        if bar_length is None:
            try:
                _columns, _rows = os.get_terminal_size()
            except OSError:
                log.error(
                    'ERROR: Could not retrive size of ProgressBar '
                    'from terminal window. Using the default of `50`. '
                    'Everything else is normal.'
                    )
                # this trick is used to guarantee 100% test coverage
                _columns = 100
            bar_length = _columns // 2

        total = int(total)
        self.prefix = prefix
        self.suffix = suffix
        self.str_format = "{0:." + str(int(decimals)) + "f}"

        # using Numpy
        self.percentages = np.linspace(0, 100, total + 1, endpoint=True)
        # 49.7 µs ± 5.34 µs per loop (7 runs, 10000 loops each)
        # Not using Numpy
        # self.percentages =  [i / total * 100 for i in range(total + 1)]
        # 974 µs ± 38.8 µs per loop (7 runs, 1000 loops each)

        self.filled_length = \
            np.round(bar_length * self.percentages / 100).astype(int)
        self.counter = 0
        self.total = total
        self.bar_length = bar_length

        assert len(self.percentages) == total + 1
        assert len(self.percentages) == len(self.filled_length)

    def __enter__(self):
        bar = '-' * self.bar_length
        percents = self.str_format.format(self.percentages[0])
        sys.stdout.write(
            f'\r{self.prefix} |{bar}| '
            f'{percents}% {self.counter}/{self.total} {self.suffix}'
            )
        self.counter += 1
        return self

    def __exit__(self, *args):
        sys.stdout.write('\n')
        sys.stdout.flush()

    def increment(self):
        """Print next progress bar increment."""
        t = self.total
        c = self.counter
        prefix = self.prefix
        suffix = self.suffix
        bl = self.bar_length
        percents = self.str_format.format(self.percentages[c])
        fl = self.filled_length[c]
        bar = f"{'█' * fl}{'-' * (bl - fl)}"
        sys.stdout.write(f'\r{prefix} |{bar}| {percents}% {c}/{t} {suffix}')
        self.counter += 1
