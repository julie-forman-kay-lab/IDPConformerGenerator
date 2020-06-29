"""Manages time."""
import functools
import os
import sys
import time

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


class ProgressWatcher:
    """Construct a Progress Watcher context."""

    def __new__(cls, items, *args, **kwargs):
        """
        Construct a progress representation.

        If possible, creates a progress bar based on the terminal
        window size and the length of items.

        If terminal window size is not measurable, mostly because
        execution is running on a server, ignores the creation of the
        bar. That is, returns a dummy representation.

        If items have no length (is a generator), returns a progression
        counter instead.

        Parameters
        ----------
        items : iterable
            Items to represent.
        """
        try:
            _cols, _rows = os.get_terminal_size()
        except OSError:
            log.warning(
                'WARNING: Could not retrieve size from terminal window. '
                'Ignoring progress watcher.'
                )
            return ProgressFake()

        try:
            litems = len(items)
        except TypeError:
            return ProgressCounter(**kwargs)
        else:
            return ProgressBar(litems, *args, bar_length=_cols // 2, **kwargs)


class ProgressBar:
    """
    Contextualize a Progress Bar.

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
        If not provided (``None``), uses 20.

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
            bar_length=20,
            ):

        self.total = int(total)
        self.prefix = prefix
        self.suffix = suffix
        self.percent_format = "{:>5." + str(int(decimals)) + "f}%"
        totals_format = "{:>" + str(len(str(self.total))) + "}"
        self.totals_format = f"{totals_format}/{totals_format}"
        self.counter = 0
        self.bar_length = bar_length

    def __enter__(self):
        bar = '-' * self.bar_length
        percents = self.percent_format.format(0)
        totals = self.totals_format.format(self.counter, self.total)
        sys.stdout.write(
            f'\r{self.prefix} '
            f'{percents} {totals} {self.suffix} |{bar}|'
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
        bl = self.bar_length
        percents = self.percent_format.format(c / t * 100)
        totals = self.totals_format.format(c, t)
        fl = int(round(bl * c // t))
        bar = f"{'█' * fl}{'-' * (bl - fl)}"
        sys.stdout.write(
            f'\r{self.prefix} {percents} {totals} {self.suffix} '
            f'|{bar}|'
            )
        self.counter += 1


class ProgressCounter:
    """Represent progression via a counter."""

    def __init__(self, suffix='', **kwargs):
        """
        Represent a progress counter.

        Used for progresses with unknown length.
        """
        self.counter = 0
        self.suffix = suffix

    def __enter__(self, *args, **kwargs):
        sys.stdout.write(f'\rRunning operations {self.suffix}: 0')
        return self

    def __exit__(self, *args):
        sys.stdout.write('\n')
        sys.stdout.flush()

    def increment(self):
        """
        Increment counter by one.

        Represent progression in terminal.
        """
        self.counter += 1
        sys.stdout.write(f'\rRunning operations {self.suffix}: {self.counter}')
        return


class ProgressFake:
    """
    Simulate the interface of ProgressBar but does nothing.

    Used in servers where ProgressBar makes no sense.
    """

    def __init__(self, *args, **kwargs):
        return

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        return

    def increment(self):
        """
        Simulate ProgressBar interface.

        Does nothing.
        """
        return
