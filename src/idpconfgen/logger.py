"""Manages operations with logging."""
import logging
from functools import wraps
from inspect import signature
from pathlib import Path
from time import time_ns
from contextlib import contextmanager


from idpconfgen import log


def titlelog(msg, *args):
    """Format a message to a title."""
    msg = msg.title()
    return f'\n* {msg.format(*args)} ...'


def subline(msg, *args, spacer=' ', indent=4):
    """Format a line to a sub log line."""
    return f'{spacer * indent}{msg.format(*args)}'


T = titlelog
S = subline


def init_files(log, logfilesname):
    """
    Initiate log files.

    Three log files are created:
        - .debug
        - .log
        - .error

    where, .debug stores all technical details related to debugging
    routines; .log referes to user directed messages, this is the log
    the user wants to save for future reference; finally, .error
    stores errors in runtime that compromise the scientific output.
    """
    debugfile = logging.FileHandler(f'{logfilesname}.debug', mode='w')
    debugfile.setLevel(logging.DEBUG)
    debugfile.setFormatter(logging.Formatter(
        "%(filename)s:%(name)s:%(funcName)s:%(lineno)d: %(message)s"
        ))
    log.addHandler(debugfile)

    infolog = logging.FileHandler(f'{logfilesname}.log', mode='w')
    infolog.setLevel(logging.INFO)
    log.addHandler(infolog)
    infolog.setFormatter(logging.Formatter('%(message)s'))

    errorlog = logging.FileHandler(f'{logfilesname}.error', mode='w')
    errorlog.setLevel(logging.ERROR)
    errorlog.setFormatter(logging.Formatter('%(message)s'))
    log.addHandler(errorlog)


#def report_on_break(exception, folder=Path.cwd(), ext='rpr_on_break'):
#    """
#    Report to a file upon exception.
#
#    This is a decorator function.
#
#    If `func` raises `exception` a detailed report of the `func`,
#    passed `args` and `kwargs` and other details are saved to a file.
#    File name is given by `time.time`. File saves in the CWD.
#
#    Propagates the `exception`.
#
#    Parameters
#    ----------
#    exception : Exception or tuple of Exceptions
#        The Exception type(s) to consider.
#
#    ext : str
#        The extension with which save the report file.
#    """
#    def decorator(func):
#        #@wraps(func)
#        def wrapper(*args, **kwargs):
#            try:
#                result = func(*args, **kwargs)
#            except exception as err:
#                sig = signature(func)
#                s = (
#                    '#Recording error for exception: {}\n\n'
#                    '#function: {}\n\n' #                    '#signature: {}\n\n'
#                    '#Args:\n\n{}\n\n'
#                    '#Kwargs:\n\n{}\n\n'
#                    ).format(
#                        exception,
#                        func.__qualname__,
#                        sig,
#                        '\n\n'.join(map(str, args)),
#                        '\n\n'.join(map(str, kwargs.items())),
#                        )
#                fout_path = Path(folder, f'{hash(s)}_{time_ns()}.{ext}')
#                fout_path.write_text(s)
#                log.error(S('saved ERROR REPORT: {}', fout_path))
#                raise err
#            return result
#        return wrapper
#    return decorator
def report_on_break(func, *args, exception=Exception, ROB_folder=Path.cwd(), ROB_ext='rpr_on_break', **kwargs):
    """
    Report to a file upon exception.

    This is a decorator function.

    If `func` raises `exception` a detailed report of the `func`,
    passed `args` and `kwargs` and other details are saved to a file.
    File name is given by `time.time`. File saves in the CWD.

    Propagates the `exception`.

    Parameters
    ----------
    exception : Exception or tuple of Exceptions
        The Exception type(s) to consider.

    ext : str
        The extension with which save the report file.
    """
    try:
        return func(*args, **kwargs)
        #yield
    except exception as err:
        sig = signature(func)
        s = (
            '#Recording error for exception: {}\n\n'
            '#function: {}\n\n'
            '#signature: {}\n\n'
            '#Args:\n\n{}\n\n'
            '#Kwargs:\n\n{}\n\n'
            ).format(
                exception,
                func.__qualname__,
                sig,
                '\n\n'.join(map(str, args)),
                '\n\n'.join(map(str, kwargs.items())),
                )
        fout_path = Path(ROB_folder, f'{hash(s)}_{time_ns()}.{ROB_ext}')
        fout_path.write_text(s)
        log.error(S('saved ERROR REPORT: {}', fout_path))
        raise err
