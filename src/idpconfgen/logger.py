"""Manages operations with logging."""
import logging
import traceback
from functools import partial
from inspect import signature
from pathlib import Path
from time import time_ns

from idpconfgen import log
from idpconfgen.core.exceptions import ReportOnCrashError


def titlelog(msg, *args):
    """Format a message to a title."""
    msg = msg.title()
    return f'{msg.format(*args)}:'


def subline(msg, *args, spacer=' ', indent=4):
    """Format a line to a sub log line."""
    return f'{spacer * indent}{msg.format(*args)}'


T = titlelog
S = subline
Snull = partial(subline, spacer='', indent=0)


def init_files(log, logfilesname, clear=False):
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
    if clear:
        log.handlers.clear()

    debugfile = logging.FileHandler(f'{logfilesname}.debug', mode='w')
    debugfile.setLevel(logging.DEBUG)
    debugfile.setFormatter(logging.Formatter(
        "[%(asctime)s]%(filename)s:%(name)s:%(funcName)s:%(lineno)d: "
        "%(message)s"
        ))
    log.addHandler(debugfile)

    infolog = logging.FileHandler(f'{logfilesname}.log', mode='w')
    infolog.setLevel(logging.INFO)
    log.addHandler(infolog)
    infolog.setFormatter(logging.Formatter('[%(asctime)s]%(message)s'))

    errorlog = logging.FileHandler(f'{logfilesname}.error', mode='w')
    errorlog.setLevel(logging.ERROR)
    errorlog.setFormatter(logging.Formatter('[%(asctime)s]%(message)s'))
    log.addHandler(errorlog)


init_clean_files = partial(init_files, clear=True)


def report_on_crash(
        func,
        *args,
        ROC_exception=Exception,
        ROC_ext='rpr_on_crash',
        ROC_folder=None,
        ROC_prefix='ROC',
        **kwargs,
        ):
    """
    Report to a file upon exception(s).

    If `func` raises `exception`(s), a detailed report of the `func` and
    passed `args` and `kwargs` and other details are saved to a file.
    File name is given by:
        * prefix
        * hash value of the text to save
        *`time.time_ns`.

    File saves to the CWD, unless otherwise stated.

    Propagates the `exception`.

    Keyword parameters have `ROC` verbose names to avoid collision with
    `func`'s `kwargs`.

    Parameters
    ----------
    ROC_exception : Exception or tuple of Exceptions
        The Exception type(s) to consider.

    ROC_ext : str
        The extension with which save the report file.

    ROC_folder : str or Path,
        The folder where to save the report.
        Defaults to the CWD.

    ROC_prefix : str
        The prefix of the report file.
    """
    ROC_folder = ROC_folder or Path.cwd()

    try:
        return func(*args, **kwargs)

    except ROC_exception as err:
        sig = signature(func)
        s = (
            '#Recording error for exception: {}\n\n'
            '#function: {}\n\n'
            '#signature: {}\n\n'
            '#Args:\n\n{}\n\n'
            '#Kwargs:\n\n{}\n\n'
            '#traceback:\n\n{}\n\n'
            ).format(
                ROC_exception,
                str(func),
                sig,
                '\n\n'.join(map(str, args)),
                '\n\n'.join(map(str, kwargs.items())),
                traceback.format_exc(),
                )

        fout_path = Path(
            ROC_folder,
            f'{ROC_prefix}_{hash(s)}_{time_ns()}.{ROC_ext}',
            )
        fout_path.write_text(s)

        log.error(S('saved ERROR REPORT: {}', fout_path))

        roc_error = ReportOnCrashError(fout_path)
        raise roc_error from err


def pre_msg(msg, sep=']'):
    def func(logmsg):
        return f'{msg}{sep}{logmsg}'
    return func
