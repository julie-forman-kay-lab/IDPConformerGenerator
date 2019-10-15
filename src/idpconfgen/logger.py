"""Manages operations with logging."""
import logging


def init_files(log, logfilesname):
    """Initiates log files.""" 
    debugfile = logging.FileHandler(f'{logfilesname}.debug', mode='w')
    debugfile.setLevel(logging.DEBUG)
    debugfile.setFormatter(
        "%(filename)s:%(name)s:%(funcName)s:%(lineno)d: %(message)s"
        )
    log.addHandler(debugfile)

    infolog = logging.FileHandler(f'{logfilesname}.log', mode='w')
    infolog.setLevel(logging.INFO)
    infolog.setFormatter('%(message)s')
    log.addHandler(infolog)

    errorlog = logging.FileHandler(f'{logfilesname}.error', mode='w')
    errorlog.setLevel(logging.ERROR)
    errorlog.setFormatter('%(message)s')
    log.addHandler(errorlog)
