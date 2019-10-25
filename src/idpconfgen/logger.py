"""Manages operations with logging."""
import logging


class TitleLog:
    """Format string to title."""

    def __init__(self, msg):
        self.msg = msg.title()
    
    def __str__(self):
        """Represent object as string."""
        return '\n* {} ...'.format(self.msg)


class SubLog:
    """
    Format string to bullet point like structure.
    
    This format performs nicely under the :class:`TitleLog` formatting.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        """Represent object as string."""
        return '    {}'.format(self.msg)


T = TitleLog
S = SubLog


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
