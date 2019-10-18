"""Manages operations with logging."""
import logging


class TitleLog:
    """Format string to title."""
    def __init__(self, msg):
        self.msg = msg.title()
    
    def __str__(self):
        return '\n* {} ...'.format(self.msg)


class SubLog:
    """
    Format string to bullet point like structure.
    
    This format performs nicely under the `TitleLog` formatting.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return '    {}'.format(self.msg)


T = TitleLog
S = SubLog


def init_files(log, logfilesname):
    """Initiates log files."""
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
