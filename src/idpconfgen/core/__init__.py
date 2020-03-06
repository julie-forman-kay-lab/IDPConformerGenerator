"""Core modules serve lib-wide."""
import string


def has_string_formatters(s):
    """
    Determine if a string has ``{}`` operators.
    
    Parameters
    ----------
    s : str
        The string to analyze.

    Returns
    -------
    bool
        ``True`` if yes, ``False`` if no.
    """
    # see: https://stackoverflow.com/questions/46161710/
    try:
        return list(string.Formatter().parse(s))[0][1] is not None
    except IndexError:
        #  when string is ''
        return False
