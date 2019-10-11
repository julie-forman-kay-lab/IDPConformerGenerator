"""Common operations for client interfaces."""

def parse_doc_params(docstring):
    """
    Parse client docstrings.

    Separates PROG, DESCRIPTION and USAGE from client main docstring.

    Parameters
    ----------
    docstring : str
        The module docstring.

    Returns
    -------
    tuple
        (prog, description, usage)
    """
    doclines = docstring.split('\n')
    prog = doclines[1]
    description = '\n'.join(doclines[3:doclines.index('USAGE:')])
    usage = '\n' + '\n'.join(doclines[doclines.index('USAGE:') + 1:])
    
    return prog, description, usage
