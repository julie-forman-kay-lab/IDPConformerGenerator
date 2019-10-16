"""Common operations for client interfaces."""
import argparse
import sys


detailed = "detailed instructions:\n\n{}"


class ArgsToTuple(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.record_name = tuple(values)


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


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
