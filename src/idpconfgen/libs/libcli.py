"""Common operations for client interfaces."""
import argparse
import sys


from idpconfgen import __version__, Path


detailed = "detailed instructions:\n\n{}"


class ArgsToTuple(argparse.Action):
    """Convert list of arguments in tuple."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        namespace.record_name = tuple(values)


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    """Custom Parser class."""

    def error(self, message):
        """Present error message."""
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
    doclines = docstring.lstrip().split('\n')
    prog = doclines[0]
    description = '\n'.join(doclines[2:doclines.index('USAGE:')])
    usage = '\n'.join(doclines[doclines.index('USAGE:') + 1:])

    return prog, description, usage


def add_subparser(parser, module):
    """
    Add a subcommand to a parser.

    Parameters
    ----------
    parser : `argparse.add_suparsers object <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_subparsers>`_
        The parser to add the subcommand to.

    module
        A python module containing the characteristics of a taurenmd
        client interface. Client interface modules require the following
        attributes: ``__doc__`` which feeds the `description argument <https://docs.python.org/3/library/argparse.html#description>`_
        of `add_parser <https://docs.python.org/3/library/argparse.html#other-utilities>`_,
        ``_help`` which feeds `help <https://docs.python.org/3/library/argparse.html#help>`_,
        ``ap`` which is an `ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_,
        and a ``main`` function, which executes the main logic of the interface.
    """  # noqa: E501
    new_ap = parser.add_parser(
        module._name,
        description=module.ap.description,
        help=module._help,
        parents=[module.ap],
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    new_ap.set_defaults(func=module.main)


def add_version(parser):
    """
    Add version ``-v`` option to parser.

    Displays a message informing the current version.
    Also accessible via ``--version``.

    Parameters
    ----------
    parser : `argparse.ArgumentParser <https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser>`_
        The argument parser to add the version argument.
    """  # noqa: E501
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        # the _BANNER contains information on the version number
        version=__version__,
        )

def add_argument_ncores(parser):
    parser.add_argument(
        '-n',
        '--ncores',
        help='Number of cores to use.',
        type=int,
        default=1,
        )


def add_parser_pdbs(parser):
    """
    Add PDBs entry to argument parser.

    Parameters
    ----------
    parser : `argparse.ArgumentParser` object
    """
    parser.add_argument(
        'pdbs',
        help='PDB file list.',
        nargs='+',
        )

def add_parser_destination_folder(parser):
    """
    Adds destination folder argument.

    Parameters
    ----------
    parser : `argparse.ArgumentParser` object
    """
    parser.add_argument(
        '-d',
        '--destination',
        help=(
            'Destination folder where PDB files will be stored. '
            'Defaults to current working directory.'
            ),
        type=Path,
        default=Path.cwd(),
        )

def add_argument_update(parser):
    parser.add_argument(
        '-u',
        '--update',
        help=(
            'Updates destination folder according to input PDB list. '
            'If not provided a comparison between input and destination '
            'is provided.'
            ),
        action='store_true',
        )

