"""Operations shared by client interfaces."""
import argparse
import sys
from os import cpu_count

from idpconfgen import Path, __version__


detailed = "detailed instructions:\n\n{}"


def load_args(ap):
    cmd = ap.parse_args()
    return cmd


def maincli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


class FolderOrTar(argparse.Action):
    """Controls if input is folder, files or tar."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Hello."""
        if values[0].endswith('.tar'):
            setattr(namespace, self.dest, values[0])
        else:
            setattr(namespace, self.dest, values)


class ArgsToTuple(argparse.Action):
    """Convert list of arguments in tuple."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        namespace.record_name = tuple(values)


def CheckExt(extensions):
    """Check extension for arguments in argument parser."""
    class Act(argparse.Action):
        """Confirms input has extension."""

        def __call__(self, parser, namespace, path, option_string=None):
            """Call on me :-)."""
            suffix = path.suffix
            if suffix not in extensions:
                parser.error(f'Wrong extension {suffix}, expected {extensions}')
            else:
                setattr(namespace, self.dest, path)
    return Act


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    """Custom Parser class."""

    def error(self, message):
        """Present error message."""
        self.print_help()
        sys.stderr.write('\nerror: %s\n' % message)
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


def add_argument_chunks(parse):
    """
    Add chunks argument.

    For those routines that are split into operative chunks.
    """
    parse.add_argument(
        '-c',
        '--chunks',
        help='Number of chunks to process in memory before saving to disk.',
        default=5_000,
        type=int,
        )


# arguments index:
# positional:
# pdb_files
# pdbids

# optional:
# -d, --destination       : destination folder
# -n, --ncores            : number of cores
# --replace               : replace (whatever shalls replace)
# -rd, --reduced          : reduces secondary structure representation
# -rn, --record-name      : PDB RECORD name
# -u, --update            : update (whatever shalls update)


def add_argument_destination_folder(parser):
    """
    Add destination folder argument.

    Accepts also a TAR file path.

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
            'Alternatively, you can provide a path to a .tar file '
            'where PDBs will be saved.'
            ),
        type=Path,
        default=Path.cwd(),
        )


def add_argument_ncores(parser):
    """Add argument for number of cores to use."""
    parser.add_argument(
        '-n',
        '--ncores',
        help='Number of cores to use.',
        type=int,
        default=1,
        const=cpu_count(),
        nargs='?',
        )


def add_argument_pdb_files(parser):
    """
    Add PDBs Files entry to argument parser.

    Parameters
    ----------
    parser : `argparse.ArgumentParser` object
    """
    parser.add_argument(
        'pdb_files',
        help=(
            'Paths to PDB files in the disk. '
            'Accepts a TAR file.'
            ),
        nargs='+',
        action=FolderOrTar,
        )


def add_argument_pdbids(parser):
    """
    Add arguments for PDBIDs.

    This differs from :func:`add_parser_pdbs` that expects
    paths to files.
    """
    parser.add_argument(
        'pdbids',
        help='PDBID[CHAIN] identifiers to download.',
        nargs='+',
        )

def add_argument_reduced(parser):
    """Add `reduced` argument."""
    parser.add_argument(
        '-rd',
        '--reduced',
        help=(
            'Reduces nomenclature for secondary structure identity '
            'to \'L\', \'H\' and \'E\'.'
            ),
        action='store_true',
        )


def add_argument_replace(parser):
    """Add argument `replace`."""
    parser.add_argument(
        '--replace',
        help='Replace existent entries',
        action='store_true',
        )
    return


def add_argument_record(parser):
    """Add argument to select PDB RECORD identifier."""
    parser.add_argument(
        '-rn',
        '--record-name',
        help='The coordinate PDB record name. Default: ("ATOM", "HETATM")',
        default=('ATOM', 'HETATM'),
        action=ArgsToTuple,
        nargs='+',
        )


def add_argument_source(parser):
    """Add `source` argument to parser."""
    parser.add_argument(
        '-sc',
        '--source',
        help=(
            'Builds the output on top of the source, '
            'instead of starting from scratch. '
            'Replaces prexisting entries'
            ),
        type=Path,
        default=None,
        )


def add_argument_update(parser):
    """Add update argument to argument parser."""
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

def add_argument_cmd(parser):
    """Add the command for the external executable."""
    parser.add_argument(
        'cmd',
        help='The path to the executable file used in this routine.',
        type=str,
        )
