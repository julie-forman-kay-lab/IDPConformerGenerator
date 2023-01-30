"""Operations shared by client interfaces."""
import argparse
import json
import sys
from os import cpu_count

from idpconfgen import Path, __version__
from idpconfgen.core.definitions import vdW_radii_dict
from idpconfgen.libs.libio import (
    is_valid_fasta_file,
    read_FASTAS_from_file_to_strings,
    )
from idpconfgen.libs.libparse import is_valid_fasta, values_to_dict


detailed = "detailed instructions:\n\n{}"


def load_args(ap):
    """Load argparse commands."""
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
        setattr(namespace, self.dest, tuple(values))


class AllParam(argparse.Action):
    """Convert list of arguments in tuple."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        if values == ['all']:
            setattr(namespace, self.dest, 'all')

        else:
            setattr(namespace, self.dest, values)


class ListOfPositiveInts(argparse.Action):
    """
    Create list of positive integers from input string.

    Raises error if non-positive integers are given.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        try:
            iv = list(map(int, values))
        except Exception:
            raise parser.error(f'{self.dest!r} can only receive integers.')
        if any(i < 1 for i in iv):
            raise parser.error(f'{self.dest!r} can not have a negative range.')
        setattr(namespace, self.dest, iv)


class CSV2Tuple(argparse.Action):
    """Convert list of arguments in tuple."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        setattr(namespace, self.dest, tuple(values.split(',')))


class ReadDictionary(argparse.Action):
    """Read to a dictionary."""

    def __call__(self, parser, namespace, value, option_string=None):
        """Execute."""
        pvalue = Path(value)
        try:
            with pvalue.open('r') as fin:
                valuedict = json.load(fin)
        except FileNotFoundError:
            valuedict = json.loads(value)

        setattr(namespace, self.dest, valuedict)


class ListOfIntsPositiveSum(argparse.Action):
    """Convert list of str to list of ints."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Call it."""
        if values is False:
            setattr(namespace, self.dest, values)
        else:
            try:
                value_int = [int(i) for i in values]
            except ValueError:
                raise parser.error(
                    f'{option_string!r} must be a list of integer numbers.'
                    )

            if sum(value_int) < 1:
                raise parser.error(
                    f'Sum of {option_string!r} probabilities cannot '
                    'be less than 1'
                    )

            if len(value_int) > 5:
                raise parser.error(
                    f'{option_string!r} takes a maximum of 5 digits'
                    )
            setattr(namespace, self.dest, [int(i) for i in values])


class SeqOrFasta(argparse.Action):
    """Read sequence of FASTA file."""

    def __call__(self, parser, namespace, value, option_string=None):
        """Call it."""
        if is_valid_fasta(value):
            seq = value
        elif is_valid_fasta_file(value):
            seqdict = read_FASTAS_from_file_to_strings(value)
            seq = list(seqdict.values())[0]
        else:
            raise parser.error('Input sequence not valid.')

        setattr(namespace, self.dest, seq)


class ParamsToDict(argparse.Action):
    """
    Convert command-line parameters in an argument to a dictionary.

    Adapted from https://github.com/joaomcteixeira/taurenmd

    Example
    -------
    Where ``-x`` is an optional argument of the command-line client
    interface.
        >>> par1=1 par2='my name' par3=[1,2,3]
        >>> {'par1': 1, 'par2': 'my name', 'par3': [1, 2, 3]}
    """

    def __call__(self, parser, namespace, values, option_string=None):
        """Execute."""
        param_dict = values_to_dict(values)

        namespace.plotvars = param_dict
        setattr(namespace, self.dest, True)


def minimum_value(minimum):
    """Define a minimum value for action."""

    class MinimumValue(argparse.Action):
        """Definies a minimum value."""

        def __call__(self, parser, namespace, value, option_string=None):
            """Call on me :-)."""
            if value < minimum:
                parser.error(
                    f'The minimum allowed value for {self.dest} is {minimum}.'
                    )
            else:
                setattr(namespace, self.dest, value)

    return MinimumValue


def CheckExt(extensions):
    """Check extension for arguments in argument parser."""
    class Act(argparse.Action):
        """Confirms input has extension."""

        def __call__(self, parser, namespace, path, option_string=None):
            """Call on me :-)."""
            suffix = Path(path).suffix
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
        sys.stderr.write(f'\nerror: {message}\n')
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
    usage = '\n' + '\n'.join(doclines[doclines.index('USAGE:') + 1:])

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
        usage=module._usage,
        # prog=module._prog,
        description=module._prog + '\n\n' + module.ap.description,
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


# arguments index:
# positional:
# db                      : database being built
# pdb_files
# pdbids

# optional:
# -c, --chunks            : number of fragments to process in memory
# -cif, --mmcif           : prioritizes downloading mmCIF files
# -d, --destination       : destination folder
# -db, -database          : the input database
# -deg, --degrees         : degrees
# -dec, --decimals        : decimal numbers
# --no_hn_term                : has H of N-terminal
# -seq, --input_seq       : the input sequence
# -m, --minimum           : minimum size
# -n, --ncores            : number of cores
# -o, --ouput             : general output string
# -of, --output-folder    : output folder
# --replace               : replace (whatever shalls replace)
# -rd, --reduced          : reduces secondary structure representation
# -rn, --record-name      : PDB RECORD name
# -rs, --random-seed      : define a random seed
# -u, --update            : update (whatever shalls update)
# -vdWt --vdW-tolerance   : vdW tolerance (vdW_A + vdW_B - tolerance)
# -vdWr --vdW-radii       : the vdW radii set
# -vdWb --vdW-bonds-apart : Maximum bond connective to ignore in vdW validation


def add_general_arg(parser, *args, **kwargs):
    """Add a general argument with args and kwargs."""
    parser.add_argument(*args, **kwargs)


def add_argument_chunks(parser):
    """
    Add fragments argument.

    For those routines that are split into operative fragments.
    """
    parser.add_argument(
        '-c',
        '--chunks',
        help='Number of fragments to process in memory before saving to disk.',
        default=5_000,
        type=int,
        )


def add_argument_cmd(parser):
    """Add the command for the external executable."""
    parser.add_argument(
        'cmd',
        help='The path to the executable file used in this routine.',
        type=str,
        )


def add_argument_cif(parser):
    """Add command to prioritize downloading CIF parser."""
    parser.add_argument(
        '-cif',
        '--mmcif',
        help="Prioritizes PDBx/mmCIF formatted structures.",
        action='store_true',
        )


def add_argument_db(parser):
    """Add argument to store path from database."""
    parser.add_argument(
        'database',
        help='The database being built. Must be JSON.',
        type=Path,
        action=CheckExt('.json'),
        )


def add_argument_idb(parser):
    """Add argument for input database."""
    parser.add_argument(
        '-db',
        '--database',
        help='The IDPConfGen database.',
        required=True,
        )


def add_argument_dloopoff(parser):
    """Add argument `--dloop-off`."""
    parser.add_argument(
        '--dloop-off',
        help=(
            'Sampling loops is active by default. '
            'Use this flag to deactivate it.'
            ),
        action="store_true",
        )


def add_argument_dhelix(parser):
    """Add argument `--dhelix`."""
    parser.add_argument(
        '--dhelix',
        help=(
            'Samples the database also for helix segments. '
            'This feature can be used in combination with --dstrand.'
            'To explore the three secondary structures, activate --dhelix and '
            '--dstrand, loop search is always active. '
            'These features need to be used in combination with the `-rd` flag '
            'in `idpconfgen sscalc`.'
            ),
        action="store_true",
        )


def add_argument_dstrand(parser):
    """Add argument `--dstrand`."""
    parser.add_argument(
        '--dstrand',
        help=(
            'Samples the database also for strand segments. '
            'See help for `--dhelix`.'
            ),
        action="store_true",
        )


def add_argument_dany(parser):
    """Add argument `--dany`."""
    parser.add_argument(
        '--dany',
        help=(
            'Samples the database based on sequence identity only. '
            'Activating this option disregards any secondary structure '
            'annotation. Requires --dloop-off.'
            ),
        action="store_true",
        )


def add_argument_duser(parser):
    """Add argument `--duser`."""
    parser.add_argument(
        '--duser',
        help=(
            'NOTE: Very advanced users only. Use this option to define your own'
            ' regular expressions for the database sampling process. '
            'You only want to use this option if you know how the code works '
            'internally. Use this option instead of --dhelix, --dstrand, '
            '--dany. Requires --dloop-off.'
            ),
        default=None,
        nargs='+',
        )


def add_argument_seq(parser):
    """Add argument for input sequence."""
    parser.add_argument(
        '-seq',
        '--input_seq',
        help='The Conformer residue sequence. String or FASTA file.',
        required=True,
        nargs='?',
        action=SeqOrFasta,
        )


def add_argument_degrees(parser):
    """Add `degree` argument to parser."""
    parser.add_argument(
        '-deg',
        '--degrees',
        help='Whether to save angles in degrees.\nUse radians otherwise.',
        action='store_true',
        )


def add_argument_decimals(parser):
    """Add decimals to parser."""
    parser.add_argument(
        '-dec',
        '--decimals',
        help='Size of decimal numbers. Default to 5.',
        default=5,
        type=int,
        )


def add_argument_nohnterm(parser):
    """Add boolean flag to H in N-terminal."""
    parser.add_argument(
        '--no_hn_term',
        help=(
            'If given, consideres no protons exist in the N-terminal Nitrogen. '
            'If not given, considers N-terminal to have `H` or `H1` labels.'
            ),
        action='store_true',
        )


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


def add_argument_minimum(parser):
    """Add argument for minimum size of fragments."""
    parser.add_argument(
        '-m',
        '--minimum',
        help='The minimum size accepted for fragment/segment.',
        type=int,
        default=2,
        )


def add_argument_ncores(parser):
    """Add argument for number of cores to use."""
    ncpus = max(cpu_count() - 1, 1)
    parser.add_argument(
        '-n',
        '--ncores',
        help=(
            'Number of cores to use. If `-n` uses all available '
            'cores except one. To select the exact number of cores '
            'use -n #, where # is the desired number.'
            ),
        type=int,
        default=1,
        const=ncpus,
        nargs='?',
        )


def add_argument_output(parser):
    """Add argument for general output string."""
    parser.add_argument(
        '-o',
        '--output',
        help=(
            'Output file. Defaults to `None`. '
            'Read CLI instructions for `None` behaviour.'
            ),
        type=str,
        default=None,
        )


def add_argument_output_folder(parser):
    """Add argument for general output string."""
    parser.add_argument(
        '-of',
        '--output-folder',
        help=(
            "The folder where to save the output. "
            "NOTE: if new files have the same name of old files in the "
            "folder, new files will replace old files."
            ),
        type=str,
        default=None,
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
            'to \'L\', \'H\', \'G\', and \'E\'. '
            'PPII segments are considered \'L\'.'
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


def add_argument_random_seed(parser):
    """Add argument to select a random seed number."""
    parser.add_argument(
        '-rs',
        '--random-seed',
        help=(
            'Define a random seed number for reproducibility. '
            'Defaults to 0.'
            ),
        default=0,
        type=int,
        )


def add_argument_source(parser):
    """Add `source` argument to parser."""
    parser.add_argument(
        '-sc',
        '--source',
        help=(
            'Updates source with the ouput generated in this CLI. '
            'Replaces prexisting entries. '
            'If None given, saves results to a new file.'
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


# TODO: this parameters must be discontinued
def add_argument_vdWr(parser):
    """Add argument for vdW radii set selection."""
    parser.add_argument(
        '-vdWr',
        '--vdW-radii',
        help=(
            'The vdW radii values set. Defaults to: `tsai1999`.'
            ),
        default='tsai1999',
        choices=list(vdW_radii_dict.keys()),
        type=str,
        )


# TODO: this parameters must be discontinued
def add_argument_vdWt(parser):
    """Add argument for vdW tolerance."""
    parser.add_argument(
        '-vdWt',
        '--vdW-tolerance',
        help=(
            'The tolerance applicable to the vdW clash validation, '
            'in Angstroms. '
            'Tolerance follows in the formula: '
            'vdWA + vdWB - tolerance - dAB = 0. '
            'Defaults to 0.4.'
            ),
        default=0.4,
        type=float,
        )


# TODO: this parameters must be discontinued
def add_argument_vdWb(parser):
    """
    Add argument for vdW bonds apart criteria.

    Further read:

    https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/findclash/findclash.html
    """
    parser.add_argument(
        '-vdWb',
        '--vdW-bonds-apart',
        help=(
            'Maximum number of bond separation on which pairs to'
            'ignore vdW clash validation. '
            'Defaults to 3.'
            ),
        default=3,
        action=minimum_value(3),
        type=int,
        )


def add_argument_plot(parser):
    """
    Add argument for plotting parameters.

    Plot kwargs that will be passed to the plotting function.
    If given, plot results. Additional arguments can be given to
    specify the plot parameters.

    Adapted from:
    https://github.com/joaomcteixeira/taurenmd/blob/6bf4cf5f01df206e9663bd2552343fe397ae8b8f/src/taurenmd/libs/libcli.py#L539-L570

    Defined by ``--plot``.
    """
    parser.add_argument(
        '--plot',
        help=(
            'Plot torsion angles for ensembles of the same protein. '
            'Using this plot for ensembles of different proteins is devoid '
            'of meaning. '
            'Additional arguments can be given to configure the '
            'plot style. '
            "Example: --plot xlabel=Sic1_Res type=omega color='b','r'. "
            'Accepted plot arguments are defined by the plotting function used.'
            ' See ploting functions in the documentation pages. '
            'Defaults to False, no plot is produced.'
            ),
        nargs='*',
        default=False,
        action=ParamsToDict,
        )
