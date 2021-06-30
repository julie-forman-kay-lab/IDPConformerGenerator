"""Define chunk probability distribution for CLI and building."""
from argparse import Action
from collections import namedtuple
from pathlib import Path

from libfuncpy import ITEX, chainf, is_not_none, pass_, vartial

from idpconfgen import log
from idpconfgen.libs.libcalc import make_seq_probabilities
from idpconfgen.libs.libcli import add_general_arg
from idpconfgen.libs.libio import read_lines
from idpconfgen.libs.libparse import convert_int_float_lines_to_dict
from idpconfgen.logger import S, T


_XmerProbs = namedtuple('XmerProbs', ['size', 'probs'])


def make_xmerprobs(sizes, probs):
    """Create the default representation of XmerProbs object."""
    ss = list(sizes)
    pp = make_seq_probabilities(list(probs), reverse=False)
    return _XmerProbs(ss, pp)


default_xmer_sizes = (1, 2, 3, 4, 5)
default_xmer_probs = (1, 1, 3, 3, 2)
default_XmerProbs = make_xmerprobs(default_xmer_sizes, default_xmer_probs)


def read_xmer_probs_from_file(fpath):
    """Read xmer probs from formatted file."""
    _ = chainf(fpath, read_lines, convert_int_float_lines_to_dict)
    xmerprobs = make_xmerprobs(_.keys(), _.values())
    logxmerprobs(xmerprobs)
    return xmerprobs


def logxmerprobs(xmerprobs, func=log.info):
    """Log XmerProbs object properties."""
    func(T('Chunk size and probabilities selected'))
    func(S('sizes of: {}', tuple(xmerprobs.size)))
    func(S('probs of: {}', tuple(xmerprobs.probs)))


def _give_defaults_with_log(*ignore, **everything):
    logxmerprobs(default_XmerProbs)
    return default_XmerProbs


read_xmer_probs_file_or_default = vartial(
    ITEX,
    read_xmer_probs_from_file,
    is_not_none,
    _give_defaults_with_log,
    )


class XmerProbsAction(Action):
    """Prepare Action for XmerProbs in CLI."""

    def __call__(self, parser, namespace, value, option_string=None):
        """Call on me :-)."""
        try:
            result = read_xmer_probs_file_or_default(value)
        except (ValueError, TypeError):
            parser.error(f'Not a compatible {self.dest!r} format.')
        except FileNotFoundError:
            parser.error(f'File {value} not found for parameter {self.dest!r}.')
        except IsADirectoryError:
            _ = f'{value} is a directory. We expect a file for {self.dest!r}.'
            parser.error(_)
        else:
            setattr(namespace, self.dest, result)


xmers_prob_help = """Size of peptide building chunks and their respective
selection probabilities.

Provide a too column file where the first column is the size of the chunks you
wish to allow and the second column is the relative probability to select each
choice. Both columns must be integers. Use "#" to comment lines. If you don't
provide any value, the default is chunk sizes: 1 to 5 with relative
probabilities (1, 1, 3, 3, 2). You can give any integer to define relative
probabilities: 10 and 90 will results 10/(10+90) and 90/(10+90) percent.
If a Proline residue follows the chunk being built, the additional Proline
angles will also be consider regardless of the selected chunk size."""

xmer_probs_args = ['-xp', '--xmer-probs']
xmer_probs_kwargs = {
    "help": xmers_prob_help,
    "type": Path,
    "default": None,
    "action": XmerProbsAction,
    }

add_xmer_arg = vartial(
    add_general_arg,
    *xmer_probs_args,
    **xmer_probs_kwargs,
    )

is_XmerProbs = vartial(isinstance, _XmerProbs)
prepare_xmer_probs = vartial(
    ITEX,
    pass_,
    is_XmerProbs,
    read_xmer_probs_file_or_default,
    )
