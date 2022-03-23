"""Implement residue replacement possibilities."""
from idpconfgen.libs.libparse import make_list_if_not
from idpconfgen.libs import libcli


residue_substitutions_cli_help = \
"""Map residue replacements.

Provide a dictionary mapping which residue exchange are allowed.
Example:
-subs '{"A": "AG"}'

will consider Ala or Glycines when searching the database. When looking
for "GHLAT" both "GHLAT" and "GHLGT" will be considered valid matches."""


def add_substitution_groups(ap):
    parser = ap.add_mutually_exclusive_group(required=False)
    parser.add_argument(
        '-usubs',
        '--user-residue-substitutions',
        help=residue_substitutions_cli_help,
        default=None,
        action=libcli.ReadDictionary,
        )

    parser.add_argument(
        '-edss50',
        '--edss50-residue-substitutions',
        nargs='+',
        default=(5, 3, 2),
        help=(
            'Uses residue substituion matrix EDSS50. '
            'Provide the indexes of the table columns to add up replacements. '
            'Example: -edss50 5 3 2'
            ),
        )


idp1 = {
    "Y": "YF",
    "F": "FY",
    "K": "KR",
    "R": "RK",
    "V": "VILA",
    "A": "AVIL",
    "I": "ILVA",
    "L": "LVIA",
    "Q": "QNST",
    "T": "TS",
    "S": "ST",
    }


EDSSMat50_idx = {5:0, 3:1, 2:2, 1:3, 0:4}

EDSSMat50_subs = {
  #   5    3     2        1     0
'A': ['',  '',   '',     'TV',  'SG'],
'C': ['',  'Y',  'W',    'F',    'SVH'],
'D': ['',  '',   'E',    '',    ''],
'E': ['',  '',   'D',    '',    ''],
'F': ['Y', '',   'WL',   'MIC', 'VH'],
'G': ['',  '',   '',     '',    'A'],
'H': ['',  '',   'YQ',   'N',   'FRC'],
'I': ['',  'VM', 'L',    'F',   'T'],
'K': ['',  '',   'R',    '',    'Q'],
'M': ['',  'I',  'VL',   'F',   'T'],
'N': ['',  '',   '',     'H',   'STD'],
'P': ['',  '',   '',     '',    ''],
'Q': ['',  '',   'H',    '',    'KR'],
'R': ['',  '',   'K',    '',    'HQW'],
'S': ['',  '',   '',     '',    'TCNA'],
'T': ['',  '',   '',     'A',   'VSMIN'],
'W': ['',  '',   'YFC',  '',    'R'],
'Y': ['F', '',   'C',    '',    'WH'],
    }


def make_EDSSMat50_subs(idx=(5, 3, 2, 1, 0)):
    """Make subs."""
    idxs = make_list_if_not(idx)
    subsd = {}
    for k, v in subs.items():
        subsd[k] = k + ''.join(v[EDSSMat50_idx[i]] for i in idxs)
    return subsd


