"""Implement residue replacement possibilities."""
import argparse

from idpconfgen.libs import libcli
from idpconfgen.libs.libparse import make_list_if_not


residue_tolerance_cli_help = \
"""Map residue replacements.

Provide a dictionary mapping which residue exchange are allowed.
Example:
-subs '{"A": "AG"}'

will consider Ala or Glycines when searching the database. When looking
for "GHLAT" both "GHLAT" and "GHLGT" will be considered valid matches."""  # noqa: E122, E501


def add_res_tolerance_groups(ap):
    """Add parameters related to residue tolerance."""
    parser = ap.add_mutually_exclusive_group(required=False)
    parser.add_argument(
        '-urestol',
        '--user-residue-tolerance',
        dest='residue_tolerance',
        help=residue_tolerance_cli_help,
        default=None,
        action=libcli.ReadDictionary,
        )

    parser.add_argument(
        '-edss50',
        '--edss50-residue-tolerance',
        dest='residue_tolerance',
        help=(
            'Uses residue tolerance matrix EDSS50. '
            'Provide the indexes of the table columns to add up tolerances. '
            'Example: -edss50 5 3 2'
            ),
        nargs='+',
        default=None,
        type=int,
        action=EDSS50_indexes,
        )


# Other interesting options for IDPs
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

_EDSSMat50_idx = {
    5: 0,
    3: 1,
    2: 2,
    1: 3,
    0: 4,
    }

EDSSMat50_subs = {
    #      5    3     2        1     0
    'A': ['',  '',   '',     'TV',  'SG'    ],  # noqa: E241, E121, E202
    'C': ['',  'Y',  'W',    'F',   'SVH'   ],  # noqa: E241, E121, E202
    'D': ['',  '',   'E',    '',    ''      ],  # noqa: E241, E121, E202
    'E': ['',  '',   'D',    '',    ''      ],  # noqa: E241, E121, E202
    'F': ['Y', '',   'WL',   'MIC', 'VH'    ],  # noqa: E241, E121, E202
    'G': ['',  '',   '',     '',    'A'     ],  # noqa: E241, E121, E202
    'H': ['',  '',   'YQ',   'N',   'FRC'   ],  # noqa: E241, E121, E202
    'I': ['',  'VM', 'L',    'F',   'T'     ],  # noqa: E241, E121, E202
    'K': ['',  '',   'R',    '',    'Q'     ],  # noqa: E241, E121, E202
    'M': ['',  'I',  'VL',   'F',   'T'     ],  # noqa: E241, E121, E202
    'N': ['',  '',   '',     'H',   'STD'   ],  # noqa: E241, E121, E202
    'P': ['',  '',   '',     '',    ''      ],  # noqa: E241, E121, E202
    'Q': ['',  '',   'H',    '',    'KR'    ],  # noqa: E241, E121, E202
    'R': ['',  '',   'K',    '',    'HQW'   ],  # noqa: E241, E121, E202
    'S': ['',  '',   '',     '',    'TCNA'  ],  # noqa: E241, E121, E202
    'T': ['',  '',   '',     'A',   'VSMIN' ],  # noqa: E241, E121, E202
    'W': ['',  '',   'YFC',  '',    'R'     ],  # noqa: E241, E121, E202
    'Y': ['F', '',   'C',    '',    'WH'    ],  # noqa: E241, E121, E202
    }
"""EDSSMat50 tolerance matrix."""


def make_EDSSMat50_subs(idx=(5, 3, 2, 1, 0)):
    """Make EDSSMat50 table com column indexes."""
    idxs = make_list_if_not(idx)
    subsd = {}
    for k, v in EDSSMat50_subs.items():
        subsd[k] = k + ''.join(v[_EDSSMat50_idx[i]] for i in idxs)
    return subsd


class EDSS50_indexes(argparse.Action):
    """Convert the number indexes to EDSS50 table."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Hello."""
        subs = make_EDSSMat50_subs(values)
        setattr(namespace, self.dest, subs)
