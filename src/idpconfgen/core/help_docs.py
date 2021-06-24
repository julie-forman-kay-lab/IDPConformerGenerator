"""
Help documentation.

Long strings and related.
"""
residue_substitutions_cli_help = \
"""Map residue replacements.

Provide a dictionary mapping which residue exchange are allowed.
Example:
-subs '{"A": "AG"}'

will consider Ala or Glycines when searching the database. When looking
for "GHLAT" both "GHLAT" and "GHLGT" will be considered valid matches."""

xmers_prob_help = \
"""Relative probability to select the building chunk size.

Chunksize can range from 1 to 5 (inclusive).
--xmers-probs is a list of five numbers corresponding to chunk sizes
from 1 to 5. For example, 1 0 1 0 0 will create IDP conformers from
chunks of size 1 (single residue) or 3 (tripeptide segments) and both
will have the same selection probability. 0 1 2 0 0, will create IDP
conformers from chunks of 2 and 3 residues where the later has twice
probability to be selected. If none is given, it defauls to 1 2 3 4 5.
You can also normalize relative probilities to 100, for example: 10 20
30 30 10. Numbers should be integers, do not give fractional numbers. If
a Proline residue follows the chunk being built, the additional Proline
angles will also be consider regardless of the selected chunk size.
Missing digits are filled with 0, that is, -xp 1 1 is equal to -xp 1 1 0
0 0.
"""

forced_helical_help = \
"""Forces helical secondary propensity in specific regions. Provide a
string of the size of `-seq` with characters from 0-9 and *, where 0 is
0 percent, 1 comprises 10-19 perc., 9 is 90-99 perc., and * is 100 perc.
Regions of `0 perc.` will be populated with the definitions from `-dr`
flag. If `-dr` also comprises helices, these will be taken into account
as well. You can provide a string in the argument or a path to a
FASTA-file like.  Defaults to FALSE, no forced propensity will be
used."""

forced_loop_help = \
"""Forces loop propensity in specific regions. Provide a string of the
size of `-seq` with characters from 0-9 and *, where 0 is 0 percent 1
comprises 10-19 perc. , 9 is 90-99 perc. , and * is 100 perc. . Regions
of `0 perc. ` will be populated with the definitions from `-dr` flag. If
`-dr` also comprises loops (very likely), these will be taken into
account as well. You can provide a string in the argument or a path to a
FASTA-file like.  Defaults to FALSE, no forced propensity will be
used."""
