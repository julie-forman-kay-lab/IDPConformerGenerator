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

