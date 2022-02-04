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

dssp_ppii_assignment_cli_help = \
    """
    Adds polyproline type-2 helix assignment to DSSP
    
    USAGE
    $0 <PDB_File> -> Read PDB_File and print DSSP output (with PPII helix) on stdout
    
    OPTIONS
    -horiz : Give output DSSP in 1D sequence fashion 
    """