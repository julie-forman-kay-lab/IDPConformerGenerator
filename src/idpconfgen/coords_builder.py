"""
Client for building IDRs on PDB files in the cartesian coordinate space.

Methodology deviates from traditional IDP or beads-on-a-string FLDR/S approach.

Name: CoPS (coordinate pathway sampling)

Ideas:
1. Build a dummy string first between start and end points
   then populate with torsion angle sampling from IDPCG
- treat folded protein (surface) as boundaries of where not to go, everywhere
  else is permissable
- analyze average (or max) length of residue backbone per AA
  (that will be the vector length for each AA)
- use sequence and average lengths to build draft paths between start
  and end points of PDB (depends on random seed for reproducibility)
- For each path, sample fragments of torsion angles from database that
  match its profile (using alignment during building reduce deviation)
- maybe have variable/flag for how far away from surface of protein
- Just an idea of random-flight model? (focus on W(x,y,z;N) being at the end
  point of the connected region for folded domain)
2. Generate IDPs using OG method starting at point A and setting clash boundary
   high when it's near the surface of the protein.
- need to see where it ends, to stop building if it'll never reach the endpoint
- this can be very computationally expensive
- can start with this?
"""
