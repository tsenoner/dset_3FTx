# Reset PyMOL and delete any previously loaded objects
reinitialize everything

# Load the PDB file
load P0DP58.pdb, protein

# Show the protein as a cartoon
show cartoon, protein


# color exons - orange, marine, firebrick, forest, salmon
# lightorange, palecyan, salmon
color raspberry, protein and resi 1-20  # SP
color grey, protein and resi 21-91  # rest
color marine, protein and resi 92-116  # MaD
#cartoon skip, protein and resi 92-116

# Select and highlight all cysteine residues
select all_cysteines, protein and resn CYS
color yellow, all_cysteines
color orange, protein and resi 26
color orange, protein and resi 33

# Highlight cystein bridges
show sticks, (cys/ca+cb+sg) and byres (cys/sg and bound_to cys/sg)

# Center the view on the protein
set_view (\
    -0.447699487,   -0.649601221,   -0.614477873,\
    -0.722383738,   -0.142253026,    0.676701486,\
    -0.526997268,    0.746849358,   -0.405574948,\
     0.000000000,    0.000000000, -266.401733398,\
    -7.221862793,    2.549882889,    1.041866302,\
  -53763.781250000, 54296.593750000,  -20.000000000 )

png fig2_P0DP58.png, ray=1