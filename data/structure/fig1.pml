# Reset PyMOL and delete any previously loaded objects
reinitialize everything

# Load the PDB file
load P0DP58.pdb, protein

# Show the protein as a cartoon
show cartoon, protein

# color exons - orange, marine, firebrick, forest, salmon
# lightorange, palecyan, salmon
color lightorange, protein and resi 1-17  # exon 1
color palecyan, protein and resi 18-52  # exon 2
color salmon, protein and resi 53-116  # exon 3

# Center the view on the protein
set_view (\
    -0.447699487,   -0.649601221,   -0.614477873,\
    -0.722383738,   -0.142253026,    0.676701486,\
    -0.526997268,    0.746849358,   -0.405574948,\
     0.000000000,    0.000000000, -266.401733398,\
    -7.221862793,    2.549882889,    1.041866302,\
  -53763.781250000, 54296.593750000,  -20.000000000 )

png fig1_P0DP58.png, ray=1