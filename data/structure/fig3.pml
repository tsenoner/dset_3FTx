# Reset PyMOL and delete any previously loaded objects
reinitialize everything

# Load the PDB file
load P60615.pdb, protein

# Show the protein as a cartoon
show cartoon, protein


# color exons - orange, marine, firebrick, forest, salmon
# lightorange, palecyan, salmon
color grey, protein and resi 1-95  # SP

# Select and highlight all cysteine residues
select all_cysteines, protein and resn CYS
color yellow, all_cysteines
color lime, protein and resi 49-59
color orange, protein and resi 50
color orange, protein and resi 54

# Highlight cystein bridges
show sticks, (cys/ca+cb+sg) and byres (cys/sg and bound_to cys/sg)

# Center the view on the protein
set_view (\
     0.589704275,    0.738407493,   -0.327120304,\
     0.677440941,   -0.672772110,   -0.297411382,\
    -0.439685315,   -0.046220966,   -0.896961272,\
     0.000000000,    0.000000000, -276.144134521,\
     2.804838181,   -4.053432465,   -2.784664154,\
   235.533142090,  316.755218506,  -20.000000000 )

png fig3_P60615.png, ray=1