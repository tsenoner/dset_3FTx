import sys
from Bio.PDB import PDBParser
from pymol import cmd

def color_protein_sections(pdb_file, output_image):
    # Read the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    # Define the sections to color (start_residue, end_residue, color)
    sections_to_color = [
        (1, 100, 'red'),
        (101, 200, 'green'),
        (201, 300, 'blue')
    ]

    # Launch PyMOL
    cmd.reinitialize()

    # Load the PDB file into PyMOL
    cmd.load(pdb_file, 'protein')

    # Color the defined sections
    for start, end, color in sections_to_color:
        cmd.select(f'section_{color}', f'(resi {start}-{end}) and protein')
        cmd.color(color, f'section_{color}')

    # Show disulfide bridges in yellow
    cmd.select('disulfide_bridges', 'cys and sg within 2.2 of cys and sg')
    cmd.show('sticks', 'disulfide_bridges')
    cmd.color('yellow', 'disulfide_bridges')

    # Save the colored 3D structure as a PNG image
    cmd.png(output_image, ray=1)

    # Quit PyMOL
    cmd.quit()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python color_pdb.py input_pdb_file output_png_image")
        sys.exit(1)

    input_pdb_file = sys.argv[1]
    output_png_image = sys.argv[2]

    color_protein_sections(input_pdb_file, output_png_image)
