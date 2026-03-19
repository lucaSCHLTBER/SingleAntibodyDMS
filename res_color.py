import pymol
import csv
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.pyplot as plt
import os


def normalize_scores(scores):

    min_score = min(scores)
    max_score = max(scores)
    if max_score == min_score:
        return [0.5] * len(scores)
    return [(s - min_score) / (max_score - min_score) for s in scores]

def generate_color_bar(filename="colorbar.png"):
    fig, ax = plt.subplots(figsize=(6, 0.5))
    # Gradient from grey to red
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    color_map = mcolors.LinearSegmentedColormap.from_list("custom", [(0.5, 0.5, 0.5), (1, 0, 0)])
    
    ax.imshow(gradient, aspect="auto", cmap=color_map)
    ax.set_xticks([0, 128, 256])
    ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_yticks([])
    ax.set_title("Normalized Enrichment Scores")

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()
    #print(f"Color bar saved as {filename}")

def color_rbd_from_csv(filename, object_name):
    residues = []
    scores = []
    #Reset the color before coloring the structure
    cmd.color_deep("gray50", 'RBD', 0)

    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            residues.append(int(row[0])) 
            scores.append(float(row[1]))  

    # Normalize scores for color mapping
    norm_scores = normalize_scores(scores)
    
    # Define color gradient (grey to red)
    color_map = mcolors.LinearSegmentedColormap.from_list("custom", [(0.5, 0.5, 0.5), (1, 0, 0)])
    
    for res, norm_score in zip(residues, norm_scores):
        color = color_map(norm_score)[:3]  
        pymol_color = f"[{color[0]},{color[1]},{color[2]}]"
        cmd.set_color(f"color_{res}", color)
        cmd.color(f"color_{res}", f"{object_name} and resi {res}")
    
    generate_color_bar("colorbar.png")


    filename_only = os.path.basename(filename)  
    filename_only = os.path.splitext(filename_only)[0]  
    output_path = os.path.join(r'C:\Users\au649453\OneDrive - Aarhus universitet\PhD\2 Projects\Deep mutational scanning of SarsCov2\Pymol\heatmaps_barcode', filename_only + '.png')

    # Save the PyMOL image correctly
    cmd.png(output_path, dpi=300, ray=1, width=1200, height=900, quiet=0)

#Defining single-letter aminoacids
aa_three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}     

def label_residues(chain,residues):
    for res in residues:
        resn_list = []
        cmd.iterate(f"chain {chain} and resi {res} and name CA", "resn_list.append(resn)", space=locals())
        one_letter_code = aa_three_to_one.get(resn_list[0], "?")  # Convert to one-letter code
        cmd.label(f"chain {chain} and resi {res} and name CA", f'"{one_letter_code}{res}"')

def label_surface_atoms(chain, residues): #This can help to better hightlight a label 
    for res in residues:
        model = cmd.get_model(f"chain {chain} and resi {res}")
        if model.atom:
            # Get the atom farthest from the molecular center (surface-exposed atom)
            most_superficial_atom = max(model.atom, key=lambda a: (a.coord[0]**2 + a.coord[1]**2 + a.coord[2]**2))
            one_letter_code = aa_three_to_one.get(most_superficial_atom.resn, "?")

            cmd.label(f"chain {chain} and resi {res} and name {most_superficial_atom.name}", f'"{one_letter_code}{res}"')



#Setting up the PDB file for consistency

cmd.fetch("6M0J") 
cmd.remove('not alt ""') #Remove all alternative conformations
cmd.remove('solvent')
cmd.create("RBD", 'Chain E')
cmd.create('ACE2', 'chain A')
cmd.delete('6m0j')

cmd.hide("everything")

#cmd.show("cartoon", "RBD")
cmd.show("surface", "RBD") 

cmd.set_view((\
    -0.966984689,   -0.237655953,    0.091937542,\
    -0.165207669,    0.310039431,   -0.936253130,\
     0.194012970,   -0.920531809,   -0.339055121,\
     0.000000000,    0.000000000, -218.305236816,\
   -32.414237976,   26.077138901,   21.105247498,\
   184.219543457,  252.390991211,  -20.000000000 ))

#Change the label settings
cmd.set("label_color", "yellow")
cmd.set("label_outline_color", "black")
cmd.set("label_size", 30)
cmd.set("label_position", (0, 0, 10))
cmd.set("depth_cue", 0)
cmd.set("ray_trace_fog", 0)

print("PDB file loaded and ready for labeling using label_residues(chain, residues)")
print("Use color_rbd_from_csv(filename, object_name) to color the RBD structure based enrichment ratios from a csv-file.")



