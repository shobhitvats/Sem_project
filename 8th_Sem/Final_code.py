import streamlit as st
import pandas as pd
from Bio.PDB import PDBParser
import tempfile
import numpy as np
import math
import matplotlib.pyplot as plt

# Streamlit app
st.set_page_config(page_title="Chalcogen Interaction Detector", layout="wide", initial_sidebar_state="expanded")



# Constants for chalcogen bond detection
VDW_S = 1.80  # Van der Waals radius of sulfur (S) in Å
VDW_O = 1.52  # Van der Waals radius of oxygen (O) in Å
VDW_N = 1.55  # Van der Waals radius of nitrogen (N) in Å
DISTANCE_THRESHOLD_SO = 3.6  # S···O threshold in Å
DISTANCE_THRESHOLD_SN = 3.6  # S···N threshold in Å
ANGLE_THETA_MIN = 115  # Minimum θ in degrees
ANGLE_THETA_MAX = 155  # Maximum θ in degrees
ANGLE_DELTA_MIN = -50  # Minimum δ in degrees
ANGLE_DELTA_MAX = 50   # Maximum δ in degrees

# Helper functions
def calculate_distance(atom1, atom2):
    coord1 = atom1.coord
    coord2 = atom2.coord
    return math.sqrt((coord1[0] - coord2[0])**2 +
                     (coord1[1] - coord2[1])**2 +
                     (coord1[2] - coord2[2])**2)

def calculate_centroid(bonded_atoms):
    coords = [atom.coord for atom in bonded_atoms]
    return np.mean(coords, axis=0)

def calculate_theta(sulfur, bonded_atoms, acceptor):
    centroid = calculate_centroid(bonded_atoms)
    s_coord = np.array(sulfur.coord)
    centroid = np.array(centroid)
    acceptor_coord = np.array(acceptor.coord)
    vec_sc = centroid - s_coord
    vec_sa = acceptor_coord - s_coord
    vec_sc_norm = vec_sc / np.linalg.norm(vec_sc)
    vec_sa_norm = vec_sa / np.linalg.norm(vec_sa)
    dot_product = np.dot(vec_sc_norm, vec_sa_norm)
    dot_product = np.clip(dot_product, -1.0, 1.0)
    return np.degrees(np.arccos(dot_product))

def calculate_delta(x_atom, bonded_atoms, sulfur, acceptor):
    centroid = calculate_centroid(bonded_atoms)
    b1 = centroid - x_atom.coord
    b2 = sulfur.coord - centroid
    b3 = acceptor.coord - sulfur.coord
    b1 = b1 / np.linalg.norm(b1)
    b2 = b2 / np.linalg.norm(b2)
    b3 = b3 / np.linalg.norm(b3)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    cos_delta = np.dot(n1, n2)
    cos_delta = np.clip(cos_delta, -1.0, 1.0)
    delta = np.degrees(np.arccos(cos_delta))
    if np.dot(np.cross(n1, n2), b2) < 0:
        delta = -delta
    return delta


def is_ligand(residue):
    """
    Check if a residue is considered a ligand.
    Ligands are defined as entities that are not part of the main backbone chain.
    Even standard amino acids can be classified as ligands if they are not part of the main chain.
    """
    # List of standard amino acid three-letter codes
    standard_amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
        "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
        "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    ]

    # Check if the residue is part of the main chain
    if "CA" in [atom.name for atom in residue]:  # Check for alpha carbon (CA) in the residue
        return False  # Residue is part of the main chain

    # If not part of the main chain, classify as a ligand
    return True

def detect_chalcogen_bonds_with_ligands(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    protein_protein_bonds = []
    protein_ligand_bonds = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'S':  # Check if the atom is sulfur
                        sulfur = atom
                        # Find two carbon atoms bonded to sulfur
                        bonded_atoms = [nbr for nbr in residue if nbr.element == 'C' and calculate_distance(sulfur, nbr) < 1.9]
                        if len(bonded_atoms) < 2:
                            continue
                        # Include sulfur in the centroid calculation
                        atoms_for_centroid = bonded_atoms + [sulfur]
                        x_atom = bonded_atoms[1]  # Reference atom for delta calculation
                        for other_model in structure:
                            for other_chain in other_model:
                                for other_residue in other_chain:
                                    for other_atom in other_residue:
                                        if other_atom.element in ['O', 'N']:
                                            distance = calculate_distance(sulfur, other_atom)
                                            if ((other_atom.element == 'O' and distance <= DISTANCE_THRESHOLD_SO) or
                                                (other_atom.element == 'N' and distance <= DISTANCE_THRESHOLD_SN)):
                                                theta = calculate_theta(sulfur, atoms_for_centroid, other_atom)
                                                delta = calculate_delta(x_atom, atoms_for_centroid, sulfur, other_atom)
                                                if (ANGLE_THETA_MIN <= theta <= ANGLE_THETA_MAX) and \
                                                   (ANGLE_DELTA_MIN <= delta <= ANGLE_DELTA_MAX):
                                                    # Simplify the description of the centroid atoms
                                                    centroid_atoms_details = [atom.name for atom in atoms_for_centroid]
                                                    bond = {
                                                        "Sulfur Residue": f"{residue.resname} {residue.id[1]} (Chain {residue.get_parent().id})",
                                                        "Acceptor Residue": f"{other_residue.resname} {other_residue.id[1]} (Chain {other_residue.get_parent().id})",
                                                        "Sulfur Coord": f"({sulfur.coord[0]:.2f}, {sulfur.coord[1]:.2f}, {sulfur.coord[2]:.2f})",
                                                        "Acceptor Coord": f"({other_atom.coord[0]:.2f}, {other_atom.coord[1]:.2f}, {other_atom.coord[2]:.2f})",
                                                        "Centroid Coord": f"({calculate_centroid(atoms_for_centroid)[0]:.2f}, {calculate_centroid(atoms_for_centroid)[1]:.2f}, {calculate_centroid(atoms_for_centroid)[2]:.2f})",
                                                        "Centroid Atoms": ", ".join(centroid_atoms_details),  # Only atom names
                                                        "Theta (θ)": theta,
                                                        "Delta (δ)": delta,
                                                        "Distance (Å)": f"{distance:.2f}",
                                                        "Theta Atoms": f"Sulfur ({sulfur.name}) → Centroid → Acceptor ({other_atom.name})",
                                                        "Delta Atoms": f"X ({x_atom.name}) → Centroid → Sulfur ({sulfur.name}) → Acceptor ({other_atom.name})"
                                                    }
                                                    if is_ligand(residue) or is_ligand(other_residue):
                                                        protein_ligand_bonds.append(bond)
                                                    else:
                                                        protein_protein_bonds.append(bond)
    return protein_protein_bonds, protein_ligand_bonds

# Custom header
st.markdown(
    """
    <style>
    .header {
        font-size: 36px;
        font-weight: bold;
        color: #4CAF50;
        text-align: center;
        margin-bottom: 10px;
    }
    .subheader {
        font-size: 18px;
        color: #555555;
        text-align: center;
        margin-bottom: 20px;
    }
    .subsubheader {
        font-size: 15px;
        color: #555555;
        text-align: center;
        margin-bottom: 20px;
    }
    </style>
    <div class="subsubheader">Developed by Shobhit Vats | SK Lab</div>
    <div class="header">Chalcogen Interaction Detector</div>
    <div class="subheader">Upload a PDB file to detect and visualize chalcogen bonds in proteins and ligands</div>
    """,
    unsafe_allow_html=True,
)

# File uploader
st.markdown("### 📂 Upload Your PDB File")
uploaded_file = st.file_uploader("Drag and drop or browse your PDB file", type=["pdb"])

st.markdown("""
### Atom Labels and Their Meanings
- **CD**: Carbon Delta (e.g., third carbon in the side chain of certain amino acids like ARG, GLU).
- **CE**: Carbon Epsilon (e.g., terminal carbon in the side chain of MET, LYS).
- **SD**: Sulfur Delta (e.g., sulfur atom in the side chain of MET).
- **SG**: Sulfur Gamma (e.g., sulfur atom in the side chain of CYS).
""")

if uploaded_file is not None:
    # Save the uploaded file to a temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_file:
        temp_file.write(uploaded_file.read())
        temp_pdb_path = temp_file.name

    # Detect chalcogen bonds
    st.write("Processing the uploaded PDB file...")
    pp_bonds, pl_bonds = detect_chalcogen_bonds_with_ligands(temp_pdb_path)

       # Display protein-protein bonds
    st.subheader("Protein-Protein Chalcogen Bonds 🧬")
    if pp_bonds:
        # Create a DataFrame and set the index to start from 1
        pp_df = pd.DataFrame(pp_bonds)
        pp_df.index = range(1, len(pp_df) + 1)  # Set index to start from 1
        st.dataframe(pp_df)
    else:
        st.write("No protein-protein chalcogen bonds found.")
    
    # Display protein-ligand bonds
    st.subheader("Protein-Ligand Chalcogen Bonds 🔗")
    if pl_bonds:
        # Create a DataFrame and set the index to start from 1
        pl_df = pd.DataFrame(pl_bonds)
        pl_df.index = range(1, len(pl_df) + 1)  # Set index to start from 1
        st.dataframe(pl_df)
    else:
        st.write("No protein-ligand chalcogen bonds found.")

    # Combine all interactions for the graph
    all_bonds = pp_bonds + pl_bonds

    if all_bonds:
        # Extract theta and delta values
        theta_values = [bond["Theta (θ)"] for bond in all_bonds]
        delta_values = [bond["Delta (δ)"] for bond in all_bonds]

        # Create a scatter plot
        fig, ax = plt.subplots(figsize=(10, 6))
        scatter = ax.scatter(
            delta_values, theta_values, c='blue', alpha=0.8, edgecolor='black', s=120
        )

        # Label each point with its interaction number
        for i, (delta, theta) in enumerate(zip(delta_values, theta_values)):
            ax.text(delta, theta, str(i + 1), fontsize=10, ha='right', color='darkred')

        # Add labels, title, and grid
        ax.set_xlabel("Delta (δ) [degrees]", fontsize=14, labelpad=10)  # X-axis is now delta
        ax.set_ylabel("Theta (θ) [degrees]", fontsize=14, labelpad=10)  # Y-axis is now theta
        ax.set_title("Delta vs Theta for Chalcogen Bond Interactions", fontsize=16, pad=15)
        ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.7)

        # Set axis limits to the specified ranges
        ax.set_xlim(-90, 90)  # Delta axis from -90 to 90
        ax.set_ylim(0, 180)  # Theta axis from 0 to 180

        # Add a rectangle to mark the boundaries
        rect = plt.Rectangle(
            (ANGLE_DELTA_MIN, ANGLE_THETA_MIN),  # Bottom-left corner
            ANGLE_DELTA_MAX - ANGLE_DELTA_MIN,  # Width
            ANGLE_THETA_MAX - ANGLE_THETA_MIN,  # Height
            linewidth=2,
            edgecolor='red',
            facecolor='none',
            linestyle='--'
        )
        ax.add_patch(rect)

        # Add a legend
        ax.legend(["Chalcogen Bond Interactions"], loc='upper right', fontsize=12)

        # Customize ticks
        ax.tick_params(axis='both', which='major', labelsize=12)

        # Display the plot in Streamlit
        st.subheader("Delta vs Theta Graph 📊")
        st.pyplot(fig)
    else:
        st.write("No chalcogen bonds found to plot.")

# Footer
st.markdown(
    """
    <style>
    .footer {
        font-size: 14px;
        color: #888888;
        text-align: center;
        margin-top: 50px;
    }
    </style>
    <div class="footer">
        Developed by Shobhit Vats | Powered by Streamlit
    </div>
    """,
    unsafe_allow_html=True,
)
