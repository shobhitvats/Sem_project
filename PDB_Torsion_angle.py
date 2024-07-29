import numpy as np
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.vectors import calc_dihedral
import large_pdb_list


def calculate_phi_psi(structure):
    phi_psi_angles = []

    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        for i, residue in enumerate(pp):
            if i == 0 or i == len(pp) - 1:
                continue

            # Calculate Phi
            prev_residue = pp[i - 1]
            n = residue["N"].get_vector()
            ca = residue["CA"].get_vector()
            c = residue["C"].get_vector()
            prev_c = prev_residue["C"].get_vector()
            phi = calc_dihedral(prev_c, n, ca, c)

            # Calculate psi
            next_residue = pp[i + 1]
            n = residue["N"].get_vector()
            ca = residue["CA"].get_vector()
            c = residue["C"].get_vector()
            next_n = next_residue["N"].get_vector()
            psi = calc_dihedral(n, ca, c, next_n)

            phi_psi_angles.append([phi, psi])
            

    return phi_psi_angles


for protein_id in large_pdb_list.filtered_pdb_list:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", "./Large_Data/filtered_pdb/" + protein_id + ".pdb")
    try:
        np.savetxt("./Large_Data/filtered_angles/" + protein_id + ".angles", calculate_phi_psi(structure))
    except:(KeyError("KeyError: 'CA'"), FileNotFoundError)