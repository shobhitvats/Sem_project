import os
import tkinter as tk
from tkinter import filedialog, messagebox
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
import matplotlib.pyplot as plt
import numpy as np

class BioinformaticsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Bioinformatics App")
        
        self.input_dir = ""
        self.output_dir = ""
        self.selected_amino_acids = []
        
        self.create_widgets()
    
    def create_widgets(self):
        frame = tk.Frame(self.root)
        frame.pack(padx=10, pady=10)
        
        self.input_button = tk.Button(frame, text="Select Input Directory", command=self.select_input_directory)
        self.input_button.grid(row=0, column=0, padx=5, pady=5)
        
        self.output_button = tk.Button(frame, text="Select Output Directory", command=self.select_output_directory)
        self.output_button.grid(row=0, column=1, padx=5, pady=5)
        
        self.amino_acid_button = tk.Button(frame, text="Select Amino Acids", command=self.select_amino_acids)
        self.amino_acid_button.grid(row=1, column=0, columnspan=2, padx=5, pady=5)
        
        self.plot_button = tk.Button(frame, text="Plot Ramachandran", command=self.plot_ramachandran)
        self.plot_button.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
        
    def select_input_directory(self):
        self.input_dir = filedialog.askdirectory(title="Select Input Directory")
        if self.input_dir:
            print(f"Selected input directory: {self.input_dir}")
    
    def select_output_directory(self):
        self.output_dir = filedialog.askdirectory(title="Select Output Directory")
        if self.output_dir:
            print(f"Selected output directory: {self.output_dir}")
    
    def select_amino_acids(self):
        amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        self.selected_amino_acids = []
        
        def on_select():
            self.selected_amino_acids.clear()
            for i, var in enumerate(vars):
                if var.get():
                    self.selected_amino_acids.append(amino_acids[i])
            root.quit()

        vars = [tk.IntVar() for _ in amino_acids]
        for i, aa in enumerate(amino_acids):
            tk.Checkbutton(root, text=aa, variable=vars[i]).pack(anchor='w')

        tk.Button(root, text="Select", command=on_select).pack()
        root.deiconify()
        root.mainloop()

        if not self.selected_amino_acids:
            messagebox.showinfo("Info", "No amino acids selected.")
            return

        print(f"Selected amino acids: {self.selected_amino_acids}")

    def extract_motif_coordinates(self, directory, selected_amino_acids):
        parser = PDBParser()
        coordinates = {aa: {1: [], 2: [], 3: [], 4: []} for aa in selected_amino_acids}

        for filename in os.listdir(directory):
            if filename.endswith(".pdb"):
                pdb_file = os.path.join(directory, filename)
                print(f"Processing file: {pdb_file}")
                structure = parser.get_structure('PDB_structure', pdb_file)
                
                for model in structure:
                    for chain in model:
                        residues = list(chain)
                        for i in range(1, len(residues) - 4):
                            for aa in selected_amino_acids:
                                if residues[i].get_resname() == 'CYS' and residues[i+3].get_resname() == aa:
                                    if all(atom in residues[i-1] for atom in ['C']) and \
                                       all(atom in residues[i] for atom in ['N', 'CA', 'C', 'SG']) and \
                                       all(atom in residues[i+1] for atom in ['N', 'CA', 'C']) and \
                                       all(atom in residues[i+2] for atom in ['N', 'CA', 'C']) and \
                                       all(atom in residues[i+3] for atom in ['N', 'CA', 'C']) and \
                                       all(atom in residues[i+4] for atom in ['N']):
                                        # Calculate the distance between the nitrogen atom of i+2 and the sulfur atom of i
                                        distance = residues[i]['SG'] - residues[i+2]['N']
                                        if distance <= 3.5:
                                            coordinates[aa][1].append((residues[i-1]['C'].get_vector(), residues[i]['N'].get_vector(), residues[i]['CA'].get_vector(), residues[i]['C'].get_vector(), residues[i+1]['N'].get_vector()))
                                            coordinates[aa][2].append((residues[i]['C'].get_vector(), residues[i+1]['N'].get_vector(), residues[i+1]['CA'].get_vector(), residues[i+1]['C'].get_vector(), residues[i+2]['N'].get_vector()))
                                            coordinates[aa][3].append((residues[i+1]['C'].get_vector(), residues[i+2]['N'].get_vector(), residues[i+2]['CA'].get_vector(), residues[i+2]['C'].get_vector(), residues[i+3]['N'].get_vector()))
                                            coordinates[aa][4].append((residues[i+2]['C'].get_vector(), residues[i+3]['N'].get_vector(), residues[i+3]['CA'].get_vector(), residues[i+3]['C'].get_vector(), residues[i+4]['N'].get_vector()))
                                            print(f"Motif {aa} found in file: {pdb_file} at position {i} with distance {distance:.2f} Å")
        return coordinates

    def calculate_phi_psi(self, coordinates):
        phi_psi_angles = {aa: {1: [], 2: [], 3: [], 4: []} for aa in coordinates}
        for aa in coordinates:
            for key in coordinates[aa]:
                for vectors in coordinates[aa][key]:
                    phi = calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])
                    psi = calc_dihedral(vectors[1], vectors[2], vectors[3], vectors[4])
                    phi_psi_angles[aa][key].append((np.rad2deg(phi), np.rad2deg(psi)))
                    print(f"Calculated phi: {np.rad2deg(phi)}, psi: {np.rad2deg(psi)} for {aa} at position {key}")
        return phi_psi_angles

    def plot_ramachandran(self, phi_psi_angles, amino_acid, position):
        if phi_psi_angles:
            phis, psis = zip(*phi_psi_angles)
            plt.scatter(phis, psis, s=3, c='blue', alpha=0.5)
            plt.title(f'Ramachandran Plot for {amino_acid} at Position {position}')
            plt.xlabel('Phi (°)')
            plt.ylabel('Psi (°)')
            plt.xlim(-180, 180)
            plt.ylim(-180, 180)
            plt.grid(True)
            plt.savefig(os.path.join(self.output_dir, f'{amino_acid}_ramachandran_plot_position_{position}.png'))
            plt.close()
            print(f"Plot saved for {amino_acid} at position {position}")
        else:
            print(f"No data to plot for {amino_acid} at position {position}")

    def plot_ramachandran(self):
        if not self.input_dir or not self.output_dir or not self.selected_amino_acids:
            messagebox.showinfo("Info", "Please select input/output directories and amino acids.")
            return

        coordinates = self.extract_motif_coordinates(self.input_dir, self.selected_amino_acids)
        print("Coordinates extraction complete.")
        phi_psi_angles = self.calculate_phi_psi(coordinates)
        print("Phi and Psi angle calculation complete.")
        
        for aa in phi_psi_angles:
            for position in range(1, 5):
                self.plot_ramachandran(phi_psi_angles[aa][position], aa, position)
                print(f"Plotting complete for {aa} at position {position}")

if __name__ == "__main__":
    root = tk.Tk()
    app = BioinformaticsApp(root)
    root.mainloop()