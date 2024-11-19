import os
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, simpledialog, Toplevel, Checkbutton, IntVar, Button, Label
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.vectors import calc_dihedral
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import requests
import re

class BioinformaticsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("PDB Structural Analysis App")
        self.root.geometry("600x400")
        self.root.configure(bg='#1e1e1e')
        self.input_dir = ""
        self.output_dir = ""
        self.selected_amino_acids = []
        self.create_widgets()

    def create_widgets(self):
        frame = tk.Frame(self.root, bg='#1e1e1e', highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=2)
        frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        self.process_button = tk.Button(frame, text="Process Ramachandran Plot", command=self.process_ramachandran, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 12, 'bold'), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.process_button.pack(pady=10)

        self.download_pdb_button = tk.Button(frame, text="Download PDB Files", command=self.download_pdb_files, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 12, 'bold'), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.download_pdb_button.pack(pady=10)

        self.ramachandran_single_button = tk.Button(frame, text="Ramachandran Plot (Single)", command=self.generate_ramachandran_single_gui, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 12, 'bold'), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.ramachandran_single_button.pack(pady=10)

        self.ramachandran_multiple_button = tk.Button(frame, text="Ramachandran Plot (Multiple)", command=self.generate_ramachandran_multiple_gui, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 12, 'bold'), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.ramachandran_multiple_button.pack(pady=10)

        self.filter_motif_button = tk.Button(frame, text="Filter Motif", command=self.filter_motif, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 12, 'bold'), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.filter_motif_button.pack(pady=10)

        self.status_text = scrolledtext.ScrolledText(frame, height=10, bg='#1e1e1e', fg='#d4af37', font=('Helvetica', 10), highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1)
        self.status_text.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

    def update_status(self, message):
        self.status_text.insert(tk.END, message + "\n")
        self.status_text.see(tk.END)
        self.root.update()

    def select_input_directory(self):
        self.input_dir = filedialog.askdirectory(title="Select Input Directory")
        if self.input_dir:
            self.update_status(f"Selected input directory: {self.input_dir}")

    def select_output_directory(self):
        self.output_dir = filedialog.askdirectory(title="Select Output Directory")
        if self.output_dir:
            self.update_status(f"Selected output directory: {self.output_dir}")

    def select_amino_acids(self):
        amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        self.selected_amino_acids = []

        def on_select():
            self.selected_amino_acids.clear()
            for i, var in enumerate(vars):
                if var.get():
                    self.selected_amino_acids.append(amino_acids[i])
            dialog.destroy()
            self.plot_ramachandran_all()

        dialog = Toplevel(self.root)
        dialog.title("Select Amino Acids")
        dialog.geometry("300x400")
        dialog.configure(bg='#1e1e1e')

        vars = [IntVar() for _ in amino_acids]
        for i, aa in enumerate(amino_acids):
            row = i % 10
            col = i // 10
            Checkbutton(dialog, text=aa, variable=vars[i], bg='#1e1e1e', fg='#d4af37', selectcolor='#1e1e1e').grid(row=row, column=col, sticky='w')

        Button(dialog, text="Proceed", command=on_select, bg='#1e1e1e', fg='#d4af37', highlightbackground='#d4af37', highlightcolor='#d4af37', highlightthickness=1).grid(row=10, columnspan=2, pady=10)

        dialog.transient(self.root)
        dialog.grab_set()
        self.root.wait_window(dialog)

        if not self.selected_amino_acids:
            messagebox.showinfo("Info", "No amino acids selected.")
            return

        self.update_status(f"Selected amino acids: {self.selected_amino_acids}")

    def download_pdb_files(self):
        input_file = filedialog.askopenfilename(title="Select Input File", filetypes=(("Text Files", "*.txt"), ("All Files", "*.*")))
        if not input_file:
            messagebox.showinfo("Info", "No input file selected.")
            return

        output_dir = filedialog.askdirectory(title="Select Output Directory")
        if not output_dir:
            messagebox.showinfo("Info", "No output directory selected.")
            return

        with open(input_file, 'r') as file:
            pdb_ids = [line.strip() for line in file if line.strip()]

        total_files = len(pdb_ids)
        downloaded_files = 0

        for pdb_id in pdb_ids:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(url)
            if response.status_code == 200:
                with open(os.path.join(output_dir, f"{pdb_id}.pdb"), 'w') as pdb_file:
                    pdb_file.write(response.text)
                downloaded_files += 1
                self.update_status(f"Downloaded {pdb_id}.pdb ({downloaded_files}/{total_files})")
            else:
                self.update_status(f"Failed to download {pdb_id}.pdb")

        self.update_status(f"Download complete. {downloaded_files}/{total_files} files downloaded.")


    
    def extract_motif_coordinates(self, directory, selected_amino_acids):
        parser = PDBParser()
        coordinates = {aa: {1: [], 2: [], 3: [], 4: []} for aa in selected_amino_acids}
        files = [f for f in os.listdir(directory) if f.endswith(".pdb")]
        total_files = len(files)
        processed_files = 0

        for filename in files:
            pdb_file = os.path.join(directory, filename)
            self.update_status(f"Processing file: {pdb_file}")
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
                                    distance = residues[i]['SG'] - residues[i+2]['N']
                                    if distance <= 3.5:
                                        coordinates[aa][1].append((residues[i-1]['C'].get_vector(), residues[i]['N'].get_vector(), residues[i]['CA'].get_vector(), residues[i]['C'].get_vector(), residues[i+1]['N'].get_vector()))
                                        coordinates[aa][2].append((residues[i]['C'].get_vector(), residues[i+1]['N'].get_vector(), residues[i+1]['CA'].get_vector(), residues[i+1]['C'].get_vector(), residues[i+2]['N'].get_vector()))
                                        coordinates[aa][3].append((residues[i+1]['C'].get_vector(), residues[i+2]['N'].get_vector(), residues[i+2]['CA'].get_vector(), residues[i+2]['C'].get_vector(), residues[i+3]['N'].get_vector()))
                                        coordinates[aa][4].append((residues[i+2]['C'].get_vector(), residues[i+3]['N'].get_vector(), residues[i+3]['CA'].get_vector(), residues[i+3]['C'].get_vector(), residues[i+4]['N'].get_vector()))
                                        self.update_status(f"Motif {aa} found in file: {pdb_file} at position {i} with distance {distance:.2f} Å")
            processed_files += 1
            self.update_status(f"Processed {processed_files}/{total_files} files")

        return coordinates

    def calculate_phi_psi(self, coordinates):
        phi_psi_angles = {aa: {1: [], 2: [], 3: [], 4: []} for aa in coordinates}
        for aa in coordinates:
            for key in coordinates[aa]:
                for vectors in coordinates[aa][key]:
                    phi = calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])
                    psi = calc_dihedral(vectors[1], vectors[2], vectors[3], vectors[4])
                    phi_psi_angles[aa][key].append((np.rad2deg(phi), np.rad2deg(psi)))
                    self.update_status(f"Calculated phi: {np.rad2deg(phi)}, psi: {np.rad2deg(psi)} for {aa} at position {key}")
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
            output_path = os.path.join(self.output_dir, f'{amino_acid}_ramachandran_plot_position_{position}.png')
            plt.savefig(output_path)
            plt.close()
            self.update_status(f"Plot saved for {amino_acid} at position {position} as {output_path}")
        else:
            self.update_status(f"No data to plot for {amino_acid} at position {position}")

    def plot_ramachandran_all(self):
        if not self.input_dir or not self.output_dir or not self.selected_amino_acids:
            messagebox.showinfo("Info", "Please select input/output directories and amino acids.")
            return

        coordinates = self.extract_motif_coordinates(self.input_dir, self.selected_amino_acids)
        self.update_status("Coordinates extraction complete.")
        phi_psi_angles = self.calculate_phi_psi(coordinates)
        self.update_status("Phi and Psi angle calculation complete.")

        for aa in phi_psi_angles:
            for position in range(1, 5):
                self.plot_ramachandran(phi_psi_angles[aa][position], aa, position)
                self.update_status(f"Plotting complete for {aa} at position {position}")

    def process_ramachandran(self):
        self.select_input_directory()
        if not self.input_dir:
            return
        self.select_output_directory()
        if not self.output_dir:
            return
        self.select_amino_acids()

    def filter_motif(self):
        amino_acid = simpledialog.askstring("Input", "Enter Amino Acid (e.g., c, d, e):", parent=self.root)
        if amino_acid:
            amino_acid = amino_acid.upper()
            file_path = filedialog.askopenfilename(title="Select Sequence File", filetypes=(("Text Files", "*.txt"), ("All Files", "*.*")))
            if file_path:
                motif_pattern = f"{amino_acid}..{amino_acid}"
                seq_list = []
                matches = []
                pdb_chain_ids = []
                current_pdb_chain_id = ""
                total_motif_count = 0
                sequences_with_motif = 0
                sequences_with_motif_more_than_once = 0
                sequences_with_2_motif = 0
                sequences_with_3_motif = 0
                sequences_with_4_motif = 0
                sequences_with_5_motif = 0

                with open(file_path, 'r') as file:
                    for line in file:
                        if line.startswith('>'):
                            current_pdb_chain_id = line[1:6]
                        else:
                            current_motifs = re.findall(motif_pattern, line)
                            if current_motifs:
                                total_motif_count += len(current_motifs)
                                pdb_chain_ids.append(current_pdb_chain_id)
                                if len(current_motifs) > 1:
                                    sequences_with_motif_more_than_once += 1
                                if len(current_motifs) == 2:
                                    sequences_with_2_motif += 1
                                if len(current_motifs) == 3:
                                    sequences_with_3_motif += 1
                                if len(current_motifs) == 4:
                                    sequences_with_4_motif += 1
                                if len(current_motifs) == 5:
                                    sequences_with_5_motif += 1

                self.update_status(f"Total motifs found: {total_motif_count}")
                self.update_status(f"PDB Chain IDs with motif: {pdb_chain_ids}")
                self.update_status(f"Sequences with motif more than once: {sequences_with_motif_more_than_once}")
                self.update_status(f"Sequences with 2 motif: {sequences_with_2_motif}")
                self.update_status(f"Sequences with 3 motif: {sequences_with_3_motif}")
                self.update_status(f"Sequences with 4 motif: {sequences_with_4_motif}")
                self.update_status(f"Sequences with 5 motif: {sequences_with_5_motif}")

                output_file_path = filedialog.asksaveasfilename(defaultextension=".py", filetypes=[("Python Files", "*.py")])
                if output_file_path:
                    with open(output_file_path, 'w') as output_file:
                        output_file.write(f"seq_list = {seq_list}\n")
                        output_file.write(f"{amino_acid} = {matches}\n")
                        output_file.write(f"pdb_chain_ids = {pdb_chain_ids}\n")

                chain_ids_file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
                if chain_ids_file_path:
                    np.savetxt(chain_ids_file_path, np.array(pdb_chain_ids, dtype=str), fmt='%s', header='Chain IDs')

                messagebox.showinfo("Success", "Output files have been saved.")

    def generate_ramachandran_single_gui(self):
        pdb_file = filedialog.askopenfilename(title="Select PDB File", filetypes=(("PDB Files", "*.pdb"), ("All Files", "*.*")))
        if pdb_file:
            self.generate_ramachandran_single(pdb_file)

    def generate_ramachandran_single(self, pdb_file):
        parser = PDBParser()
        structure = parser.get_structure('PDB', pdb_file)
        ppb = PPBuilder()
        phi_psi = []
        for pp in ppb.build_peptides(structure):
            phi_psi.extend(pp.get_phi_psi_list())

        phi_psi = [angles for angles in phi_psi if None not in angles]
        phi, psi = zip(*phi_psi)
        phi = np.rad2deg(phi)
        psi = np.rad2deg(psi)

        plt.scatter(phi, psi, s=1)
        pdb_code = os.path.basename(pdb_file).replace('.pdb', '')
        default_filename = f"{pdb_code}_ramachandran_plot.png"

        root = tk.Tk()
        root.withdraw()

        output_file_path = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(pdb_file),
            title="Save Ramachandran Plot",
            filetypes=[("PNG files", "*.png")],
            defaultextension=".png",
            initialfile=default_filename
        )

        if output_file_path:
            plt.title(f"Ramachandran Plot for {pdb_code}")
            plt.xlabel('Phi')
            plt.ylabel('Psi')
            plt.xlim(-180, 180)
            plt.ylim(-180, 180)
            plt.grid(True)
            plt.savefig(output_file_path)
            plt.close()
            self.update_status(f"Ramachandran plot saved as {output_file_path}")
        else:
            plt.close()

        messagebox.showinfo("Success", f"Ramachandran plot saved as {output_file_path}")

    def generate_ramachandran_multiple_gui(self):
        directory = filedialog.askdirectory(title="Select Directory Containing PDB Files")
        if directory:
            self.generate_ramachandran_multiple(directory)

    def generate_ramachandran_multiple(self, directory):
        root = tk.Tk()
        root.withdraw()

        pdb_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".pdb")]
        num_files = len(pdb_files)
        parser = PDBParser()
        ppb = PPBuilder()

        plt.figure(figsize=(15, 15))

        for pdb_file in pdb_files:
            structure = parser.get_structure('PDB', pdb_file)
            for pp in ppb.build_peptides(structure):
                phi_psi = pp.get_phi_psi_list()
                phi_psi = [angles for angles in phi_psi if None not in angles]
                phi, psi = zip(*phi_psi)
                phi = np.rad2deg(phi)
                psi = np.rad2deg(psi)
                if num_files <= 5:
                    plt.scatter(phi, psi, s=1, label=os.path.basename(pdb_file))
                else:
                    plt.scatter(phi, psi, s=1)

        plt.title("Ramachandran Plot for Multiple PDBs")
        plt.xlabel('Phi')
        plt.ylabel('Psi')
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.grid(True)

        if num_files <= 5:
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.tight_layout()

        default_filename = f"ramachandran_plot_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.png"

        output_file_path = filedialog.asksaveasfilename(
            initialdir=directory,
            title="Save Ramachandran Plot",
            filetypes=[("PNG files", "*.png")],
            defaultextension=".png",
            initialfile=default_filename
        )

        if output_file_path:
           plt.savefig(output_file_path)
           plt.close()
           self.update_status(f"Ramachandran plot saved as {output_file_path}")
        else:
            plt.close()

        messagebox.showinfo("Success", f"Ramachandran plot saved as {output_file_path}")

if __name__ == "__main__":
    root = tk.Tk()
    app = BioinformaticsApp(root)
    root.mainloop()