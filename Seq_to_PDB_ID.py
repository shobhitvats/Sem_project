from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import filter_sequence

def get_pdb_id_from_sequence(protein_sequence):
    """Get PDB IDs from a protein sequence using NCBI BLAST."""
    # Perform BLAST search against the PDB database
    result_handle = NCBIWWW.qblast("blastp", "pdb", protein_sequence)
    
    # Parse the BLAST results
    blast_record = NCBIXML.read(result_handle)
    
    pdb_ids = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # Extract PDB ID from the title, assuming the format is always consistent
            # Example title: "pdb|1TUP|A Chain A, ..."
            pdb_id = alignment.title.split("|")[1]
            pdb_ids.append(pdb_id)
            break  # Assuming you only want the first HSP per alignment for simplicity
    with open("filtered_pdbs.txt", "w") as output:
        output.write(str(pdb_ids)) 
    print (pdb_ids)
    return list(set(pdb_ids))  # Return unique PDB IDs


protein_sequence = [filter_sequence.dvd4,filter_sequence.dvd5,filter_sequence.dvd6]
# get_pdb_id_from_sequence(protein_sequence)
# print(pdb_ids)

def get_pdb_ids_from_sequences():
    for seq in protein_sequence:
        pdb_ids = get_pdb_id_from_sequence(seq)

    print(pdb_ids)
    return pdb_ids

for seq in protein_sequence:
    get_pdb_id_from_sequence(seq)
