import re

from Bio import SeqIO

def extract_sequences(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

fasta_file = "1.fasta"
protein_sequences = extract_sequences(fasta_file)

print(protein_sequences)


def filter_protein_sequences(sequences):
    # """Filter proteins that have two cysteines with any two amino acids in between."""
    pattern = re.compile(r'C..C')  # Regular expression for the pattern
    filtered_sequences = [seq for seq in sequences if pattern.search(seq)]
    return filtered_sequences

filtered_sequences = filter_protein_sequences(protein_sequences)
# print(f"Filtered Sequences: {filtered_sequences}")
# print(filtered_sequences)

better_list = "\n\n".join(filtered_sequences)

# print (better_list)

with open("all_sequences.txt", "w") as output:
    output.write(str(better_list ))
    
dvd1 = "\n\n".join(filtered_sequences[:357])
dvd2 = "\n\n".join(filtered_sequences[357:714])
dvd3 = "\n\n".join(filtered_sequences[714:1071])
dvd4 = "\n\n".join(filtered_sequences[1071:1228])
dvd5 = "\n\n".join(filtered_sequences[1228:1385])
dvd6 = "\n\n".join(filtered_sequences[1385:1542])    
dvd7 = "\n\n".join(filtered_sequences[1542:2142])

with open("dvd1.txt", "w") as output:
    output.write(str(dvd1 ))
    
with open("dvd2.txt", "w") as output:
    output.write(str(dvd2))
    
with open("dvd3.txt", "w") as output:
    output.write(str(dvd3))    
    
with open("dvd4.txt", "w") as output:
    output.write(str(dvd4 ))    
    
with open("dvd5.txt", "w") as output:
    output.write(str(dvd5 ))
    
with open("dvd6.txt", "w") as output:
    output.write(str(dvd6 ))
    
with open("dvd7.txt", "w") as output:
    output.write(str(dvd7 ))   


