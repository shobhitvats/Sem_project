
import requests

def download_pdb_file(pdb_chain_id):
    pdb_id = pdb_chain_id[:4] # PDB IDs are case-insensitive and usually provided in lowercase
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(f"./Large_Data/PDB/{pdb_chain_id}.pdb", "w") as file:
            file.write(response.text)
        print(f"Downloaded {pdb_chain_id}.pdb")
    else:
        print(f"Failed to download PDB file for {pdb_chain_id}. Status code: {response.status_code}")

# opening the file in read mode 
my_file = open("./Large_Data/len_list.txt", "r") 
  
# reading the file 
data = my_file.read() 
  
# replacing end of line('/n') with ' ' and 
# splitting the text it further when '.' is seen. 
data_into_list = data.split("\n") 
  
# printing the data 
# print(data_into_list) 
# my_file.close() 



for pdb_chain_id in data_into_list:
    download_pdb_file(pdb_chain_id)