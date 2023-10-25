from Bio import Entrez
from Bio import SeqIO
from urllib.request import HTTPError
from urllib.request import urlopen
import time
Entrez.email = "valeriia.dccclxiv@gmail.com"

def get_taxonomy(species):
    taxonomy_list = []
    while len(taxonomy_list) == 0: 
        try:
            handle = Entrez.esearch(db="taxonomy", term=f"{species}[Orgn]")
            record = Entrez.read(handle)
            identifier = record["IdList"][0]
            time.sleep(2)
            handle = Entrez.efetch(db="taxonomy", id=identifier, retmode='xml')
            info = Entrez.read(handle)[0] 
            for i in info['LineageEx']:
                taxonomy_list.append(i['ScientificName'])
            taxonomy = ":".join(taxonomy_list)
            time.sleep(2)
        except:
            print('Connection error, trying again')
    return taxonomy
    
taxonomy_present = []
taxonomy_missed =[]
for p in present:
    print(p)
    taxonomy_present.append(get_taxonomy(p))

for m in missed:
    print(m)
    taxonomy_missed.append(get_taxonomy(m))
    
present_dict = {}
for p, tax_p in zip(present, taxonomy_present):
    present_dict[tax_p] = p
    
    
def taxonomic_distance(taxonomy1, taxonomy2, delimiter=':'):
    levels1 = taxonomy1.split(delimiter)
    levels2 = taxonomy2.split(delimiter)
    
    for i, (level1, level2) in enumerate(zip(levels1, levels2)):
        if level1 != level2:
            return i
    return min(len(levels1), len(levels2))

def closest_taxonomy(new_taxonomy, existing_taxonomies):
    distances = [taxonomic_distance(new_taxonomy, t) for t in existing_taxonomies]
    return existing_taxonomies[distances.index(max(distances))]
