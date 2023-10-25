from Bio import Entrez
from Bio import SeqIO
from urllib.request import HTTPError
from urllib.request import urlopen
import time
from utils import get_shortname
Entrez.email = "valeriia.dccclxiv@gmail.com"

import pandas as pd

targets_path = 'test_data/'
targets_filename = 'test_set.csv'

urlpart = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'\
             '&id={}&rettype=fasta&format=text&from={}&to={}'
             
             targets_df = pd.read_csv(targets_path + targets_filename)
target_species = targets_df['species']
target_genes = targets_df['genes']
target_splist = target_species.dropna().to_list()
target_glist = target_genes.dropna().to_list()

   
def get_gene_sequence(identifier, organism):
    handle = Entrez.efetch(db="gene", id=identifier, retmode='xml')
    info = Entrez.read(handle)[0]
    locus = info['Entrezgene_locus'][0]
    acc_id = locus.get('Gene-commentary_accession')
    seq_int = locus['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
    coord = []
    coord.append(seq_int['Seq-interval_from'])
    coord.append(seq_int['Seq-interval_to'])
    with urlopen(urlpart.format(acc_id, coord[0], coord[1])) as webpage:
        sequence = webpage.read().decode("utf-8")
        lines = sequence.split('\n')
        sequence = '>' + get_shortname(organism) + '\n' + '\n'.join(lines[1:])
    return sequence
    
record_dict = {}


for organism in target_splist:
    inner_dict = {}
    #record_dict[organism] = {'shortname': record_dict[organism]}
    for gene in target_glist:
        handle = Entrez.esearch(db="gene", term=f"{organism}[Orgn] AND {gene}[Gene]")
        record = Entrez.read(handle)
        if 'ErrorList' in record:
            inner_dict = 0
        elif record["IdList"]:
            inner_dict[gene] = record["IdList"][0]    
        else:
            inner_dict[gene] = 0
    record_dict[organism] = inner_dict
    
for gene in target_glist:
    for species, info in record_dict.items():
        if isinstance(info, dict) and gene in info and info[gene] != 0:
#            print(f"{species} - {gene}: {info[gene]}")
            print(get_gene_sequence(info[gene], species))
