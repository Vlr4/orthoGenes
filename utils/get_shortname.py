from Bio import Entrez
from Bio import SeqIO
from urllib.request import HTTPError
from urllib.request import urlopen
import time
Entrez.email = "valeriia.dccclxiv@gmail.com"
targets_path = '/home/vlr/Documents/Projects/Zoology/Programs_codes/Blast_gene_search/Orthogenes/test_data/'
targets_filename = 'test_set.csv'
import pandas as pd



def get_shortname(name):
    genus, species = name.split()
    short_genus = genus[:3]
    short_species = species[:3]
    short_name = f"{short_genus}_{short_species}"
    return short_name
