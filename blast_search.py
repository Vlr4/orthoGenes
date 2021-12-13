import os
import glob
from blast_parser import parse
    
def extract_genomes():
    cmd = "gunzip ./subjects/*.fna.gz"
    os.system(cmd)


def blast_search(file):
    cmd = "blastn -word_size 25 -query test_query.fna -subject " + file + " -out out.txt -outfmt 6"
    os.system(cmd)

def convert(x):
    alt_map = {'N': '0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
    for k,v in alt_map.items():
        x = x.replace(k,v)
    bases = list(x) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def linearize_genome(genome):
    g = open(genome, "r")
    linear_genome = open("linear_genome.txt", "w+")
    first_line = g.readline()
    linear_genome.write(first_line.rstrip("\n") + "sequenceStart")
    for line in g:
        header = ""
        if line.startswith(">"):
            stripped_line = line.rstrip("\n")
            header += "\n" + stripped_line + "sequenceStart"
            linear_genome.write(header)
        else:
            stripped_line = line.rstrip("\n")
            linear_genome.write(stripped_line)
     
    linear_genome.close()
    g.close()

def cut_genes(qid, sid, start, end, compl, name):
    start = start-1
    genes = open("genes.txt", "a+")
    with open ("linear_genome.txt", "r") as l_genome:
       lines = l_genome.readlines()
       for line in lines:
           if line.startswith(">" + sid):
               splitted = line.split("sequenceStart")
               seq = splitted[1]
               gene = seq[start:end]
               if compl == True:
                   gene = convert(gene)
               else:
                   gene = gene.upper()
               genes.write(">" + qid + "\t" + name + "  " + sid + "\n" + gene + "\n")
    genes.close()
    l_genome.close()
    
def define_coordinates(out, name):
    for blast_record in parse(out):
        hit_num = 1
        
        print('query id: {}'.format(blast_record.qid))
        for hit in blast_record.hits:
            start_points = []
            end_points = []
            for hsp in hit:
                start_points.append(hsp.sstart)
                end_points.append(hsp.send)
            print(start_points)
            print(end_points)
            if max(start_points) > max(end_points):
                rev_com = True
                start = min(end_points) 
                end = max(start_points) 
            else:
                start = min(start_points) 
                rev_com = False
                end = max(end_points) 
            sid = hsp.sid
            qid = hsp.qid
            print(sid)
            print(start, end)

            cut_genes(qid, sid, start, end, rev_com, name)
    out.close()
    
def main():
    extract_genomes()
    for file in glob.glob("./subjects/*"):
        genome = file
        name = str(genome[11:])
        blast_search(genome)
        out = open("out.txt")
        linearize_genome(genome)
        define_coordinates(out, name)
        rename_res = "mv genes.txt ortol_" + name + ".txt"
        os.system(rename_res)
        
main()




