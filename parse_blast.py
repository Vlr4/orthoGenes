import os
import gzip
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from Bio.Blast import NCBIXML

def count_upper_mis(row):
    return sum([1 for i, j in zip(row['upper_ref'], row['upper_part']) if i != j])
def count_lower_mis(row):
    return sum([1 for i, j in zip(row['lower_ref'], row['lower_part']) if i != j])
    
def get_best_value(row):
    if row['upper_mismatches'] >= row['lower_mismatches']:
        return row['lower_part']
    else:
        return row['upper_part']

def parse_results(blast_output):
    result_handle = open(TMP_DIR + blast_output)
    blast_record = NCBIXML.read(result_handle)
    filename = blast_record.alignments[0].title
    alignment_dict = {}
    count = 1
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.score > 50 and hsp.query_start != 316 and hsp.query_start != 309:
                feature_dict = {}
                feature_dict['title'] = alignment.title
                feature_dict['start'] = hsp.query_start
                feature_dict['end'] = hsp.query_end
                feature_dict['score'] = hsp.score
                feature_dict['query'] = hsp.query
                feature_dict['subject'] = hsp.sbjct
                alignment_dict[count] = feature_dict
                count += 1
    with open (RES_DIR + blast_output + '.txt', 'a') as res:
        df = pd.DataFrame(alignment_dict)
        df = df.T
        if len(df.index) > 1:
            df.sort_values(by=['start', 'end'], ascending=[True, False], inplace=True)
            df.drop_duplicates(inplace=True, keep='first') 
            df.reset_index(inplace=True, drop=True)
            rows_to_remove = []
            for i, row in df.iterrows():
                start = row['start']
                end = row['end']
                for j in range(i + 1, len(df)):
                    next_start = df.iloc[j]['start']
                    next_end = df.iloc[j]['end']
                    if start <= next_start and end >= next_end:
                        rows_to_remove.append(j)
            df = df.drop(rows_to_remove)
            df.reset_index(inplace=True, drop=True)
            df['upper_diff'] = 0
            df['lower_diff'] = 0
            for i in range(0, len(df)-1):
                    df['lower_diff'][i+1] = (df['end'][i] - df['start'][i+1])+1
            for i in range(0, len(df)-1):
                df['upper_diff'][i] = (df['end'][i] - df['start'][i+1])+1
            df['upper_diff'].clip(lower=0, inplace=True)
            df['lower_diff'].clip(lower=0, inplace=True)
            df['upper_part'] = df.apply(lambda x: x['query'][-x['upper_diff']:] if x['upper_diff'] != 0 else '', axis=1)
            df['upper_ref'] = df.apply(lambda x: x['subject'][-x['upper_diff']:] if x['upper_diff'] != 0 else '', axis=1)
            df['lower_part'] = df.apply(lambda x: x['query'][:x['lower_diff']] if x['lower_diff'] != 0 else '', axis=1)
            df['lower_ref'] = df.apply(lambda x: x['subject'][:x['lower_diff']] if x['lower_diff'] != 0 else '', axis=1)
            
            df[['lower_part', 'lower_ref']] = df[['lower_part', 'lower_ref']].shift (-1)
            df.fillna('', inplace=True)
            df['upper_mismatches'] = df.apply(count_upper_mis, axis=1)
            df['lower_mismatches'] = df.apply(count_lower_mis, axis=1)
            
            df['best_part'] = df.apply(get_best_value, axis=1)

            df['Normalized'] = df.apply(lambda x: x['subject'][x['lower_diff']:-x['upper_diff']] \
                                        if x['upper_diff'] != 0 else x['subject'][x['lower_diff']:], axis=1)
            df['Final'] = df['Normalized'] + df['best_part']  
            out = ''.join(df['Final'])
            res.write('>' + filename + '\n' + out + '\n')
        elif len(df.index) == 0:
            pass
        else:
            out = ''.join(df['subject'])
            res.write('>' + filename + '\n' + out + '\n')
        df.to_csv(RES_DIR + blast_output + '.csv')
