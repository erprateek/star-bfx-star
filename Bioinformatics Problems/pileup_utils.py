from io import StringIO
from collections import Counter
import pandas as pd
import re

def parse_pileup_row(row):
    if row:
        pileup_row = StringIO(row)
        pileup_df = pd.read_csv(pileup_row, header=None,sep="\t")
        pileup_df.columns = ['CHROM','POS','REF','NREADS','pileup_string','qual_string']
        return pileup_df, Counter(pileup_df['pileup_string'].tolist()[0])
    else:
        return pd.DataFrame(), Counter(row)

def parse_pileup_row2(row, col_prefix=''):
    if row:
        pileup_row = StringIO(row)
        pileup_df = pd.read_csv(pileup_row, header=None,sep="\t")
        pileup_df.columns = ['CHROM','POS','REF',"_".join([col_prefix, 'NREADS']),"_".join([col_prefix, 'pileup_string']),"_".join([col_prefix,'qual_string'])]
        return pileup_df, Counter(pileup_df["_".join([col_prefix,'pileup_string'])].tolist()[0])
    else:
        return pd.DataFrame(), Counter(row)


def get_indel_from_pileup_string(s):
    pattern = re.compile(r'[+-]\d+[ATGCatgc]*')
    indel_patterns = re.findall(pattern, s)
    unique_indel_patterns = [item.upper() for item in indel_patterns]
    return (Counter(indel_patterns), Counter(unique_indel_patterns))
