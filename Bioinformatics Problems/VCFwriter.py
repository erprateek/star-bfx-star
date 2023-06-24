#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import pysam
from Bio import SeqIO


class VCFwriter(object):
    def __init__(self):
        pass

    def create(fasta, sample, input_file, output_vcf):
        df = pd.read_csv(input_file, delimiter='\t')
        vcfh = pysam.VariantHeader()     
        vcfh.add_sample(sample)
        fasta_sequences = SeqIO.parse(open(fasta),'fasta')
        contigOrder = []
        
        for fasta in fasta_sequences:
            contig_name, contig_len = fasta.id, len(fasta.seq)
            contigOrder.append(contig_name)
            vcfh.add_meta('contig', items=[('ID', contig_name),
                                       ('length', contig_len)])
        fasta_sequences.close()
        contigOrderIndex = dict(zip(contigOrder,range(len(contigOrder))))
        df['#CHROM_rank'] = df['#CHROM'].map(contigOrderIndex)
        df.sort_values(['#CHROM_rank', 'POS'], \
                           ascending = [True, True], inplace = True)
        df.drop('#CHROM_rank', 1, inplace = True)
        vcfh.add_meta('FORMAT',items=[('ID',"GT"),('Number',1),('Type','String'),('Description','Genotype')])
        vcf = pysam.VariantFile(output_vcf, "w", header=vcfh)
        for i, r in df.iterrows():
            x=vcf.new_record(contig=str(r['#CHROM']), start=r['POS']-1, stop=r['POS'],
                             alleles=(r['REF'],r['ALT']),filter='PASS')
            x.samples[sample]['GT'] = (0,1)
            vcf.write(x)
        vcf.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create a VCF file')
    parser.add_argument('input_file', metavar='<input_file>', type=str,
                                        help='Path to text file with basic variant information')
    parser.add_argument('--sample', dest='sample', action='store',
                        default="<sample_unset>", required=False,
                        help='Sample ID (necessary when choosing by "FORMAT" value)')
    parser.add_argument('--fasta', dest='fasta', action='store',
                        default="<sample_unset>", required=True,
                        help='Reference genome to set VCF header against')

    options = parser.parse_args()
    v = VCFwriter()
    v.create(options.fasta, options.sample, options.input_file, '-')
