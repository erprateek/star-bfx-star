# Class for managing BAM files and associated data
import os
import sys
import pysam
import argparse
import pandas as pd
import logging
import subprocess as sp
import shlex
from utils import setup_logger
#from resources import *

class BAM(object):
    classname = __qualname__
    def __init__(self, bam_path, verbose=False):
        """ new BAM object from the BAM/SAM at bam_path
        """
        self.logger = setup_logger(self.classname, verbose)
        self.path = bam_path

    def get_reads_at_position(self, chrom, pos_start, pos_end):
        """ Gets reads at the specified position from the BAM file.
        Returns a list of reads.
        
        Arguments:
            chrom {String} -- Chromosome
            pos_start {Int} -- Start Position
            pos_end {Int} -- End position
        """
        reads = []
        with pysam.AlignmentFile(self.path, 'rb') as samfile:
            for read in samfile.fetch("%s" % chrom, pos_start, pos_end):
                reads.append(read)
            if reads == []:
                self.logger.info("No reads found")
        return reads
    
    def generate_stats(self, bam, SAMTOOLS=''):
        ## TODO: Implement this
        #self.logger.info("Generating stats for {b}".format(b=bam))
        command = '{samtools} stats {bam}'
        return command
        
        
    def get_pileup(self, sample_name_sitelist_bam_tuple):
        sample_name = sample_name_sitelist_bam_tuple[0]
        sample_sitelist = sample_name_sitelist_bam_tuple[1]
        sample_bam = sample_name_sitelist_bam_tuple[2]
        ref_fasta = REF
        outfile = os.path.join(self.pileup_dir,sample_name+".pileup")
        
        if not os.path.exists(outfile):
            pileup_command = '{samtools} mpileup {bam} -f {ref} -l {sample_sitelist} -o {output_pileup}'.format(
                samtools=SAMTOOLS, ref=ref_fasta, sample_sitelist=sample_sitelist, output_pileup=outfile,bam=sample_bam
            )
            self.logger.info("Executing command: %s"%(pileup_command))
            with open(os.devnull, 'w') as devnull:
                #ret_code = sp.check_output(pileup_command, shell=True).split('\n')
                p = sp.Popen(pileup_command, shell=True)
                return p
        else:
            self.logger.info("Pileup exists. Skipping generation.")
        return None

    def get_pileup_for_variant(self, var):
        ref = REF
        pileup_command = f"{SAMTOOLS} mpileup {self.path} -f {ref} -r {var}"
        pileup_command = shlex.split(pileup_command)
        s=sp.check_output(pileup_command)
        return s.decode("utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","bam_file",help="Input bam file")
    args = parser.parse_args()
    BAM(args.bam_file)
    logging.info("BAM object created successfully.")

if __name__ == '__main__':
    main()

