import pandas as pd
import argparse
import bisect
import gzip
import logging
import os
import subprocess
import time
import ntpath
from utils import setup_logger
from collections import Counter, defaultdict
from functools import partial

import numpy as np
import pandas as pd

from VCF_utils import *

# VCF I/O Functions

class VCF:
    log_file_name = __qualname__
    def __init__(self, path, verbose=False):
        self.path = path
        self.vcf_file = ntpath.basename(path)
        self.logger = setup_logger(self.log_file_name, verbose)
        self.vcf_df = self.load_vcf(self.path, verbose=verbose)
        self.vcf_df_verbose = self.vcf_df.copy(deep=True)
        self.vcf_df_verbose = separate_alleles(self.vcf_df)
        self.vcf_df_verbose = self.explode_mnvs(self.vcf_df_verbose)
        self.vcf_df_pass    = self.vcf_df[self.vcf_df['FILTER'] == 'PASS']
        self.vcf_df_pass_verbose = self.vcf_df_verbose[self.vcf_df_verbose['FILTER'] == 'PASS']

    def get_header(self, path, header_indicator="##"):
        """[Returns list of header lines from file]

        :param path: [File path]
        :type path: [str]
        :param header_indicator: [Lines starting with this string will denote the header. Must be continuous]
        :type header_indicator: [str]
        :return: [List of lines in the header]
        :rtype: [list]
        """
        header = []
        opener = gzip.open if path.endswith(".gz") else open
        open_method = "rb" if path.endswith(".gz") else "r"
        with opener(path, open_method) as infile:
            for line in infile:
                line = line.decode() if type(line) == bytes else line
                if line.startswith(header_indicator):
                    header.append(line)
                else:
                    break
        return header

    def get_sample_names(self, vcf_path):
        """[Gets a list of sample names for the given vcf]

        :param vcf_path: [vcf path]
        :type vcf_path: [str]
        :return: [list of sample names]
        :rtype: [list]
        """
        columns = self.get_header(vcf_path, "#")[-1]
        #vcf_specification = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        vcf_specification = ["#CHROM","CHROM","POS","ID","REF","ALT","QUAL","INFO","FILTER","FORMAT","ALT_string","allele_index","index"]
        sample_names = [i.strip() for i in columns.split("\t") if i not in vcf_specification]
        return sample_names

    def load_vcf(self, vcf_path, verbose=True):
        """[Loads vcf into pandas dataframe]

        :param vcf_path: [Path to VCF]
        :type vcf_path: [str]
        :param verbose: [Print things], defaults to True
        :param verbose: bool, optional
        :return: [Vcf dataframe]
        :rtype: [pandas dataframe]
        """
        start = time.time()
        df = pd.read_csv(vcf_path, header=len(self.get_header(vcf_path, "##")), sep="\t", dtype={"#CHROM": str, "POS": int})
        df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
        df.dropna(subset=["REF", "ALT"], inplace=True)
        if "ALT" not in df.columns:
            df['ALT'] = df['ALT_string']
        if verbose:
            print("VCF: {0} loaded in {1} seconds".format(vcf_path, time.time() - start))
        return df

    def explode_mnvs(self,vcf_df):
        output_df = vcf_df.copy(deep=True)
        new_rows = []
        for idx, row in vcf_df.iterrows():
            if ((len(row.REF) == len(row.ALT)) and len(row.REF) > 1):
                # Is a MNV
                self.logger.debug("Found an MNV. Here is what it is: {mnv}".format(mnv=row))
                for i in range(len(row.ALT)):
                    cur_new_row = row.copy(deep=True)
                    new_pos = cur_new_row.POS+i
                    new_ref = cur_new_row.REF[i]
                    new_alt = cur_new_row.ALT[i]
                    cur_new_row.POS = new_pos
                    cur_new_row.REF = new_ref
                    cur_new_row.ALT = new_alt
                    new_rows.append(cur_new_row)
                self.logger.debug("Exploded the MNV into {n} variants".format(n=len(new_rows)))
        new_rows_df = pd.DataFrame(new_rows)
        return output_df.append(new_rows_df)

    def write_vcf(self, vcf_df, output_path, header=None, sample_columns=None, verbose=True):
        """[Writes a VCF from the vcf dataframe, to 4.2 specification (["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]) + sample_columns]

        :param vcf_df: [VCF dataframe]
        :type vcf_df: [pandas dataframe]
        :param output_path: [Path to write VCF to]
        :type output_path: [str]
        :param header: [List of header lines to write before writing VCF], defaults to None
        :param header: [list], optional
        :param sample_columns: [List of sample names in the vcf to write sample columns for, None = no sample columns will be written], defaults to None
        :param sample_columns: [list], optional
        :param verbose: [Print things], defaults to True
        :param verbose: bool, optional
        """
        start = time.time()
        vcf_df.rename(columns={"CHROM": "#CHROM"}, inplace=True)
        vcf_df.rename(columns={"ALT_string": "ALT"}, inplace=True)
        print(vcf_df.columns)
        vcf_specification = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        vcf_specification = vcf_specification + sample_columns if sample_columns else vcf_specification
        vcf_df = vcf_df[vcf_specification]
        if header:
            with open(output_path, "w") as fout:
                for line in header:
                    fout.write(line.strip() + "\n")
        with open(output_path, "a") as fout:
            vcf_df.to_csv(fout, sep="\t", header=True, index=False)
        if verbose:
            print("VCF: {0} written in {1} seconds".format(output_path, time.time() - start))

    def get_variant_info(self, pos):
        return self.vcf_df[self.vcf_df['POS'] == pos]

    def __eq__(self, other):
        ## TODO: Implement this
        raise NotImplementedError

    def get_concordance_with(self, other):
        return vcf_concordance_sets(self.vcf_df, other.vcf_df)

    def intersect_bed(self, bed, output_path):
        bedtools_intersect_vcf(self.vcf_df, bed, output_path = output_path)
        self.bed_intersect = output_path
        return bed_intersect

    def normalize(self, output_path=""):
        output_path = ntpath.dirname(self.path)
        input_fname = ntpath.basename(self.path)
        input_split = input_fname.split(".")
        input_fname_string = input_split[0]
        input_extension_string = input_split[1:]
        output_filename = ".".join("_".join([input_fname_string,"normalized"]))
        bcftools_normalize_vcf(self.vcf_df, output_path)

    def subtract(self, other, keys=['CHROM','POS','REF','ALT']):
        """[This method fetches the set of variants that are not in the second VCF callset]
        
        Arguments:
            other {[VCF]} -- [VCF object]
        """
        df1 = self.vcf_df.copy(deep=True)

        # Subtract out the rows
        df1_subset = df1[keys]
        df2_subset = other.vcf_df[keys]
        df1_subset['variant'] = df1_subset[keys[0]].astype(str)+"|"+df1_subset[keys[1]].astype(str)+"|"+df1_subset[keys[2]].astype(str)+"|"+df1_subset[keys[3]].astype(str)
        df1_subset.set_index(['variant'],inplace=True)
        df2_subset['variant'] = df2_subset[keys[0]].astype(str)+"|"+df2_subset[keys[1]].astype(str)+"|"+df2_subset[keys[2]].astype(str)+"|"+df2_subset[keys[3]].astype(str)
        df2_subset.set_index(['variant'],inplace=True)
        diff = df1_subset[~df1_subset.index.isin(df2_subset.index)].dropna()
        self.subtracted_df = diff
        return diff
    
    def get_snvs_df(self, df):
        snvs = df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
        return snvs
    
    def get_indels_df(self, df):
        indels = df[((df['REF'].str.len() > 1) | (df['ALT'].str.len() > 1)) & (df['REF'].str.len() != df['ALT'].str.len())]

    def print_af_stats(self, df, af):
        """[This method outputs various metrics for the variant attributes]
        
        Arguments:
            df {[pd.DataFrame]} -- [Input dataframe derived from the VCF. If there is no AF column, it will be generated and used.
            Also, if there are > 1 samples in the VCF, only the first one will be used.]
            af {[float]} -- [Allele frequency to filter the VCF on]
        """
        if "AF" not in list(df.columns):
            sample_names = self.get_sample_names(self.path)
            ad_values, dp_values = parse_sample(df, sample_name=sample_names[0], keys=['AD','DP'], parse_alleles=True, coerce=True)
            ad_values_array = np.array(ad_values)
            dp_values_array = np.array(dp_values)
            df['AF'] = ad_values_array/dp_values_array

        snps = df.loc[((df['REF'].str.len() == df['ALT'].str.len()) & (df['REF'].str.len()==1))]
        indels = df.loc[df['REF'].str.len() != df['ALT'].str.len()]
        mnvs = df.loc[((df['REF'].str.len() == df['ALT'].str.len()) & (df['REF'].str.len()==2))]
        nsnps_gt_af = snps[snps['AF'] >= af].shape[0]
        nindels_gt_af = indels[indels['AF'] >= af].shape[0]
        nmnvs_gt_af = mnvs[mnvs['AF'] >= af].shape[0]
        self.logger.info("######################### AF  = {af} #########################".format(af=af))
        self.logger.info("SNPs > {af} = {nsnps}".format(af=af, nsnps=nsnps_gt_af))
        self.logger.info("Indels > {af} = {nindels}".format(af=af, nindels=nindels_gt_af))
        self.logger.info("MNVs > {af} = {nmnvs}".format(af=af, nmnvs=nmnvs_gt_af))
        self.logger.info("#"*60)

    def print_stats(self):
        """[Output to screen the relevant stats]
        """
        print("VCF contains the following: ")
        af_list = [
            0.1,
            0.05,
            0.2
        ]
        # This makes a call using the af list and can slice and dice metrics for a VCF
        vcf_df_with_af = add_af(self.vcf_df_verbose, self.get_sample_names(self.path)[0], coerce=True)
        for af in af_list:
            self.print_af_stats(df = vcf_df_with_af, af=af)
        #map(partial(self.print_af_stats, df=vcf_df_with_af), af_list)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--vcf_path',
                        help= "Input path to the VCF under consideration")
    args = parser.parse_args()
    vcf_obj = VCF(args.vcf_path)
    vcf_obj.print_stats()

if __name__ == '__main__':
    main()
