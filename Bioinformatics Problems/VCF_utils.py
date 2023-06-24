import bisect
import gzip
import os
import time
import numpy as np
import pandas as pd
import subprocess
from collections import Counter, defaultdict
from functools import partial

def separate_alleles(vcf_df, verbose=True):
    """[Explodes multi-allele calls into individual rows in a pandas dataframe, adds new column "allele_index" for sub-indexing of data in the format/info fields]

    :param vcf_df: [VCF dataframe]
    :type vcf_df: [pandas dataframe]
    :return: [VCF dataframe with exploded allele records]
    :rtype: [pandas dataframe]
    """
    start = time.time()

    # Format dataframe
    vcf_df.rename(columns={"ALT": "ALT_string"}, inplace=True)

    # Explode alts
    exploded_alts = pd.DataFrame(vcf_df["ALT_string"].str.split(',', expand=True).stack().reset_index(level=1, drop=True), columns=["ALT"])

    # Annotate allele sub-index
    cnt = Counter()
    counts = []
    for i in exploded_alts.index:
        counts.append(cnt[i] + 1)
        cnt[i] += 1
    exploded_alts["allele_index"] = counts

    # Join
    vcf_df = vcf_df.join(exploded_alts, how="left")
    vcf_df['ALT'] = vcf_df['ALT_string']

    #if verbose:
    #    print("VCF records separated in {0} seconds".format(time.time() - start))
    return vcf_df

def remove_ref_calls(vcf_df, remove_N=True, verbose=True):
    """[Returns a vcf dataframe with all reference calls removed (denoted by alt = <*>), optionally all 'N' calls removed as well]

    :param vcf_df: [VCF dataframe, would typically be a naive callset / gvcf]
    :type vcf_df: [pandas dataframe]
    :param remove_N: [Remove any ref/alt pair record with an N as well], defaults to True
    :param remove_N: bool, optional
    :param verbose: [Print things], defaults to True
    :param verbose: bool, optional
    :return: [VCF dataframe with ref calls and potentially N calls removed]
    :rtype: [pandas dataframe]
    """
    start = time.time()
    if remove_N:
        vcf_df = vcf_df[(vcf_df.REF != "N") & (vcf_df.ALT != "N") & (vcf_df.ALT != "<*>")]
    else:
        vcf_df = vcf_df[vcf_df.ALT != "<*>"]
    if verbose:
        print("Removed ref calls{0} in {1} seconds.".format(" and N calls" if remove_N else "", time.time() - start))
    return vcf_df

def parse_info(vcf_df, keys, parse_alleles=False, coerce=False, verbose=True):
    """[Parse information from the VCF's INFO field. An example usage would be:
    ad_values, dp_values = parse_info(vcf_df, keys=['AD', 'DP'], parse_alleles=True)
    vcf_df['AF'] = ad_values / dp_values]

    :param vcf_df: [Dataframe representation of the VCF]
    :type vcf_df: [pandas dataframe]
    :param keys: [List of keys to parse from the INFO field of the VCF]
    :type keys: [list]
    :param parse_alleles: [If the dataframe has had mutliallelic sites exploded, and the key value is separated by commas for each allele, this will attempt to correctly identify which value belongs to which allele for a given key in the INFO field. If dataframe does not have allele_index column (generated by the separate_alleles method), then the default allele column will be -1], defaults to False
    :param parse_alleles: bool, optional
    :param coerce: [If key isnt present, coerce will fill in with NaN as the default value, otherwise error], defaults to False
    :param coerce: bool, optional
    :param verbose: [Print things], defaults to True
    :param verbose: bool, optional
    :return: [Tuple of lists of parsed key values, where the tuple is in the order of the keys argument. each list in the tuple is ready to be assigned to the dataframe as a column]
    :rtype: [tuple]
    """
    start = time.time()
    if "allele_index" not in vcf_df.columns:
        vcf_df["allele_index"] = -1
    if parse_alleles:
        temp_series = vcf_df.allele_index.astype(str) + "@" + vcf_df.INFO
    else:
        temp_series = vcf_df.INFO
    annotation_data = zip(*temp_series.apply(partial(_parse_info, keys=keys, parse_alleles=parse_alleles, coerce=coerce)))
    if verbose:
        print("Parsed INFO field for keys {0} in {1} seconds.".format(keys, time.time() - start))
    return annotation_data

def _parse_info(x, keys, parse_alleles=True, coerce=False):
    """[Parses the INFO field in a VCF for specific key values, designed to be used in a vcf_df.INFO.apply(partial(annotate_info, keys=["AD", "ADF", etc..]))]

    :param x: [Row to be filled in by pandas Series APPLY method]
    :type x: [row]
    :param keys: [Elements to be parsed from the info field]
    :type keys: [list]
    :param parse_alleles: [Attempts to use alt_index for slicing data out of comma delimited fields], defaults to True
    :param parse_alleles: bool, optional
    :param coerce: [Force NaNs for unseen values, rather than erroring out], defaults to False
    :param coerce: bool, optional
    :raises KeyError: [If coerce is false and key fails to exist]
    :return: [List of values that have been parsed from the INFO field]
    :rtype: [list]
    """
    if parse_alleles:
        temp_x = x.split("@")
        allele_index = int(temp_x[0])
        x = "@".join(temp_x[1:])
    key_dict = dict([i.split("=") for i in x.split(";") if "=" in i])
    key_values = []
    for key in keys:
        try:
            if parse_alleles and "," in key_dict[key]:
                key_values.append(key_dict[key].split(",")[allele_index])
            else:
                key_values.append(key_dict[key])
        except KeyError:
            if coerce:
                key_values.append(np.nan)
            else:
                raise KeyError("{0} does not exist in row {1}. Please use coerce=True if you want to fill NaN for this value".format(key, x))
    assert len(keys) == len(key_values)
    return key_values

def parse_sample(vcf_df, keys, sample_name, parse_alleles=False, coerce=False, verbose=True):
    """[Parse information from the VCF's SAMPLE_NAME field. An example usage would be:
    ad_values, dp_values = parse_info(vcf_df, keys=['AD', 'DP'], sample_name="HCC1187")
    vcf_df['AF'] = ad_values / dp_values]

    :param vcf_df: [Dataframe representation of the VCF]
    :type vcf_df: [pandas dataframe]
    :param keys: [List of keys to parse from the INFO field of the VCF]
    :type keys: [list]
    :param sample_name: [Name of sample column to parse]
    :param sample_name: str
    :param coerce: [If key isnt present, coerce will fill in with NaN as the default value, otherwise error], defaults to False
    :param coerce: bool, optional
    :param verbose: [Print things], defaults to True
    :param verbose: bool, optional
    :return: [Tuple of lists of parsed key values, where the tuple is in the order of the keys argument. each list in the tuple is ready to be assigned to the dataframe as a column]
    :rtype: [tuple]
    """
    start = time.time()
    if "allele_index" not in vcf_df.columns:
        vcf_df["allele_index"] = -1
    if parse_alleles:
        temp_series = vcf_df.allele_index.astype(str) + "@" + vcf_df.FORMAT.astype(str) + "@" + vcf_df[sample_name]
    else:
        temp_series = vcf_df.FORMAT.astype(str) + "@" + vcf_df[sample_name]
    annotation_data = zip(*temp_series.apply(partial(_parse_sample, keys=keys, parse_alleles=parse_alleles, coerce=coerce)))
    if verbose:
        print("Parsed SAMPLE field for keys {0} in {1} seconds.".format(keys, time.time() - start))
    return annotation_data

def add_af(vcf_df, sample_name, coerce=False):
    vcf_df['FORMAT'] = vcf_df['FORMAT'].str.replace('FA','AF')
    ad, dp = parse_sample(vcf_df, keys=['AD','DP'],sample_name=sample_name, parse_alleles=True, coerce=coerce)
    ad = [float(item) for item in list(ad)]
    dp = [float(item) for item in list(dp)]
    vcf_df['ADDPAF'] = np.array(ad)/np.array(dp)
    raw_af = parse_sample(vcf_df, keys=['AF']     ,sample_name=sample_name, parse_alleles=True, coerce=coerce)
    vcf_df['AF'] = [float(item) for item in list(list(raw_af)[0])]
    return vcf_df

def add_var_type(vcf_df):
    pass


def _parse_sample(x, keys, parse_alleles=True, coerce=False):
    """[Parses the SAMPLE field in a VCF for specific key values, designed to be used in a vcf_df.parse_sample.apply(partial(annotate_info, keys=["AD", "DP", etc..]))]

    :param x: [Row to be filled in by pandas Series APPLY method]
    :type x: [row]
    :param keys: [Elements to be parsed from the SAMPLE field, must exist in format field]
    :type keys: [list]
    :param coerce: [Force NaNs for unseen values, rather than erroring out], defaults to False
    :param coerce: bool, optional
    :raises ValueError: [If coerce is false and key fails to exist]
    :return: [List of values that have been parsed from the SAMPLE field]
    :rtype: [list]
    """
    if parse_alleles:
        temp_x = x.split("@")
        allele_index = int(temp_x[0])
        x = "@".join(temp_x[1:])
    format_keys, data = map(partial(str.split, sep=":"), x.split("@"))
    assert(len(data) == len(format_keys))
    key_values = []
    for key in keys:
        try:
            value = data[format_keys.index(key)]
            if parse_alleles and "," in value:
                key_values.append(value.split(",")[allele_index])
            else:
                key_values.append(value)
        except ValueError:
            if coerce:
                key_values.append(np.nan)
            else:
                raise ValueError("{0} does not exist in row {1}. Please use coerce=True if you want to fill NaN for this value".format(key, x))
    assert len(keys) == len(key_values)
    return key_values

def explode_mnvs(vcf_df):
    output_df = vcf_df.copy(deep=True)
    new_rows = []
    for idx, row in vcf_df.iterrows():
        if ((len(row.REF) == len(row.ALT)) and len(row.REF) > 1):
            # Is a MNV
            for i in range(len(row.ALT)):
                cur_new_row = row.copy(deep=True)
                new_pos = cur_new_row.POS+i
                new_ref = cur_new_row.REF[i]
                new_alt = cur_new_row.ALT[i]
                cur_new_row.POS = new_pos
                cur_new_row.REF = new_ref
                cur_new_row.ALT = new_alt
                new_rows.append(cur_new_row)
    new_rows_df = pd.DataFrame(new_rows)
    return output_df.append(new_rows_df)

def vcf_concordance_sets(vcf_df1, vcf_df2):
    """[Generates concordance sets of variants]

    :param vcf_df1: [vcf_df1]
    :type vcf_df1: [pandas dataframe]
    :param vcf_df2: [vcf_df2]
    :type vcf_df2: [pandas dataframe]
    :return: [set of discordant variants in vcf1, set of discordant variants in vcf2, set of concordant variants]
    :rtype: [list of sets]
    """
    variants_1 = set(vcf_df1[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None))
    variants_2 = set(vcf_df2[["CHROM", "POS", "REF", "ALT"]].itertuples(index=False, name=None))
    concordant = variants_1.intersection(variants_2)
    discordant_1 = variants_1 - concordant
    discordant_2 = variants_2 - concordant
    return discordant_1, discordant_2, concordant

def bcftools_normalize_vcf(vcf_path, output_path, multi_allelic_mode="-both", reference="/reference/science/data/genomes/reference/hs37d5/raw/hs37d5.fa", bcftools_dir="/reference/env/dev/apps/bcftools/1.9"):
    """[Left - align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows.]

    :param vcf_path: [path to vcf]
    :type vcf_path: [str]
    :param output_path: [path to normalized vcf file]
    :type output_path: [str]
    :param multi_allelic_mode: [split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). Can be None for no splitting behavior. If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any. Can be snps|indels|both|any]], defaults to "-both"
    :type multi_allelic_mode: str, optional
    :param reference: [reference fasta], defaults to "/reference/science/data/genomes/reference/hs37d5/raw/hs37d5.fa"
    :type reference: str, optional
    :param bcftools_dir: [bcftools location], defaults to "/reference/env/dev/apps/bcftools/1.9"
    :type bcftools_dir: str, optional
    """
    if multi_allelic_mode:
        multi_allelic_arg = ["-m", multi_allelic_mode]
    else:
        multi_allelic_arg = []

    command = [os.path.join(bcftools_dir, "bcftools"), "norm", "-f", reference, "-O", "v", "-c", "s"] + multi_allelic_arg + [vcf_path]

    # perform normalization
    print("Normalizing {0}".format(vcf_path))
    with open(output_path, 'wb') as out:
        p1 = subprocess.Popen(command, stdout=out, shell=False)
    p1.wait()
    if p1.returncode != 0:
        raise Exception("Normalization failed")

def bedtools_intersect_vcf(vcf_path, bed_file, output_path, bedtools_dir="/reference/env/dev/apps/bedtools/2.27.1/bin"):
    """[Use bedtools to perform intersection of vcf and bed file]

    :param vcf_path: [path to input vcf]
    :type vcf_path: [str]
    :param bed_file: [path to bed file]
    :type bed_file: [str]
    :param output_path: [path to intersected vcf file]
    :type output_path: [str]
    :param bedtools_dir: [bedtools bin dir - should contain intersectBed executable], defaults to "/reference/env/dev/apps/bedtools/2.27.1/bin"
    :param bedtools_dir: str, optional
    :raises IOError: [executable does not exist]
    :raises Exception: [intersection command line failed]
    """
    executable = os.path.join(bedtools_dir, "intersectBed")
    if not os.path.exists(executable):
        raise IOError(executable + "does not exist.")

    # perform intersection
    print("Intersecting {0} with {1}".format(vcf_path, bed_file))
    with open(output_path, 'wb') as out:
        p1 = subprocess.Popen([executable, '-a', vcf_path, '-b', bed_file, '-header'], stdout=out, shell=False)
    p1.wait()
    if p1.returncode != 0:
        raise Exception("Intersection failed")

def intersect_vcf(vcf_df, bed_file):
    """[Intersects VCF dataframe with bed file, returns column of True / False values if intersect for each index in the vcf dataframe]

    :param vcf_df: [Loaded VCF dataframe]
    :type vcf_df: [pandas dataframe]
    :param bed_file: [Bed file path]
    :type bed_file: [str]
    :return: [True or False (intersected or not) values in order of vcf_df.index. can be used to set a column of original dataframe]
    :rtype: [list]
    """
    return list(map(partial(_intersect, bed_dict=_load_bed_file(bed_file)), vcf_df.index))

def _load_bed_file(path):
    df = pd.read_csv(path, sep='\t', header=None, dtype={0: str})
    df = df.iloc[:, 0:3]  # just take chrom, start, and stop positions from bed file
    df.columns = ['chrom', 'start', 'stop']
    df['chrom'] = df['chrom'].str.replace("chromosome", "").replace("chrom", "").replace("chr", "")
    # build chrom / pos data for that specific bed
    bed_dict = defaultdict(list)
    for index, chrom, start, stop in df[['chrom', 'start', 'stop']].itertuples():
        bed_dict[str(chrom)].append((int(start) + 1, int(stop)))  # start +1 because BED files have start 0 based, whereas stop is 1 based.
    for chrom in bed_dict:
        bed_dict[chrom] = sorted(bed_dict[chrom])
    return bed_dict

def _intersect(index, bed_dict={}):
    chrom = str(index[0])
    if chrom not in bed_dict:
        return False
    pos = int(index[1])
    ref_length = len(index[2])
    bindex = bisect.bisect_left(bed_dict[chrom], (pos,))
    if ref_length > 1:
        pos2 = pos + ref_length - 1
        bindex2 = bisect.bisect_left(bed_dict[chrom], (pos2,))
    else:
        bindex2 = None
    if bindex < len(bed_dict[chrom]):
        if (bed_dict[chrom][bindex - 1][0] <= pos <= bed_dict[chrom][bindex - 1][1]) or (bed_dict[chrom][bindex][0] <= pos <= bed_dict[chrom][bindex][1]):
            return True
    else:
        if (bed_dict[chrom][bindex - 1][0] <= pos <= bed_dict[chrom][bindex - 1][1]):
            return True
    if bindex2 and bindex2 < len(bed_dict[chrom]):
        if (bed_dict[chrom][bindex2 - 1][0] <= pos2 <= bed_dict[chrom][bindex2 - 1][1]) or (bed_dict[chrom][bindex2][0] <= pos2 <= bed_dict[chrom][bindex2][1]):
            return True
        elif bindex2 and (bed_dict[chrom][bindex2 - 1][0] <= pos2 <= bed_dict[chrom][bindex2 - 1][1]):
            return True
    return False
