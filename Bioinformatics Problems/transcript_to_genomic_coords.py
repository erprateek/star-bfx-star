# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import glob
import logging
import ntpath
import re

from cigar import CIGAR

__author__ = 'Prateek Tandon'
__email__  = 'prateektandon@alumni.cmu.edu'

logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

## DONE Error handling for querying beyond CIGAR - Script returns -1
## DONE Error handling for querying unknown transcripts - Script returns -1 with empty string for Chromosome

def create_cigar_obj_for_row(row):
    """Helper function that generates a CIGAR object for a row of the dataframe

    Args:
        row (pandas.Series): A row of the the input 4 column transcripts dataframe

    Returns:
        CIGAR: A CIGAR object used for generating genomic coordinates
    """
    start = int(row['Pos'])
    cigar_string = row['Cigar']
    return CIGAR(cigar_string, start=start)

def qpos_to_gpos(row):
    """Helper method that processes a row of the query dataframe using a pre-processed row of the input transcript dataframe. 
    It uses prepopulated CIGAR object for generating genomic coordinates for the query position. 
    If the query position cannot be found with the transcript coordinates, returns -1

    Args:
        row (pandas.Series): A row of the 2 column input query dataframe

    Returns:
        integer : The genomic coordinate corresponding to the query position on the given transcript
    """
    logger.debug(f"Currently processing row: {row}")
    transcript_id = row['Transcript']
    qpos = row['Pos']
    transcript_row = input_base_df.loc[input_base_df['Transcript'] == transcript_id]
    if transcript_row.shape[0] > 0:
        logger.debug("Successfully retrieved transcript info.")
    else:
        logger.warning(f"Did not find any transcripts by the name: {transcript_id}")
        return -1
    cigar_obj = transcript_row['cigar_obj'].iloc[0]
    curr = cigar_obj.start
    curt = 0
    for idx, item in enumerate(cigar_obj.cigar_scorer):
        # Use the scoring matrix from the CIGAR object to advance/not advance the position on ref/transcript
        # per the CIGAR op one by one
        transcript_score = item[0]
        ref_score = item[1]
        logger.debug(f"Current ref position: {curr}, current transcript position: {curt}")
        if curt == qpos:
            logger.debug(f"For transcript position {curt}, the corresponding genomic position is: {curr}")
            return curr
        curr = curr+ref_score
        curt = curt+transcript_score
    return -1

def f_reader_df(input_fname, cols=None):
    """File reader function to isolate each component of the task

    Args:
        input_fname (string): Path to the file to be read
        cols (list, optional): List of column names to be used for parsing the file. Defaults to None.

    Returns:
        pandas.DataFrame: File read as a pandas dataframe
    """
    df = pd.DataFrame()
    try:
        df = pd.read_csv(input_fname, sep="\t", header=None)
    except:
        logger.error("A problem occurred while reading the file into a dataframe.")
    if cols:
        df.columns = cols
    return df

def generate_output_file(df, output_fname='output.txt'):
    """Helper function to help write the output dataframe to a file. 
    Encapsulates the formatting of the output file per the requirement. 
    i.e. No header, No column names etc.

    Args:
        df (pandas.DataFrame): Output dataframe to be written to a file.
        output_fname (str, optional): Output filename. Defaults to 'output.txt'.
    """
    output_df = df
    output_df['Pos']  = output_df['Pos'].astype(int)
    output_df['gPos'] = output_df['gPos'].astype(int)
    try:
        df.to_csv(output_fname, sep="\t", index=False, header=None)
    except:
        logger.error("Data generation successful but a problem occurred while writing results to file.")
        sys.exit(1)
    logger.info(f"Output file: {output_fname} was created successfully.")

def check_args(num_args=3):
    """Helper function to check whether valid number of arguments have been passed to the script.
    Modularized to take in changes to requirements if they come in the future.

    Args:
        num_args (int, optional): Number of arguments to use for validation. Defaults to 3.
    """
    if len(sys.argv) != num_args:
        logger.error("Usage: python transcript_to_genomic_coords.py <TRANSCRIPTS_FILE> <QUERY_FILE>")
        logger.error(f"Received: {len(sys.argv[1:])} arguments instead. Exiting!")
        sys.exit(1)
    logger.debug(f"Verified number of arguments passed == {num_args}")

def get_chromosome(row):
    """Gets the chromosome that the transcript resides on. If not found, returns empty string. 

    Args:
        row (pandas.Series): Row of input 2 column query dataframe

    Returns:
        str: Chromosome that the transcript corresponds to. If not found, returns empty string.
    """
    transcript_rows = input_base_df[input_base_df['Transcript'] == row['Transcript']]
    if transcript_rows.shape[0] > 0:
        return transcript_rows['Chromosome'].tolist()[0]
    return ""

def main():
    check_args(num_args=3)
    global input_base_df
    input_base_file            = sys.argv[1]
    input_q_file               = sys.argv[2]
    output_fname               = "output_from_script.txt"
    input_base_df              = f_reader_df(input_base_file,cols=['Transcript','Chromosome','Pos','Cigar'])
    input_q_df                 = f_reader_df(input_q_file   ,cols=['Transcript','Pos'])
    input_base_df['cigar_obj'] = input_base_df.apply(lambda x: create_cigar_obj_for_row(x), axis=1)
    input_q_df['Chromosome']   = input_q_df.apply(get_chromosome, axis =1)
    input_q_df['gPos']         = input_q_df.apply(qpos_to_gpos, axis=1)
    generate_output_file(input_q_df, output_fname=output_fname)

if __name__ == '__main__':
    main()

