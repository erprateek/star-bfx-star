# utils and such
import logging
import sys
import pandas as pd
import sys
import ntpath
import glob
import os
import re
from os import access, R_OK
from os.path import isfile
from functools import partial
from functools import reduce

def setup_logger(log_file=None, debug=False):
    logger = logging.getLogger()
    logger.propagate = False
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # handlers
    if not logger.handlers:
        st_handler = logging.StreamHandler()
        log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        st_handler.setFormatter(log_format)
        logger.addHandler(st_handler)

    # optional logging to file                                                                                                                                                                                                              
    if log_file:
        cwd = os.getcwd()
        mkdir(os.path.join(cwd,'logs'))
        f_handler = logging.FileHandler(os.path.join(cwd, 'logs',log_file+".log"))
        log_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        f_handler.setFormatter(log_format)
        logger.addHandler(f_handler)
    return logger

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

def get_best_dims_for_subplots(datasize):
    factor_list = factors(datasize)
    factor_list = [item for item in factor_list if item > 1 and item < datasize]
    if len(factor_list) == 2:
        return (int(max(factor_list)), int(min(factor_list)))
    else:
        ### max columns is 5
        max_fact_lteq5  = max([item for item in factor_list if item <= 5])
        cols = int(max_fact_lteq5)
        rows = int(datasize/cols)
        return (rows, cols)

def is_valid_file(input_file="", num_cols=3, sep="\t"):
    # Check if the path is valid
    # Check if the file is readable
    # Check if the number of columns is as expected
    # Check if the separator in the file is as supplied
    # Check if the number of rows correspond with the wc -l output from bash
    if not isfile(input_file) or not access(input_file, R_OK):
        print("File {} doesn't exist or isn't readable".format(input_file))
        return False
    file_df = pd.read_csv(input_file, sep=sep)
    nrows, ncols = file_df.shape
    if os.stat(input_file).st_size == 0:
        print(input_file+" File size is zero, Exiting!")
        return False
    if ncols != num_cols:
        print("Number of columns after splitting and using delimiter="+sep+" equals - "+str(ncols)+" expected "+str(num_cols))
        return False
    return True

def get_unique_file_name(path):
    if not os.path.isfile(path):
        return path
    else:
        filebasename = ntpath.basename(path)
        filepath = ntpath.dirname(path)
        filename, extension = os.path.splitext(filebasename)
        fname_split = filename.split("_")
        fno = 0
        new_filebasename = filename
        ## If there is a digit, increment and append
        ## if /a/b/c/utils_1.py -> /a/b/c/utils_2.py
        if fname_split[-1].isdigit():
            fno = int(fname_split[-1])
            new_filebasename = "_".join(fname_split[:-1])
        ## if there is no digit, append
        ## if /a/b/c/utils.py -> /a/b/c/utils_1.py
        fno_incremented = fno+1
        new_fname_no_ext = "_".join([new_filebasename, str(fno_incremented)])
        ## reconstruct new path
        new_filename = "".join([new_fname_no_ext, extension])
        new_path = os.path.join(filepath, new_filename)
        unique_path = get_unique_file_name(new_path)
    return unique_path
        

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        pass
        #print("Directory exists. Skipping creation")

def extract_float_from_string(s):
    '''Helper method that returns the first float element that is extractable from the string'''
    #print("Extracting float from string")
    float_from_string = re.findall(r"\d+\.\d+", s)
    if float_from_string:
        print("Found float in string: %s" %(s))
        return float_from_string[0]
    #print("Could not find any float in the string %s"%(s))
    return []

def is_empty_list(x):
    #print("Checking if the list is empty")
    if len(x) == 0:
        return True
    elif len(x) != 0:
        print("Type of list: %s"%(type(x[0])))
        if isinstance(x[0], list):
            print("This item has sublist")
            flat_list = [item for sublist in x for item in sublist]
            flat_list = [item for item in flat_list if item != ""]
            print("Flattened list: %s"%(flat_list))
            return len(flat_list) == 0
        else:
            print("Removing unwanted elements: empty strings and None elements to see if the list is still empty")
            clean_list = [item for item in x if item != "" and item is not None]
            return len(clean_list) == 0

def merge_dataframes(df_list):
    """[Takes in a list of dataframes and return a combined dataframe out of it.]
    
    Arguments:
        df_list {[list]} -- [List containing dataframes containing same structure so that they are easily concatenated]
    
    Returns:
        [pd.DataFrame] -- [Dataframe concatenated from the list]
    """
    result_df = pd.DataFrame()
    for item in df_list:
        result_df = result_df.append(item, ignore_index=True)
    return result_df

