from __future__ import print_function
from Bio import Entrez
from Bio import SeqIO
import re
import os

def fetch_sequence(refseqid):
    """Communicates through the NCBI Entrez API to fetch the id corresponding 
    to the refseq id so as the fetch the CDS for the purpose of the task"""
    Entrez.email = "prateektandon@cmu.edu"
    handle = Entrez.esearch(db="nucleotide", term=refseqid)
    info = Entrez.read(handle)
    # checks if it succeeds
    if not info:
        print("Fetching did not succeed")
        return None
    seq_id = info['IdList'][0]
    handle = Entrez.efetch(db="nuccore", id=seq_id, rettype="fasta_cds_na", retmode="text")
    # check if it succeeds
    record = handle.read()
    if not record:
        print("Could not fetch the CDS")
        return None
    # writing the fasta file to be easily parsed while reading by the SeqIO module
    fname = refseqid+".fasta"
    out_handle = open(fname, 'w')
    out_handle.write(record.rstrip('\n'))
    out_handle.close()
    # check if the file was successfully created
    if not os.path.exists(fname):
        print("File not created")
        return None
    return fname

def read_fetched(fname):
    # isolate this part so that it can be used without API communication
    seqs = SeqIO.parse(open(fname),'fasta')
    sequence = []
    for seq in seqs:
        sequence.append(seq.seq)
    if not len(sequence) == 1:
        print("The CDS was not obtained successfully. Please check the generated file.")
        return None
    return sequence[0]

def mutation_info(mutation):
    """Parses the mutation pattern provided and segregates the 
    original allele, mutated allele and the position at which the
    mutation takes place.
    Original allele is required only for the verification that the 
    position actually has the allele mentioned in the position.
    The position returned is 0 based index so that it can be directly
    used in the sequence indexing.
    """
    if not is_valid_pattern(mutation):
        print("Please check the pattern")
        return None
    #Finds the position by getting the numeber in the pattern
    loc = re.findall(r'\d+', mutation)[0]   
    loc_int = int(loc)
    #index based identification of the original allele index
    orig_allele = mutation[mutation.index(loc)+len(loc)]
    mutated_allele = mutation[-1]
    return orig_allele,mutated_allele,(loc_int-1)

def is_valid_pattern(pattern):
    """ Checks for the validity of the pattern provided """
    if not pattern.startswith("c."):
        return False
    if not ">" in pattern:
        return False
    arrow_index = pattern.index(">")
    orig_index = arrow_index-1
    mutated_index = arrow_index+1
    alphabet = "ATCG"
    if pattern[orig_index] not in alphabet or pattern[mutated_index] not in alphabet:
        return False
    locs = re.findall(r'\d+', pattern)
    if not len(locs) == 1:
        return False
    return True

def main():
    mutation = 'c.239A>G'
    refseqid = "NM_000022.2"
    # If internet communication needs to be omitted, change the following line to be 
    # fname = "Name.fasta" where Name is the actual filename
    # Rest everything should work just fine
    fname = fetch_sequence(refseqid)
    cds = read_fetched(fname)
    wildtype,mutatedtype,pos = mutation_info(mutation)
    # mutating the sequence in the required position
    cds_mutated = cds.tomutable()
    if cds_mutated[pos] == wildtype:
        cds_mutated[pos] = mutatedtype
    else:
        print("The original allele does not occur at the position specified")
    cds_mutated = cds_mutated.toseq()
    cds_protein = str(cds_mutated.translate())
    print(cds_protein)

if __name__ == '__main__':
    main()
