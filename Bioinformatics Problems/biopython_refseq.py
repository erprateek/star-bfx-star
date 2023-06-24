import os
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#Downloaded the sequence from RefSeq by following the CDS location in the list of features/
#https://www.ncbi.nlm.nih.gov/nuccore/NM_000022.2?from=129&to=1220&report=fasta
#Save it to the same directory as this file
if len(sys.argv) < 2:
    print "Usage: python Ques2.py filename"
    print "Please enter the file name too and try again"
    sys.exit(1)
else:
    filename = sys.argv[1]

#Reading in the sequence file
cds=open(filename,'r')
mutation_description = 'c.239A>G'
position = 239
ref_allele = 'A'
alt_allele = 'G'
cds_seq = cds.readlines()

#Parsing the CDS sequence file
header = cds_seq[0]
sequence_string_list = cds_seq[1:]
sequence_string_list = [line.strip() for line in sequence_string_list]
sequence = ''.join(sequence_string_list)
if not sequence.startswith('ATG'):
    print "Sequence does not start with the start codon i.e. not CDS"
    sys.exit(1)

#Generating the mutated sequence
mutated_sequence = list(sequence)
original_allele = mutated_sequence[position-1]
if original_allele != ref_allele:
    print "Ref allele does not match the allele at the sequence position."
    sys.exit(1)

mutated_sequence[position-1] = alt_allele
mutated_sequence = ''.join(mutated_sequence)

#Translation using Biopython
coding_dna = Seq(mutated_sequence, generic_dna)
original_protein = Seq(sequence, generic_dna).translate()
mutated_protein = coding_dna.translate()

print "Original protein: "+original_protein
print "Mutated protein: "+mutated_protein
