"""
Naive SNV caller from pileups:

Given a reference sequence and a read pileup, make one of the following variant calls at each position in the reference:
1. "nocall": not enough evidence to make a call
2. "ref": only the reference allele is present
3. "[ACGT]hom: an Alt allele (A,C,G, or T) is present in >75% of reads
4. "[ACGT]het: an Alt allele (A,C,G, or T) is present in 25%--75% of reads
5. "[ACGT]low: an Alt allele (A,C,G, or T) is present in <25% of reads

If more than 1 Alt allele is present, the variant call should contain all of them, comma delimited

The call_variants function takes a min_depth argument; if the total read depth is less than this, 
the position is "nocall"
"""

## Imports
from __future__ import print_function
from collections import Counter

## Inputs
test_ref = 'TTTAGAGCGC'
test_pileup = ['....,...,,,,.,...,',
          ',..,,,.c..C.,,,',
          'AAaAAaaaA.,AaaaAa',
          '.C.,.g,gG.G...,G.',
          '..,C,..C,,CCCC.C,c',
          '.,..tT,,.t',
          '',
          '.Taat..,.,...,',
          '.cc',
          'aA.,,A.,...,a.,Aa...,']


## Variant Caller
def call_variants(pileup, ref, min_depth):
    # TODO: Make a variant call at each position on ref, using the pileup
    call_set = []
    call = ""
    for idx, item in enumerate(pileup):
        ## If depth is lesser than the threshold, no need to process this position
        if len(item) < min_depth:
            call = "nocall"
            call_set.append(call)
            continue
        ## Since in this problem upper and lower alphabets mean the same,
        ## Simplifying by converting all lowercase to uppercase and
        ## Simplifying the problem by converting . to , since strand specificity is not needed
        item = item.upper()
        item = item.replace('.',',')
        
        ## Create counter object from the pileup item 
        pileup_counter = Counter(item)
        if ',' in pileup_counter and len(pileup_counter) == 1:
            call = "ref"
            call_set.append(call)
            continue
        
        ## Process each pileup item
        pileup_call = {}
        for k, v in pileup_counter.items():
            ## Do not need to process matches with ref
            if k == ',':
                continue
                
            ## Calculate percentage and assign suffix based on the percentage value
            percent = v*100/len(item)
            if percent > 75:
                suffix = "hom"
            if percent >= 25 and percent <= 75:
                suffix = "het"
            if percent < 25:
                suffix = "low"
            ## Added space for readability while joining
            pileup_call[k] = k+" "+suffix
        call_set.append(','.join(pileup_call.values()))
    return call_set

def main():
    call_set = call_variants(test_pileup, test_ref, 5)

    result = [print(str(idx+1)+"\t"+str(item)) for idx, item in enumerate(call_set)]


if __name__ == "__main__":
    main()