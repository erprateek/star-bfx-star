import re

class CIGAR(object):
    cigar_regex = r'\d+[IDMPX]+'
    pattern = re.compile(cigar_regex)
    cigar_dict = {    
        'M': [1, 1],    # Represents a MATCH in sequence
        'I': [1, 0],    # Represents an insertion 
        'D': [0, 1],    # Represents a deletion 
        'N': [0, 1],    # Represents Reference Skip
        'S': [1, 0],    # Represents Soft Clip
        'H': [0, 0],    # Represents Hard Clip
        'P': [0, 0],    # Represents padding. Mostly for denovo assemblers. Does not advance positions.
        '=': [1, 1],    # Represents more granular version of MATCH
        'X': [1, 1],    # Represents mismatch. M segments can be decomposed into = and X segments
    }
    def __init__(self, cstring, start=0):
        """Constructor for CIGAR class. Takes in the input CIGAR string, separates individual CIGAR
        operations and then creates a expanded score list for the sequence based on the CIGAR.
        Primarily used for querying reference positions corresponding to given transcript positions 

        Args:
            cstring ([type]): [description]
            start (int, optional): [description]. Defaults to 0.
        """
        self.cstring = cstring
        self.cparsed = self.separate_cigar()
        self.start = start
        self.cigar_scorer = self.expand_cigar()
        
    def separate_cigar(self):
        """Helper method to separate all the CIGAR operations into individual operations
        Example: 4M3I6D -> ['4M','3I','6D']

        Returns:
            list: separated CIGAR operations in the form of a list. List is in the original order as input.
        """
        cigar_parsed = re.findall(self.pattern, self.cstring)
        return cigar_parsed
    
    def expand_cigar(self):
        """Helper public method to generate a per position score 
        for transcript and reference based on the CIGAR Operation 
        Example: 4M -> [[1,1],[1,1],[1,1],[1,1]]

        Returns:
            list: List containing an expanded score version of the CIGAR string
        """
        cigar_dict = self.cigar_dict
        ref_q_cigar_seq_incrementer = []
        for item in self.cparsed:
            clen = int(re.findall(r'\d+',item)[0])
            cop = re.findall(r'[A-Z]+', item)[0]
            score_for_op = cigar_dict[cop]
            cur_interval_scores = [score_for_op for i in range(clen)]
            ref_q_cigar_seq_incrementer+=cur_interval_scores
        return ref_q_cigar_seq_incrementer