import re
from collections import Counter

class Pileup(object):
    def __init__(self, pileup_string):
        self.pileup_string = pileup_string

    def normalize_pileup(self):
        pileup_string = self.pileup_string
        pileup_normalized = pileup_string.upper() \
                             .replace(',','.') \
                             .replace('<','>')
        return pileup_normalized

    def get_indel_from_pileup_string(self, s):
        pattern = re.compile(r'[+-]\d+[ATGCatgc]*')
        indel_patterns = re.findall(pattern, s)
        #unique_indel_patterns = [item.upper() for item in indel_patterns]
        return Counter(indel_patterns)

    def parse_pileup(self):
        s = self.normalize_pileup()
        indels_counter = self.get_indel_from_pileup_string(s)
        ps_copy = s
        for i in indels_counter:
            ps_copy = ps_copy.replace(i,'')
        remaining_counts = Counter(ps_copy)
        combined = dict(list(indels_counter.items()) + list(remaining_counts.items()))
        return combined, remaining_counts