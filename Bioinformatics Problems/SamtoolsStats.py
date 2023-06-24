import os
import json
import argparse

class SamtoolsStats:
    def __init__(self, path, verbose=False):
        self.path = path

    @staticmethod
    def parse_samtools_stats(stats_file):
        output_dict = {}
        with open(stats_file) as cur_file:
            for line in cur_file:
                line = line.strip()
                if "SN" in line:
                    if "reads mapped:" in line:
                        mapped_reads = line.split(":")[1].strip()
                        output_dict['reads_mapped'] = mapped_reads
                    elif "reads mapped and paired" in line:
                        mapped_and_paired = line.split(":")[1].strip().split("\t")[0]
                        output_dict['reads_mapped_and_paired'] = mapped_and_paired
                    elif "reads unmapped" in line:
                        unmapped_reads = line.split(":")[1].strip()
                        output_dict['unmapped_reads'] = unmapped_reads
        return output_dict
    
    def print_stats(self):
        output_dict = self.parse_samtools_stats(self.path)
        print(json.dumps(output_dict, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--stats_file",type=str,required=True, help="Output of samtools stats")
    args = parser.parse_args()
    stats_file = SamtoolsStats(args.stats_file)
    stats_file.print_stats()
