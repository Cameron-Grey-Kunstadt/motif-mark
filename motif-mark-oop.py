# Cameron Kunstadt
# UO import debugpy, platform
# 2/27/2025

#TODO: function to find exons

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta",)  
parser.add_argument('-m', "--motif")
args = parser.parse_args()


class FastaRecord:
    def __init__(self, header, seq, unique_motif_list):
        self.header = header
        self.seq = seq
        self.length = len(seq)
        self.gene = self.parse_gene(self, header)
        self.needs_revcomp = self.parse_needs_revcomp(self, header)

        if self.needs_revcomp == True:
            self.seq = self.reverse_compliment(self, seq)

        self.unique_motif_list = unique_motif_list
        self.motif_object_list = []
        self.exon_start, self.exon_stop = self.get_exon_range(seq)


        

    def reverse_compliment(self, seq):
        compliments = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'G': 'C',
                       'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'g': 'c',
                       'N':'N', 'n':'n'}
        rev_seq = ""

        for base in seq:
            rev_seq = compliments[base] + rev_seq

        return rev_seq
    
    def parse_gene(self, header):
        gene_match = re.search(">[A-Z0-9]+", header)
        gene = gene_match.group()
        return gene[1:]
    
    def parse_needs_revcomp(self, header):
        rev_comp = re.search("reverse complement", header)
        
        try:
            if rev_comp.group():
                return True
        except AttributeError:
            return False
        
    def get_exon_range(self, seq):

        capital_letters = list(re.finditer(r'[A-Z]', seq))

        exon_start = capital_letters[0].start()
        exon_stop = capital_letters[-1].start()

        return exon_start, exon_stop
    
    def find_motifs(self, seq):
        



        

# class MotifClass:
#     def __init__(self, motif, color, seq_length, exon_start, exon_stop):
#         self.motif = motif
#         self.color = color
#         self.seq_length = seq_length
#         self.exon_start = exon_start
#         self.exon_stop = exon_stop

#     def draw_line(self):
#         return None



# unique_motif_list = []

# with open(args.motif, 'r') as motif_file:
#     for line in motif_file:
#         unique_motif_list.append(line.strip())

# with open(args.fasta, 'r') as fasta:
#     while True:
#         header = fasta.readline()
#         seq = fasta.readline()

#         record = FastaRecord(header, seq, unique_motif_list)



