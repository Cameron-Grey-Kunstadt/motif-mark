import argparse
import re, os
import bioinfo
import cairo
import random
# Cameron Kunstadt
# UO BGMP
# 2/27/2025

#TODO: 
# - deal with the weird ys in the unique motifs list
# - basically scale nt to pixels
# - draw exons to scale
# - add labels and key
# - find better ways to add color
# - ignore case in motif sequences

#TODONE:
# - ignore case in motif sequences, DONE
# - ignore reverse compliment, DONE
# - allow for the printout of multiple record graphs DONE, may need tweaks
# - Convert all Us to Ts before doing anything, DONE, need to validate

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta",)  
parser.add_argument('-m', "--motif")
args = parser.parse_args()


class FastaRecord:
    '''A FastaRecord object serves to parse and hold useful information from a FASTA header
    and a FASTA sequence, and to find motifs within the sequence, based on the given unique_motif_list.
    Motifs are create as Motif objects, and saved in a list.'''
    def __init__(self, header, seq, unique_motif_list, record_num):
        self.header = header
        self.seq = self.convert_U_to_T(seq)
        self.length = len(seq)
        self.gene = self.parse_gene(header)
        self.unique_motif_list = unique_motif_list
        self.motif_object_list = self.find_motifs(self.seq, self.unique_motif_list, record_num)
        self.exon_start, self.exon_stop = self.get_exon_range(seq)
        self.record_num = record_num

    
    def parse_gene(self, header):
        '''Pulls gene out of fasta header, assumes its the first letters & numbers
        until the first space'''
        gene_match = re.search(">[A-Z0-9]+", header)
        gene = gene_match.group()
        return gene[1:]
    
    def convert_U_to_T(self, seq):
        '''Converts all Us in the file to Ts'''
        for i in range(0, len(seq)):
            if seq[i] == "U": seq[i] = "T"
            elif seq[i] == "u": seq[i] = "t"
        return seq
    
        
    def get_exon_range(self, seq):
        '''Gets the exon from the given sequence provided it has one, determines
        this by getting the first and last capital letter in the sequence'''
        capital_letters = list(re.finditer(r'[A-Z]', seq))
        exon_start = capital_letters[0].start()
        exon_stop = capital_letters[-1].start()

        return exon_start, exon_stop
    
    def build_regex(self,motif):
        motif = motif.lower()
        list_motif = list(motif)
        '''Returns regex of given motif, will handle Ys accordingly, ignores case'''
        for i, base in enumerate(list_motif):
            if base == "y":
                list_motif[i] = "[ct]"

        str_motif = ''.join(list_motif)
        print(str_motif)
        return str_motif
    
    def find_motifs(self, seq, unique_motif_list, record_num):
        '''Searches for all the instances of motifs in the given sequence
        and creates motif objects for each of them, case insensitive'''
        seq = seq.lower()
        motif_objects_list = []

        for unique_motif in unique_motif_list:
            unique_motif = unique_motif.lower()
            motif_instances = list(re.finditer(self.build_regex(unique_motif), seq))
            for motif_instance in motif_instances:
                motif_span = motif_instance.span()
                motif_start = motif_span[0]
                motif_stop = motif_span[1]

                motif_object = Motif(unique_motif, unique_motif_color_dict[unique_motif], len(seq), motif_start, motif_stop, record_num)
                motif_objects_list.append(motif_object)

        return motif_objects_list
    
    def draw_graph(self, record_num):
        # Set main line as black, stroke
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(25,(175 + (350 * record_num)))
        ctx.line_to(1375, (175 + (350 * record_num)))
        ctx.set_line_width(5)

        ctx.stroke()

        # Set Exon, fill, stroke
        ctx.set_source_rgb(0.4, 0.9, 0.4)
        ctx.rectangle(200, (150 + (350 * record_num)), 300, 50)

        ctx.fill()
        ctx.stroke()
        return 0


class Motif:
    '''A Motif object holds important information about the motifs found in a given
    FASTA record, and allows for easy drawing on a pycairo context.'''
    def __init__(self, motif, color, seq_length, motif_start, motif_stop, record_num):
        self.motif = motif
        self.color = color
        self.seq_length = seq_length
        self.motif_start = motif_start
        self.motif_stop = motif_stop
        self.record_num = record_num

    def draw_line(self, x, y, record_num):
        r, g, b = self.color

        ctx.set_source_rgb(r/255, g/255, b/255)
        ctx.rectangle(x, (y + (record_num * 350)), 5, 50)
        ctx.fill()
        ctx.stroke()

# I will probably delete this, this isn't my favorite
def random_color_generator():
    r = random.randint(0, 255)
    g = random.randint(0, 0)
    b = random.randint(0, 255)
    return (r, g, b)

def convert_bp_to_pixels(basepairs):
    return None

def get_length_of_file(filepath):
    with open("temp_oneline.fasta", 'r') as file:
        return len(file.readlines())

def determine_length_of_surface(line_num):
    return ((line_num / 2) * 350) + 175


# Create unique list of motifs, assign each a unique random color and put into dictionary
unique_motif_list = []
unique_motif_color_dict = {}

with open(args.motif, 'r') as motif_file:
    for line in motif_file: unique_motif_list.append(line.strip().lower())

for unique_motif in unique_motif_list: unique_motif_color_dict[unique_motif] = random_color_generator()


# Not assuming the incoming file is oneline, so a temp one-line fasta file is created, we also get the number of lines to define the length of the image
bioinfo.oneline_fasta(args.fasta, "temp_oneline.fasta")
file_numlines = get_length_of_file("temp_oneline.fasta")

# Set surface, paint white
surface = cairo.ImageSurface(cairo.FORMAT_RGB24, 1400, int(determine_length_of_surface(file_numlines)))
ctx = cairo.Context(surface)
ctx.set_source_rgb(255, 255, 255)
ctx.paint()


with open("temp_oneline.fasta", 'r') as fasta:
    record_num = 0
    while True:
        header = fasta.readline()
        if header == "":
            break
        seq = fasta.readline()

        record = FastaRecord(header, seq, unique_motif_list, record_num)
        record.draw_graph(record.record_num)

        for motif in record.motif_object_list:
            motif.draw_line(motif.motif_start,150, motif.record_num)
        
        record_num += 1


os.remove('temp_oneline.fasta')
surface.write_to_png('test_rect_and_line.png')


# def build_regex(motif):
#         motif = motif.lower()
#         list_motif = list(motif)
#         '''Returns regex of given motif, will handle Ys accordingly, ignores case'''
#         for i, base in enumerate(list_motif):
#             if base == "y":
#                 list_motif[i] = "[ct]"

#         return ''.join(list_motif)

# x = build_regex("atcyatcyyaaa")
# print(x)
    