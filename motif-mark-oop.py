import argparse
import re, os
import bioinfo
import cairo
import itertools
# Cameron Kunstadt
# UO BGMP
# 2/27/2025

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta",)  
parser.add_argument('-m', "--motif")
args = parser.parse_args()

# GLOBAL
screen_width = 2000
margin_width = 25
color_palette = [
    (230, 0, 73),
    (11, 180, 255),
    (80, 233, 145),
    (230, 216, 0),
    (155, 25, 245),
    (255, 163, 0),
    (220, 10, 180),
    (179, 212, 255),
    (0, 191, 160)
]

class FastaRecord:
    '''A FastaRecord object serves to parse and hold useful information from a FASTA header
    and a FASTA sequence, and to find motifs within the sequence, based on the given unique_motif_list.
    Motifs are create as MotifInstance objects, and saved in a list.'''
    def __init__(self, header, seq, unique_motif_list, record_num):
        self.header = header
        self.seq = self.convert_U_to_T(seq)
        self.gene_object = Gene(header, record_num)
        self.unique_motif_list = unique_motif_list
        self.motif_object_list = self.find_motifs(self.seq, self.unique_motif_list, record_num)
        self.multiplier = (screen_width - (2 * margin_width)) / len(seq)
        self.record_num = record_num
        self.exon_object = Exon(self.seq, self.multiplier, self.record_num)
    
    def convert_U_to_T(self, seq):
        '''Converts all Us in the file to Ts'''
        seq = seq.replace("U", "T")
        seq = seq.replace("u", "t")
        return seq
    
    def build_regex(self,motif):
        '''Returns regex of given motif, will handle Ys accordingly, ignores case'''
        motif = motif.lower()
        motif = motif.replace("y", "[ct]")
        return motif
    
    def find_motifs(self, seq, unique_motif_list, record_num):
        '''Searches for all the instances of motifs in the given sequence
        and creates motif objects for each of them, case insensitive'''
        seq = seq.lower()
        motif_objects_list = []

        for unique_motif_object in unique_motif_list:
            unique_motif = unique_motif_object.motif.lower()
            unique_motif = self.convert_U_to_T(unique_motif)
            motif_instances = list(re.finditer(self.build_regex(unique_motif), seq))
            for motif_instance in motif_instances:
                motif_span = motif_instance.span()
                motif_start = motif_span[0]
                motif_stop = motif_span[1]

                motif_object = MotifInstance(unique_motif, unique_motif_object.color, len(seq), motif_start, motif_stop, record_num)
                motif_objects_list.append(motif_object)

        return motif_objects_list
    
    def draw_graph(self):
        '''Draws horizontal line to represent sequence'''
        # Set main line as black, stroke
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(margin_width,(175 + (350 * self.record_num)))
        ctx.line_to(screen_width-25, (175 + (350 * self.record_num)))
        ctx.set_line_width(5)
        ctx.stroke()


class Exon:
    def __init__(self, seq, multiplier, record_num):
        self.exon_start, self.exon_stop = self.get_exon_range(seq)
        self.exon_length = self.exon_stop - self.exon_start
        self.multiplier = multiplier
        self.record_num = record_num

    def get_exon_range(self, seq):
        '''Gets the exon from the given sequence provided it has one, determines
        this by getting the first and last capital letter in the sequence'''
        capital_letters = list(re.finditer(r'[A-Z]', seq))
        exon_start = capital_letters[0].start()
        exon_stop = capital_letters[-1].start()

        return exon_start, exon_stop

    def draw_exon(self):
        '''Draws larger rectangle on the horizontal sequence line to represent exon'''
        # Set Exon, fill, stroke
        r, g, b = color_palette[8]
        ctx.set_source_rgb(r/255, g/255, b/255)
        ctx.rectangle((self.exon_start * self.multiplier) + margin_width, (150 + (350 * self.record_num)), (self.exon_length * self.multiplier), 50)
        ctx.fill()
        ctx.stroke()

class Gene:
    def __init__(self, header, record_num):
        self.gene_name = self.parse_gene(header)
        self.record_num = record_num
    
    def parse_gene(self, header):
        '''Pulls gene out of fasta header, assumes its the first letters & numbers
        until the first space'''
        gene_match = re.search(">[A-Z0-9]+", header)
        gene = gene_match.group()
        return gene[1:]

    def write_gene(self):
        '''Writes the gene name in appropriate place over corresponding gene graph'''
        # Color, size, font, position, text
        ctx.set_source_rgb(0, 0, 0) 
        ctx.set_font_size(40) 
        ctx.select_font_face("Papyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
        ctx.move_to(50, 50 + (self.record_num * 350)) 
        ctx.show_text(self.gene_name) 
        ctx.stroke() 


class UniqueMotif:
    '''A Motif object holds important information about a given motif, and allows for writing'''
    def __init__(self, motif, motif_num, color):
        self.motif = motif
        self.motif_num = motif_num
        self.color = color

    def write_key(self):
        '''Writes the gene name for the key in the bottom right of the graph'''
        ctx.set_source_rgb(0, 0, 0) 
        ctx.set_font_size(30) 
        ctx.select_font_face("Papyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL) 
        ctx.move_to(50, (screen_height - 250) + (self.motif_num * 50)) 
        ctx.show_text(self.motif) 
        ctx.stroke()

    def draw_color_square(self):
        '''Draws a color square corresponding to the motif, for the key'''
        r, g, b = self.color

        ctx.set_source_rgb(r/255, g/255, b/255)
        ctx.rectangle(250 , (screen_height - 260) + (self.motif_num * 50), 25, 25)
        ctx.fill()
        ctx.stroke()

class MotifInstance(UniqueMotif):
    '''A MotifInstance inherates the Motif class, but holds additional important information about the motifs found in a given
    FASTA record, and allows for easy drawing on a pycairo context.'''
    def __init__(self, motif, color, seq_length, motif_start, motif_stop, record_num):
        UniqueMotif.__init__(self, motif, 666, color)
        self.seq_length = seq_length
        self.motif_start = motif_start
        self.motif_stop = motif_stop
        self.record_num = record_num

    def draw_line(self,y):
        '''Draws small vertical line at corresponding points in the sequence graph to represent a found motif'''
        r, g, b = self.color
        ctx.set_source_rgb(r/255, g/255, b/255)
        ctx.rectangle(convert_bp_to_pixels(self.motif_start, self.seq_length), (y + (self.record_num * 350)), 5, 50) # record num adjusts which graph it needs to be placed on, on the page
        ctx.fill()
        ctx.stroke() 



def convert_bp_to_pixels(bp_loc, seq_length):
    '''Based on screen with, size of the margin, and sequence length, this converts a given basepair
    location within the sequence to where it should appear on the actual screen'''
    multiplier = (screen_width - (2 * margin_width)) / seq_length
    return (bp_loc * multiplier) + margin_width

def determine_length_of_surface(filepath):
    '''Determines the needed length of the surface based on how many records there are'''
    with open(filepath, 'r') as file:
        num_lines = len(file.readlines())
    return ((num_lines / 2) * 350) + 175

def cycle_colors():
  '''Cycles through a list of colors indefinitely.'''
  yield from itertools.cycle(color_palette)

color_generator = cycle_colors()


# Create unique list of motifs, assign each a unique random color and put into dictionary
unique_motif_list = []

with open(args.motif, 'r') as motif_file:
    for i, line in enumerate(motif_file): 
        motif = UniqueMotif(line.strip().lower(), i, next(color_generator))
        unique_motif_list.append(motif)


# Not assuming the incoming file is oneline, so a temp one-line fasta file is created, we also get the number of lines to define the length of the image
bioinfo.oneline_fasta(args.fasta, "temp_oneline.fasta")

screen_height = int(determine_length_of_surface("temp_oneline.fasta"))

# Set surface, paint white
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, screen_width, screen_height)
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
        record.draw_graph()
        record.exon_object.draw_exon()
        record.gene_object.write_gene()

        for motif in record.motif_object_list:
            motif.draw_line(150)
        
        record_num += 1

    for unique_motif in unique_motif_list:
        unique_motif.write_key()
        unique_motif.draw_color_square()


os.remove('temp_oneline.fasta')

prefix = args.fasta.replace(".fasta", "")
prefix = prefix.replace(".fa", "")
surface.write_to_png(f'{prefix}.png')

    