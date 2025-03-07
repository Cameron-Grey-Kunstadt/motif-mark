# OOP Motif Mark

### Objective

This code generates a png image that allows for visualization of found motifs on fasta records.

### Inputs
-f, --fasta, Input FASTA file, must be unzipped, can be 1seq per line but does not have to be, introns must be lowercase and EXONS uppercase

-m, --motif, Input motif file, expecting one motif per line, case insensitive

### Outputs
Single png file, filename is based on the input fasta file

### Code
**motif-mark-oop.py** is the main file to run. This file operates in approximately the following steps:

- Takes in FASTA and motif file (also converts FASTA file to a oneline-per-seq version if not already)
- The FASTA file is used to determine the necessary height for the output image, and the pycairo surface and context are created
- The motif file is used to make a list of UniqueMotif objects, which store the motif, some information, and have functions to draw the motif on the final image
- The FASTA file is used to create a FastaRecord object for every record in the file, this object parses out the header info and saves it, as well as the seq
- Each FastaRecord object generates a motif_object_list upon initialization. This list is created by in-class functions parsing the FASTA sequence for each of the possible UniqueMotifs. Each found motif gets created as a MotifInstance object (a class that inherites UniqueMotif).
- Each FastaRecord object uses its draw function to draw the horizontal line, exon, label, and all the motifs (represented as vertical lines)
- Each UniqueMotif object uses a draw function to create the key at the bottom, along with the color squares
- The image is written out

**bioinfo.py** contains some useful functions, primarily converting the FASTA file to a oneline-per-seq version

**Example_Image.png** an example of how the graphs can look

**Funny_Example_Image.png** an example of how the graphs look if you like to have fun

#### Notes
- This code used to have a section where it could put a photo in the background of the graph, but as that required the use of additional packages, it was removed
- A color cycler is used to run through whatever color pallete you define at the beginning, feel free to make better use of the colors and fonts of this tool.

