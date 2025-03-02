#!/usr/bin/env python

# Author: Cameron Grey Kunstadt, cameron.kunstadt@gmail.com

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Calculates average phred score for a string of character phred scores'''
    sum_score = 0
    for letter in phred_score:
        score = convert_phred(letter)
        sum_score += score
    return sum_score / len(phred_score)

def reverse_compliment(seq):
        '''This function generates the reverse compliment of given sequence,
        will strip newline characters'''

        seq = seq.strip()
        assert validate_base_seq(seq)

        compliments = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'G': 'C',
                       'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'g': 'c',
                       'N':'N', 'n':'n'}
        rev_seq = ""

        for base in seq:
            rev_seq = compliments[base] + rev_seq

        return rev_seq


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content(seq: str):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq)
    gc = 0
    seq = seq.upper()
    for letter in seq:
        if letter == 'G' or letter == 'C':
            gc += 1

    return (gc / len(seq))
    

def calc_median(lst: list) -> float:
    '''Calculates and returns the median of a 1 dimensional list that is already sorted'''

    if len(lst) % 2 == 0: 
        median1 = lst[len(lst)//2] 
        median2 = lst[len(lst)//2 - 1] 
        median = (median1 + median2)/2
    else: 
        median = lst[len(lst)//2]
        
    return median

def oneline_fasta(input_file, output_file):
    '''This function takes in an input FASTA file and output file name, and will read the input file, and if if has
       sequences on multiple lines, it will join them together, this new 'one-line sequence' FASTA file
       is written out to the output file'''
    with open(input_file,"r") as fh_in, open(output_file,"w") as fh_out:
        cur_seq = ""
        while True:
            line = fh_in.readline()
            if not line:
                fh_out.write(cur_seq)
                break
            elif line[0] == '>':
                if len(cur_seq) > 0:
                    fh_out.write(cur_seq + '\n')
                fh_out.write(line)
                cur_seq = ""
            else:
                cur_seq = cur_seq + line.strip('\n')


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    #print("Your convert_phred function is working! Nice job")

    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    #print("Passed DNA and RNA tests")
