#!/usr/bin/env python

"""
This script compares the sequences derived from the old and the new gtf annotation files, as a validation step. If the nucleotide 
sequences differ,then information about that annotaion is returned.
"""

from Bio import SeqIO,SearchIO
import subprocess
import os
import csv
import sys

def bash_cmd(command):
    """
    Calls a bash command and returns the standard output.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, cwd=os.curdir)
    proc_stdout = process.communicate()[0].strip()
    print("BASH_cmd_stdout: {}".format(proc_stdout))


def get_gtf_info(gtf_filename):
    # open a gtf file, transfer lines to list
    gtf_lines = []
    with open(gtf_filename,'r') as openhandle:  # open old gtf file
        reader = csv.reader(openhandle,delimiter='\t')
        for i in reader:
            gtf_lines.append(i) # create a data structure that is [['line1_item0','line1_item1'...],['line2_item0','line2_item2'...],...]
    
    # transfer scaffold, start, stop information to dictionary
    gtf_dict = {}
    for line in gtf_lines:
        name = line[8].split('"')[1]  # get name of annotation
        scaffold = line[0] # get the scaffold
        start = line[3]  # get the start position
        end = line[4]  # get the end position
        gtf_dict[name]=[scaffold,start,end]

    return gtf_dict
    

def return_fasta_dict(input_fa):
    # transfer the fasta biopython generator handles to dictionaries (hash tables)

    # load fa with biopython
    fasta = SeqIO.parse(input_fa, "fasta")

    fa_dict= {}
    for i in fasta:
        fa_dict[str(i.id)] = str(i.seq)
    return fa_dict


def compare_sequence_files(sequences_old_filename, sequences_new_filename, old_gtf_file, new_gtf_file):
    print("\nChecking for differences between the sequences produced from:")
    print("the old gtf version with the old genome, versus")
    print("the new gtf version (BLAT made file) with the new genome.\n")

    dict_old = return_fasta_dict(sequences_old_filename)
    dict_new = return_fasta_dict(sequences_new_filename)

    old_gtf_dict = get_gtf_info(old_gtf_file)
    new_gtf_dict = get_gtf_info(new_gtf_file)


    annotations_still_not_the_same = []
    for anot in dict_old.keys():
        if dict_old[anot] != dict_new[anot]:
            annotations_still_not_the_same.append(anot)

    print("{} annotations still do not have the same sequences and should be looked at by hand for verification\n".format(len(annotations_still_not_the_same)))

    for annotation in sorted(annotations_still_not_the_same):
        print(annotation)
        print("...length of annotation from old genome: {} bp".format(len(dict_old[annotation])))
        print("...length of annotation from new genome: {} bp".format(len(dict_new[annotation])))
        print("...coordinates from old genome: {} {}-{}".format(old_gtf_dict[annotation][0], old_gtf_dict[annotation][1], old_gtf_dict[annotation][2]))
        print("...coordinates from new genome: {} {}-{}\n".format(new_gtf_dict[annotation][0], new_gtf_dict[annotation][1], new_gtf_dict[annotation][2]))


if __name__ == "__main__":
    #inputs
    genome_filename_new = sys.argv[1]
    old_sequences = "sequences_oldgtf_oldgenome.fa"
    new_sequences = "sequences_newgtf_newgenome.fa"
    old_gtf = sys.argv[3]
    new_gtf = sys.argv[2]
    

    # extract sequences from the new genome using the updated gtf
    bash_cmd("gffread {} -g {} -w sequences_newgtf_newgenome.fa".format(new_gtf, genome_filename_new))
    
    compare_sequence_files(old_sequences, new_sequences, old_gtf, new_gtf)

    print("compare_and_verify.py Done.")