#!/usr/bin/env python

"""
This script determines the annotations that need to be updated and uses BLAT to generate the psl output file for each annotation. 
"""

from Bio import SeqIO
import subprocess
import os
import sys


def bash_cmd(command):
    """
    Calls a bash command and returns the standard output.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, cwd=os.curdir)
    proc_stdout = process.communicate()[0].strip()
    print("BASH_cmd_stdout: {}".format(proc_stdout))


def get_annotations_to_update(fasta1_filename,fasta2_filename):
    """
    Returns a dictionary with the differences between two multifasta files that have all the same header
    lines but different sequences. The dictionary has sequences from the first file (fasta1_filename), which should be from the outdated
    version of the genome.
    """

    # load in the files
    fasta1 = SeqIO.parse(fasta1_filename, "fasta")
    fasta2 = SeqIO.parse(fasta2_filename, "fasta")

    # transfer the fasta biopython generator handles to dictionaries (hash tables)
    dic_1,dic_2 = {},{}
    for annotation in fasta1:
        dic_1[str(annotation.id)] = str(annotation.seq)
    for annotation in fasta2:
        dic_2[str(annotation.id)] = str(annotation.seq)

    # get annotations that need to be updated
    annotations_to_update = {}
    for annotation in dic_1.keys():
        if annotation not in dic_2.keys():  # if the annotation is not present using the new genome (e.g. the chromosome was removed)
            annotations_to_update[annotation] = dic_1[annotation]  
        elif dic_1[annotation] != dic_2[annotation]:  # if the sequences do not match for that annotation
            annotations_to_update[annotation] = dic_1[annotation]

    print("\n{} annotations that need to be updated were found".format(len(annotations_to_update)))

    print("\nThe following annotations need to be updated")
    for annotation_key in annotations_to_update.keys():
        print(annotation_key)
        print("...length of annotation in old genome: {} bp".format(len(dic_1[annotation_key])))
        if annotation_key in dic_2:
            print("...length of annotation in new genome: {} bp\n".format(len(dic_2[annotation_key])))
        else:
            print("...the annotation was not in the new genome (the chromosome may have been removed, and may be found via BLAT)\n")
    return annotations_to_update


def generate_query_fas_and_BLAT_result_psls(annotations_to_update, genome_to_blat):
    for incorrect_annotation_key in annotations_to_update.keys():

        query_filename ="{}_query.fa".format(incorrect_annotation_key)
        print("\nWriting the query fasta file {}".format(query_filename))
        file_lines=['>'+incorrect_annotation_key+'\n',annotations_to_update[incorrect_annotation_key]]
        with open(query_filename, 'w+') as file_handle:
            file_handle.writelines(i for i in file_lines)

        blat_result_filename = "{}_blat_result.psl".format(incorrect_annotation_key)
        print("Running BLAT and writing psl file {}".format(blat_result_filename))
        bash_cmd("blat {} {} -t=dna -q=dna {}".format(genome_to_blat, query_filename, blat_result_filename))


if __name__ == "__main__":
    #inputs
    genome_filename_old = sys.argv[1]
    genome_filename_new = sys.argv[2]
    gtf_filename_old = sys.argv[3]

    # make fa files of the sequences from either genome
    bash_cmd("gffread {} -g {} -w sequences_oldgtf_oldgenome.fa".format(gtf_filename_old, genome_filename_old))
    bash_cmd("gffread {} -g {} -w sequences_oldgtf_newgenome.fa".format(gtf_filename_old, genome_filename_new))

    # get dictionary of annotations that need to be updated {annotation:sequence from old genome, ...}
    annotations_to_update = get_annotations_to_update("sequences_oldgtf_oldgenome.fa","sequences_oldgtf_newgenome.fa")
    
    # make querry fa and BLAT and write output to psl file
    generate_query_fas_and_BLAT_result_psls(annotations_to_update, genome_filename_new)

    print("blat.py Done.")