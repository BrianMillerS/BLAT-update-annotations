#!/usr/bin/env python

"""
This script uses the blat output files (.psl files) to get the single top blat hit (according to the 'number of nucleotide matches as a percentage of querry length').
Then a new gtf file is made using the scaffold, start, and stop information from each of these top hits.
""" 

from Bio import SeqIO,SearchIO
import Bio.SearchIO.BlatIO as BlatIO
from operator import itemgetter
import os
import csv
import sys


def get_psl_file_paths():
    """
    gathers the full paths of all psl files in /BLAT_psl_files
    """
    psl_files = []
    psl_dir = os.getcwd() + "/BLAT_psl_files"  # get directory name with all the psl files
    cwdir_files = os.listdir(psl_dir)  # get all the files in that directory

    for file in cwdir_files:
        if file.endswith(".psl"):  # make sure all files are psl files
            psl_files.append(psl_dir + "/" + file)  # add to list

    return(psl_files)


def strand(int_input):
    if int_input == 1:
        return "+"
    elif int_input == -1:
        return "-"
    else:
        return "ISSUE, no strand was identified"


def get_coordinates_from_best_hit(psl_files):
    print("\nReporting info on best result from each psl file")

    updated_anntotations = {}  # initilize new dir for output

    n_best_hits_under_25 = 0  # initilize number of best hits (one per annotation) that have a 'number of matches as a percentage of querry length' less than 25%

    for psl_file in psl_files:

        annotation_name = psl_file.split("/")[-1][0:-16]  # get the name of the annoatation corresponding the the psl file

        # parse BLAT .psl file and add all hits to list
        blat_results = []
        psl_result = SearchIO.read(psl_file, 'blat-psl')
        for hit in psl_result:  # for every scaffold (hit)
            for HSP in hit:  # for every collection of fragments (HSP)
                # [scaffold, start position, stop position, percent identity (UCSC formula), number of matches, querry length, number of matches as a percentage of querry length, strand hit was from, strand querry was from]
                for subHSP in HSP:
                    hit_to_add = [hit.id, HSP.hit_start + 1, HSP.hit_end, HSP.ident_pct, HSP.ident_num, psl_result.seq_len, 100 * (float(HSP.ident_num) / psl_result.seq_len), strand(subHSP.hit_strand)]
                    blat_results.append(hit_to_add)

        # order the blat_results list with the greatest 'number of matches as a percentage of querry length' first
        blat_results_ordered = sorted(blat_results, key=itemgetter(6), reverse=True)

        # get the top most blat result
        best_hit = blat_results_ordered[0]  # get the best hit according to the 'number of matches as a percentage of querry length' (0 = best, 1 = second best, 2 = 3rd best, etc...)

        # print some results about the best hit
        print("\n{}".format(annotation_name))
        print("...Best hit: Scaffold: {}".format(best_hit[0]))
        print("...Best hit: Start: {}".format(best_hit[1]))
        print("...Best hit: End: {}".format(best_hit[2]))
        # print("...Best hit: Strand: {}".format(best_hit[7]))
        print("...Best hit: Percent identity (UCSC's formula): {}%".format(best_hit[3]))
        print("...Best hit: Number of matches: {} nts".format(best_hit[4]))
        print("...Query Length: {} nts".format(best_hit[5]))
        print("...Number of matches as a percent of query length: {}%".format(best_hit[6]))

        # add the scaffold, start, and end information to the dict with new annotations
        updated_anntotations[annotation_name] = [best_hit[0], best_hit[1], best_hit[2], best_hit[7]]  # [scaffold, start, end, strand] added to dict with new annotations

        # report other stats from psl file
        if best_hit[6] < 25:  # if 'number of matches as a percentage of querry length' is less than 25%
            n_best_hits_under_25 += 1  # add one to counter

    # potentially report issue about top best hits
    n_pslfiles = len(psl_files)
    if n_best_hits_under_25 > 0:
        print("WARNING: {} of the {} psl files have a best BLAT hit with a 'number of matches as a percentage of querry length' less than 25%".format(n_pslfiles, n_best_hits_under_25))

    print("\nNew coordinates gathered\n")
    return(updated_anntotations)


def make_new_gtf_with_updated_coordinants(old_gtf_filename, new_gtf_filename, updated_annotations):

    # load in older gtf
    gtf_lines = []
    with open(old_gtf_filename, 'r') as openhandle:  # open old gtf file
        reader = csv.reader(openhandle, delimiter='\t')
        for i in reader:
            gtf_lines.append(i)  # create a data structure that is [['line1_item0','line1_item1'...],['line2_item0','line2_item2'...],...]

    # alter gtf with new coordinates
    for new_anot_key in updated_annotations:  # for every annotation to update
        for line in gtf_lines:
            if new_anot_key == line[8].split('"')[1]:  # if annotation name in gtf = annotation name in hash table
                line[0] = updated_annotations[new_anot_key][0]  # change the scaffold
                line[3] = str(updated_annotations[new_anot_key][1])  # change the start position
                line[4] = str(updated_annotations[new_anot_key][2])  # change the end position

    # write changes to a new gtf
    with open(new_gtf_filename, 'w+') as file:
        for line in gtf_lines:
            file.writelines("\t".join(line) + "\n")
    print("The new gtf file with the updated coordinates was written to {}".format(new_gtf_filename))


if __name__ == "__main__":
    # inputs
    gtf_filename_old = sys.argv[1]
    gtf_filename_new_to_make = sys.argv[2]

    psl_file_paths = get_psl_file_paths()

    print("\nParse each psl output and gather coordinates from primary (best) blat hit")
    new_coordinates = get_coordinates_from_best_hit(psl_file_paths)

    print("Making a new .gtf annotations file with the new coordinates:\n")
    make_new_gtf_with_updated_coordinants(gtf_filename_old, gtf_filename_new_to_make, new_coordinates)

    print("parse_generate_new_annotations.py Done.")
