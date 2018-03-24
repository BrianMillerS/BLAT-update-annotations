#!/usr/bin/env python

from Bio import SeqIO,SearchIO
import Bio.SearchIO.BlatIO as BlatIO
import subprocess
import os
import csv

def bash_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, cwd=os.curdir)
    proc_stdout = process.communicate()[0].strip()
    print("BASH_cmd_stdout: " + proc_stdout + '\n')


def get_annotations_to_update(fasta1_filename,fasta2_filename):
    """
    Returns a dictionary with the differeces between to multifasta files 
    that have all the same headerlines but maybe different sequences
    The dictionnary has sequences from the first file
    """
    print("...determining which annotations to update")
    fasta1 = SeqIO.parse(fasta1_filename, "fasta")
    fasta2 = SeqIO.parse(fasta2_filename, "fasta")

    # transfer the fasta biopython generator handles to dictionaries (hash tables)
    dic_1,dic_2 = {},{}
    for i in fasta1:
        dic_1[str(i.id)] = str(i.seq)
    for i in fasta2:
        dic_2[str(i.id)] = str(i.seq)
    
    # get annotations that need to be updated
    annotations_to_update = {}
    for i in dic_1.keys():
        if dic_1[i] == dic_2[i]:
            continue
        else:
            annotations_to_update[i] = dic_1[i]
    print("...{} annotations to update were found".format(len(annotations_to_update)))
    return annotations_to_update


annotations_to_update = get_annotations_to_update("Couvillion_Ranges_CDS_FROM2006.fa", "Couvillion_Ranges_CDS_FROM2014.fa",)

print("Starting BLAT & gathering now coordinates...")
for incorrect_annotation_key in annotations_to_update.keys():
    blat_result_filename = "{}_blat_result.psl".format(incorrect_annotation_key)

    # write the BLAT querry to a file (must be done)
    file_lines=['>'+incorrect_annotation_key+'\n',annotations_to_update[incorrect_annotation_key]]
    querry_filename ="{}_querry.fa".format(incorrect_annotation_key)
    with open(querry_filename, 'w+') as file_handle:
        file_handle.writelines(i for i in file_lines)

    # run BLAT (writes output to file)
    bash_cmd("~/opt/bin/blat T_thermophila_June2014_assembly.fa {} -t=dna -q=dna {}".format(querry_filename,blat_result_filename))  # BLAT 
    
    #parse BLAT file
    updated_annotations = {}  # new dictionary for annotation name and the new coordinates
    qresult = SearchIO.read(blat_result_filename, 'blat-psl')
    for hit in qresult:
        if len(hit) == 1:  # if alignment has 1 HSP(alignment with similar scores/details)
            if len(hit[0]) == 1:  # if HSP has 1 fragment
                if hit[0].match_num == len(annotations_to_update[incorrect_annotation_key]):  # if the number of matched bases equals the total querry length
                    HSP = hit[0]
                    scaffold = hit.id
                    start = HSP.hit_start
                    end = HSP.hit_end
    updated_annotations[incorrect_annotation_key] = [scaffold,start,end]
    break
print("...new coordinates gathered")


print("Editing GTF files...")
# load in older gtf
gtf_lines=[]
with open("all_Couvillion_Ranges_for_RNAi_FROM2006.gtf",'r') as openhandle:  # open old gtf file
    reader = csv.reader(openhandle,delimiter='\t')
    for i in reader:
    	gtf_lines.append(i)  # create a data structure that is [['line1_item0','line1_item1'...],['line2_item0','line2_item2'...]]
# alter gtf with new coordinates
for new_anot_key in updated_annotations:  # for every annotation to update
    for line in gtf_lines:
        if new_anot_key == line[8].split('"')[1]:  # if annotation name in gtf = annotation name in hash table
            line[0]=updated_annotations[new_anot_key][0]  # change the scaffold
            line[3]=str(updated_annotations[new_anot_key][1])  # change the start position
            line[4]=str(updated_annotations[new_anot_key][2])  # change the end position

# write changes to a new gtf
filename = 'all_Couvillion_Ranges_for_RNAi_FROM2014_made.gtf'
with open(filename, 'w+') as file:
    for line in gtf_lines:
        file.writelines("\t".join(line)+"\n")


# make fasta from gtf and check

print("Script Done.")
