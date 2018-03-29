#!/usr/bin/env python
"""
This script takes in three inputs (an older version of a genome, a newer version of a genome,
and a gtf file made from the older genome). This script makes a new gtf file with the 
updated annotations using the newer version of the genome. This script uses Biopython and
the command line BLAT ('blat - Standalone BLAT v. 36x2 fast sequence search command line tool').
"""

from Bio import SeqIO,SearchIO
import Bio.SearchIO.BlatIO as BlatIO
import subprocess
import os
import csv

def bash_cmd(command):
    """
    Calls a bash command and returns the standard output.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, cwd=os.curdir)
    proc_stdout = process.communicate()[0].strip()
    print("BASH_cmd_stdout: {}".format(proc_stdout))


def get_annotations_to_update(fasta1_filename,fasta2_filename):
    """
    Returns a dictionary with the differeces between to multifasta files 
    that have all the same headerlines but maybe different sequences
    The dictionnary has sequences from the first file
    """
    print("\nDetermining which annotations to update")
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
        if i not in dic_2.keys():
            annotations_to_update[i] = dic_1[i]
        else:
            if dic_1[i] != dic_2[i]:
                annotations_to_update[i] = dic_1[i]

    print("...{} annotations to update were found".format(len(annotations_to_update)))
    print("...the following annotations need to be updated:")
    for i in annotations_to_update.keys():
        print(i)
    print("...the genomic sequences for these annotations were extracted from the 2006 genome\n")
    return annotations_to_update


# Step1: Determine which annotation coordinates need to be updated
old_gtf_filename = "Couvillion_sRNA_precursor_annotations_for_2006_genome.gtf"
print("\nMaking CDS files for each genome version using the gtf file")
print("...making CDS fa using 2006 genome")
bash_cmd("~/opt/bin/gffread {} -g TTA1_200603_genome_altered_headerlines.fa -w Couvillion_Ranges_CDS_oldGTF_FROM2006.fa".format(old_gtf_filename))
print("...making CDS fa using 2014 genome")
bash_cmd("~/opt/bin/gffread {} -g T_thermophila_June2014_assembly.fa -w Couvillion_Ranges_CDS_oldGTF_FROM2014.fa".format(old_gtf_filename))
annotations_to_update = get_annotations_to_update("Couvillion_Ranges_CDS_oldGTF_FROM2006.fa", "Couvillion_Ranges_CDS_oldGTF_FROM2014.fa")

print("Starting BLAT & gathering now coordinates")
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
print("...new coordinates gathered\n")


print("Editing GTF files\n")
# load in older gtf
gtf_lines=[]
with open(old_gtf_filename,'r') as openhandle:  # open old gtf file
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
new_gtf_filename = "all_Couvillion_Ranges_for_RNAi_FROM2014_made.gtf"
with open(new_gtf_filename, 'w+') as file:
    for line in gtf_lines:
        file.writelines("\t".join(line)+"\n")
print("...new GTF file with the updated coordinates was written to {} \n".format(new_gtf_filename))

# make fasta from gtf and check
print("To verify that BLAT functioned as anticipated a new fa file from the new annotations and the 2014 genome...")
print("will be made and compared to the fa from the 2006 genome and the original annotations.")
bash_cmd("~/opt/bin/gffread {} -g T_thermophila_June2014_assembly.fa -w Couvillion_Ranges_CDS_FROM2014_check_after_update.fa".format(new_gtf_filename))

dif_dic_after_update = get_annotations_to_update("Couvillion_Ranges_CDS_oldGTF_FROM2006.fa","Couvillion_Ranges_CDS_FROM2014_check_after_update.fa")

print("\nPrinting annotations that resulted in a different sequence")
for i in dif_dic_after_update.keys():
	print(i)
	print(dif_dic_after_update[i])
	print("\n")


print("\nMoving new files to /blat_querry_files folder")
bash_cmd("mkdir blat_querry_files")
bash_cmd("mv *_querry.fa *_blat_result.psl blat_querry_files")

print("\nScript Done.")   