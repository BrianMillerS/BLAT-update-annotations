#!/bin/bash

# Author: Brian Miller
#
##DESCRIPTION##
# This script takes in three inputs (an older version of a genome, a newer version of a genome,
# and a gtf file made from the older genome). This script makes a new gtf file with the 
# updated annotations using the newer version of the genome. 
#
##USAGE##
# ./master_BLAT_and_generate_updated_annotations.sh [outdated version of genome] [new version of genome] [annotation file derived from outdated genome] [new annotation file to be made]
#
# For example:
# ./master_BLAT_and_generate_updated_annotations.sh old_genome.fa new_genome.fa outdated_annotations.gtf updated_annotations.gtf
#
# please note that an .fai file will be created for both versions of the genome and thus the genomes should be placed in the copy of the cloned repository
#
##DEPENDENCIES: programs## (add the programs to your $PATH so that they can be run as seen below)
# samtools faidx
# blat (Standalone BLAT v.36x2)
# gffread (from Cufflinks, http://cole-trapnell-lab.github.io/cufflinks/)
#
# NOTE that the scripts might throw an error when trying to run the executables, this can be prevented by replacing the above executables with their full paths in the scripts.
#
##OUTPUTS##
# This script will output a gtf file with the updated compatible with the new genome using the name provided in the [new annotation file to be made] field.
# Along with this there will be three files in the /terminal_outputs that describe the steps taken and any annotations that should be confirmed by hand.
# blat_terminal_output.txt has terminal outputs from when blat was run on each annotation sequences
# parse_generate_new_annotations_terminal_output.txt has terminal outputs from when the best blat hit coordinates were extracted from the blat output file
# compare_and_verify_terminal_output.txt has terminal outputs about any annotation whos sequence changed after the update
#
##DEPENDENCIES: python packages##
# Biopython (v1.70)


# make index files for the genomes
samtools faidx $1
samtools faidx $2

# run BLAT, write outputs to /BLAT_psl_files
python blat.py $1 $2 $3 |& tee blat_terminal_output.txt

mkdir query_fasta_files
mkdir BLAT_psl_files
mv *_blat_result.psl BLAT_psl_files
mv *_query.fa query_fasta_files

# parse the output files and extract the new coordinates, make a new gtf file
python parse_generate_new_annotations.py $3 $4 |& tee parse_generate_new_annotations_terminal_output.txt

# compare the old and new gtf file, and compare the sequences (sequences derived from 2006 gtf from and 2006 genome vs sequences derived from 2014 gtf and 2014 genome)
python compare_and_verify.py $2 $4 $3 |& tee compare_and_verify_terminal_output.txt

mkdir terminal_outputs
mv *terminal_output.txt terminal_outputs

mkdir sequence_processing_files
mv sequences_* sequence_processing_files

echo "Done."