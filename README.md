Synteny_ortholog_finder
=======================

find_othologe_annotation.py is a python script that is supposed compare two Bacterial genomes against each other. By taking the neighborhood and the sequence homology the script find pairs of genes that are the same in the genomic context and in the sequence.

As an input the faa sequence and a gff3 file of both genomes is needed.

#usage: 
python find_othologe_annotation.py ref_protein1.faa protein2.faa ref_annotation1.gff annotation2.gff strain1name.strain2name


