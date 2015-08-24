SynOrF
=======================

SynOrF.py is a python script that determines shared orthologs between two bacteria genomes. The algorithm takes the neighborhood and the sequence homology into account and thereby finds orthologs that are conserved on a sequence and genomic location level. 

As an input the amino acid sequence (.faa file) sequence and a gff3 file of both genomes is needed.

#usage: 
python SynOrF.py ref_protein1.faa protein2.faa ref_annotation1.gff annotation2.gff strain1name.strain2name


