SynOrF
=======================

SynOrF.py is a python script that determines shared orthologs between two bacteria genomes. The algorithm takes the neighborhood and the sequence homology into account and thereby finds orthologs that are conserved on a sequence and genomic location level. 

As an input the amino acid sequence (.faa file) sequence and a gff3 file of both genomes is needed.

The algorithm is described in: http://journal.frontiersin.org/article/10.3389/fmicb.2017.00218/full

#Requirements:
python (Was tested with version 2.7)

Blast+ (Was tested with version 2.2.28+, available ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

NetworkX (Was tested with version 1.9.1, available https://networkx.github.io/)

#Usage: 
python SynOrF.py strain1_protein.faa strain2_protein.faa strain1_annotation.gff strain2_annotation.gff strain1name.strain2name > orthologs.tab


