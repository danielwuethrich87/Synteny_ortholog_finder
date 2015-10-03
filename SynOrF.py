#bin/python

import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import math
import random, string
import operator


#usage: ref_protein1.faa protein2.faa ref_annotation1.gff annotation2.gff strain1.strain2

inputOptions = sys.argv[1:]

def main():


	run_blast(inputOptions)
	proteins=read_blast(inputOptions)
	proteins=construct_protein_list(inputOptions[3], proteins)
	
	proteins_query=read_proteins(inputOptions[0])
	proteins_query=construct_protein_list(inputOptions[2], proteins_query)	

	genome_graph=construct_reference_graph(proteins)
	genome_graph_query=construct_reference_graph(proteins_query)


	selected_protein_pairs=select_proteins(proteins, proteins_query, genome_graph, genome_graph_query)


	select_clusters_for_length(selected_protein_pairs, genome_graph, genome_graph_query)
#	check_neighborhood(valid_connection, proteins, genome_graph, genome_graph_query) # planed control step

	
def check_neighborhood(valid_connection, proteins, genome_graph, genome_graph_query):
	proteins_ordered=sorted(proteins.values(), key=lambda Protein: Protein.start)

	ordering_valid_connection=list()
	
	for protein in proteins_ordered:
		found_valid_neigbor=bool(0)
		if protein.name in valid_connection.keys():
			ref_neihgbors=genome_graph.neighbors(protein.name)

			for neihgbor in ref_neihgbors:
				if protein.name in valid_connection.keys() and neihgbor in valid_connection.keys():
					if len(nx.shortest_path(genome_graph_query,valid_connection[protein.name],valid_connection[neihgbor]))==2:	
						found_valid_neigbor=bool(1)

		if found_valid_neigbor==bool(1):

			ordering_valid_connection.append(protein.name+"\t"+valid_connection[protein.name])

	for i in ordering_valid_connection:
		print i


def select_clusters_for_length(selected_protein_pairs, genome_graph, genome_graph_query):
	# defines the shorteste length of a synteny otholog cluster
	pairs={}
	for n in selected_protein_pairs:	
		pairs[n[0]]=n[1]


	valid_graph=nx.Graph()
	edges=list()
	for edge in genome_graph.edges():


		if edge[0] in pairs.keys() and edge[1] in pairs.keys():
			valid_graph.add_edge(edge[0],edge[1])

	for cluster in nx.connected_components(valid_graph):
		if len(cluster)>4:
			query_nodes={}
			valid_query_graph=nx.Graph()
			for node in cluster:
				
				query_nodes[pairs[node]]=node


			for query_edge in genome_graph_query.edges():
				if query_edge[0] in query_nodes.keys() and query_edge[1] in query_nodes.keys():
					
					valid_query_graph.add_edge(query_edge[0],query_edge[1])
			
			for query_cluster in nx.connected_components(valid_query_graph):
				if len(query_cluster)>4:
					for query_node in query_cluster:
						print query_nodes[query_node]+"\t"+query_node
			

	


def select_proteins(proteins, proteins_query, genome_graph, genome_graph_query):
	# this method is double sorting, it only accepts othologues that are also neigbors in the query genome: and it also filters for the best connected orologue in the query genomes

	selected_proteins={}
	orthologues_of_proteins={}
	orthologues={}
	for query in proteins_query:
		orthologues[query]=list()	
	
	for protein in proteins:
		orthologues_of_proteins[protein]={}
		for neighbor in genome_graph.neighbors(protein):
			for blast_hit in proteins[protein].valid_blast_hits:
				for blast_hit_neighbor in proteins[neighbor].valid_blast_hits:
					if nx.has_path(genome_graph_query,blast_hit.name,blast_hit_neighbor.name):
						if len(nx.shortest_path(genome_graph_query,blast_hit.name,blast_hit_neighbor.name))==2:	
							if blast_hit.name in orthologues_of_proteins[protein]:
								orthologues_of_proteins[protein][blast_hit.name]+=1
							else:
								orthologues_of_proteins[protein][blast_hit.name]=1


		sorted_orthologues=sorted(orthologues_of_proteins[protein].iteritems(), key=operator.itemgetter(1),reverse=True)

		if len(sorted_orthologues)>1:

			if sorted_orthologues[0][1] != sorted_orthologues[1][1]:

				selected_proteins[sorted_orthologues[0][0]]=protein
				orthologues[sorted_orthologues[0][0]].append(sorted_orthologues[0])



		elif (len(sorted_orthologues)>0):
			
			selected_proteins[sorted_orthologues[0][0]]=protein		
			orthologues[sorted_orthologues[0][0]].append(sorted_orthologues[0])


	query_database_edges=list()
	for othologue in orthologues.keys():
		sorted_orthologues=sorted(orthologues[othologue], key=lambda tup: tup[1],reverse=True)
		if len(sorted_orthologues)>1:
			if sorted_orthologues[0][1] != sorted_orthologues[1][1]:
				edge=selected_proteins[sorted_orthologues[0][0]],sorted_orthologues[0][0]
				query_database_edges.append(edge)
		elif (len(sorted_orthologues)>0 ):
			edge=selected_proteins[sorted_orthologues[0][0]],sorted_orthologues[0][0]
			query_database_edges.append(edge)

	
	return query_database_edges

def construct_reference_graph(proteins):

	graph=nx.Graph()

	
	proteins_ordered=sorted(proteins.values(), key=lambda Protein: Protein.start)	

	all_contigs=set([proteins[n].contig for n in proteins.keys()])
	for contig in all_contigs:
		proteins_ordered_contig=list()
		for protein in proteins_ordered:
			if protein.contig==contig:
				proteins_ordered_contig.append(protein)

		for i in range(0,len(proteins_ordered_contig)):
			k=i+1
			z=i+2
			l=i+3
			
			if i == len(proteins_ordered_contig)-3:
				k=0
				z=1
				l=2
			elif i == len(proteins_ordered_contig)-2:
				k=0
				z=0
				l=1	
			elif i == len(proteins_ordered_contig)-1:
				k=0
				z=0			
				l=0

			graph.add_edge(proteins_ordered_contig[i].name,proteins_ordered_contig[k].name, weight=1)
			graph.add_edge(proteins_ordered_contig[i].name,proteins_ordered_contig[z].name, weight=0.5)
			graph.add_edge(proteins_ordered_contig[i].name,proteins_ordered_contig[l].name, weight=0.25)
	

	return graph

	

def construct_protein_list(in_file, proteins):
	# connects protein objects to genomic locations
	gff_annotations = [n for n in open(in_file,'r').read().replace("\r","").split("\n") if len(n)>0]
	genome_length=0
	for line in gff_annotations:
		if len(line.split("\t"))==9:
			if line.split("\t")[2]=="CDS":

				protein = proteins[line.split("ID=")[1].split(";")[0]]
				protein.start=int(line.split("\t")[3])
				protein.end=int(line.split("\t")[4])
				protein.genome_length=genome_length
				protein.contig=line.split("\t")[0]
				if line.split("\t")[6]=="-":
					protein.is_foreward=bool(0)

		elif line.find("##sequence-region")!=-1:
			if line.split(" ")[2]=="1":
				genome_length=int(line.split(" ")[3])




	return proteins

def read_blast(inputOptions):
	# parses the blast ouput for homologs
	proteins=read_proteins(inputOptions[1])
	proteins_of_blast=read_proteins(inputOptions[0])
	blast_output = [n for n in open("result/"+inputOptions[4]+".blast.annot.tab",'r').read().replace("\r","").split("\n") if len(n)>0]
	for line in blast_output:
		blast_hit = Blast_hit(line.split("\t")[1])
		start=int(line.split("\t")[8])
		end=int(line.split("\t")[9])
		if start > end:
			blast_hit.start=end-1
			blast_hit.end=start
			blast_hit.is_foreward=bool(0)
			
		else:
			blast_hit.start=start-1
			blast_hit.end=end
			blast_hit.is_foreward=bool(1)
		
		blast_hit.protein=proteins_of_blast[blast_hit.name]
		blast_hit.alignment_length=int(line.split("\t")[3])
	
		blast_hit.similarity=float(line.split("\t")[2])
		proteins[line.split("\t")[0]].blast_hits.append(blast_hit)
	
	proteins=filter_blast(proteins)	

	return proteins
	
def filter_blast(proteins):
	for key in proteins.keys():
		
		protein= proteins[key]

		for blast_hit in protein.blast_hits:
 
			if blast_hit.similarity>50 and blast_hit.alignment_length>= 0.5*float(len(protein.sequence)) and blast_hit.alignment_length>= 0.5*float(len(blast_hit.protein.sequence)):
				
				protein.valid_blast_hits.append(blast_hit)
				

	return proteins			
	
def run_blast(inputOptions):
	# runs the external blast search
	os.system("mkdir -p result")
	os.system("makeblastdb -in "+inputOptions[0]+" -dbtype prot -out result/reference" + " > database_info")
	os.system("blastp -query "+inputOptions[1]+" -db result/reference -outfmt 6 > result/"+inputOptions[4]+".blast.annot.tab")

def read_proteins(fasta_file):
	# constructs protein objects out of fasta sequences
	protein_fasta = [n for n in open(fasta_file,'r').read().replace("\r","").split("\n") if len(n)>0]
	proteins={}
	for line in protein_fasta:
		if line[0:1]==">":
			key=line[1:].split(" ")[0]
			protein=Protein(key)
			proteins[key]=protein

		else:
			proteins[key].sequence+=line.replace("*","")
			
	return proteins		


def graph_overlap(graph1, graph2):
	# creates the intersection of two graphs
	temp_graph=nx.Graph()


	edges=list()
	for edge in graph1.edges():
		try:


			if len(nx.shortest_path(graph2,edge[0],edge[1]))==2:		
				edges.append(edge)

		except :
    			pass
	temp_graph.add_edges_from(edges)
			
	return temp_graph


class Protein:
	def __init__(self, name):
		self.name=name
		self.sequence=""
		self.blast_hits=list()
		self.valid_blast_hits=list()
		self.start=0
		self.end=0
		self.is_foreward=bool(1)
		self.genome_length=0
		self.edges=list()
		self.contig=""

class Blast_hit:
	def __init__(self, name):
		self.name=name
		self.similarity=0
		self.start=0
		self.end=0
		self.is_foreward=bool(1)
		self.protein_length=0
		self.edge_name=""
		self.alignment_length_ratio=0
		self.best_alignment=bool(1)
		self.protein=0

class Edge:
	def __init__(self, name):
		self.name=name
		self.distance=0
		self.blast_hits=0
		self.valid=bool(1)


main()

