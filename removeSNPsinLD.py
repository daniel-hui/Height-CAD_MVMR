# Takes a list of pairs of SNPs in LD
# and a file with p-value of each SNP
# and gives a list of SNPs to remove such that none of the pairs of
# SNPs in LD remains.

# Usage: python removeRelatedIndividuals.py <LD_file> <file with P-values> <rsid column in p-file> <p-value column in p-file> > out.remove

# Garrett Wong
# 19 Dec 2014

# updated by Daniel Hui
# 11/13/17
# changed to:
# work with python 2.7+
# compare missingness for node removal with >2 but equal degree

#updated by Daniel Hui
#October 11, 2019
#changed to p-values instead of missingness, and SNPs instead of individuals
#for Voight lab at Penn
#note this script is originally from the Patsopoulos lab affl. Brigham and Women's, Harvard Medical School, Broad Insitute

#october 23, 2019
#changed rsid column and p-value column to input args

import sys
import networkx as nx

ldfile = sys.argv[1] #get_final_snp_list/LD_almostFinal/all_LD_sorted_greater.05.txt
pfile = sys.argv[2] #get_final_snp_list/Final_SNPlist_v1.txt
rs_col = int(sys.argv[3]) - 1 #eg. 1
p_col = int(sys.argv[4]) - 1 #eg. 9
			      #+1 since they would be columns from `cut` originally

# We use the following two functions to choose which SNPs to remove.
# Modelling the SNPs and the LD relation as nodes and edges,
# ideally we would find a maximum independent set. Unfortunately, this problem
# is NP-complete, and so we find a maximal independent set instead using a 
# greedy approach: we begin with the set of all nodes and remove the node
# with highest degree until the set is independent. For ties in degree,
# the node with the lower p-value is kept.
def disconnect(graph):
	'''adds a list of nodes which, when removed, leaves a maximal
	independent set to global list toRemove. Destroys the graph.'''
	# It's much easier to deal with connected subgraphs, so we simply
	# use this function as a shell to call disconnectCore, which deals
	# with a single connected subgraph.
	for component in nx.connected_component_subgraphs(graph):
		disconnectCore(component)

def disconnectCore(graph):
	'''Given a connected graph, adds nodes which must be removed to disconnect
	its nodes to toRemove. Recursive...'''
	if len(graph) == 1:
		return
	elif len(graph) == 2:
		# With a subgraph of two nodes, one must be removed. We always remove
		# the node with higher p-value.
		left, right = graph.nodes()
		if pD[left] < pD[right]:
			toRemove.append(right)
		else:
			toRemove.append(left)
	else:
		# If the graph has three or more nodes, we remove the node with highest
		# degree and recurse on each connected subgraph.
		degreeD = graph.degree()


		#daniel_edit:
		#finds node with max degree
		#this block did not work with python 2.7+
		#also modified to compare missingness
		max_deg = 0
		max_miss = 1
		nodeToRemove = None
		for node in degreeD:
			#python 2.6
			if(type(degreeD) is dict):
				if degreeD[node] > max_deg:
					max_deg = degreeD[node]
					nodeToRemove = node
					max_miss = pD[node]
				elif (degreeD[node] == max_deg and pD[node] > max_miss):
					nodeToRemove = node
					max_miss = pD[node] 
			#python 2.7+ and python 3
			else:
				if node[1] > max_deg:
					max_deg = node[1]
					nodeToRemove = node[0]
					max_miss = pD[node[0]]
				elif (node[1] == max_deg and pD[node[0]] > max_miss):
					nodeToRemove = node[0]
					max_miss = pD[node[0]] 
		#end daniel_edit


		graph.remove_node(nodeToRemove)
		toRemove.append(nodeToRemove)

		for component in nx.connected_component_subgraphs(graph):
			disconnectCore(component)


# We make the dictionary of SNP p-values:
pD = {} # keyed by rsid; value a p-value float
with open(pfile) as fp:
	fp.readline()
	for line in fp:
		line = line.rstrip().split()
		rs = line[rs_col]
		p = float(line[p_col])
		pD[rs] = p

# Then, we build a graph. Nodes are rs, and edges are LD


#and avoid checking for header each line
G = nx.Graph()
with open(ldfile, "r") as fp:
	fp.readline()
	for line in fp:
		fields = line.rstrip().split()
		rs1, rs2 = fields[2], fields[5]
		G.add_edge(rs1, rs2)


# Finally, we find the nodes to remove and print them: 
toRemove = []
disconnect(G)
for node in toRemove:
	print(node)

