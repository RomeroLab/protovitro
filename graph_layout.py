from sys import argv
import numpy
import sequence_tools
import matplotlib.pyplot as plt
import networkx as nx


filename = 'BM3heme_uniref50.aln'


# get the query sequence (the first seq)
with open(filename) as infile:
    for line in infile:
        if line[0]!='#' and len(line)>10:
            name,queryseq = line.split()
            break


# the positions where the query doesn't have a gap
querypos = [i for i in range(len(queryseq)) if queryseq[i]!='-']
queryseq = ''.join([queryseq[i] for i in querypos])



coverage_cut = 0.75 # want sequences to be at least 75%
identity_cut = 0.20

sequences = []
names = []
k = 0
with open(filename) as infile:
    for line in infile:
        if len(line)>10 and line[0]!='#':
            name,seq = line.split()
            seq = ''.join([seq[i] for i in querypos])
            cov = 1 - seq.count('-')/len(seq)
            iden = len([p for p in range(len(queryseq)) if queryseq[p]==seq[p]])/(cov*len(queryseq))
            if cov>coverage_cut and iden>identity_cut:
                sequences.append(seq)
                names.append(name)
                k += 1
                #if k==1000: break



align = [p for p in zip(*sequences) if set(p)!=set(['-'])]
sequences = [''.join(s) for s in zip(*align)]


G = nx.Graph()
#G.add_nodes_from(range(len(sequences)))


similarities = numpy.zeros((len(sequences),len(sequences)))
for i,s1 in enumerate(sequences):
    for j,s2 in enumerate(sequences):
        if s1<s2:
            print(i)
            HD = sequence_tools.hamming_dist(s1,s2)
            similarities[(i,j)] = HD
            similarities[(j,i)] = HD
            if HD<200:
                G.add_edge(i,j)



# write to GraphML file for loading into Cytoscape
nx.write_graphml(G,'cytoscape_graph50.xml')


# plot graph using Graphviz neato
pos = nx.nx_pydot.graphviz_layout(G,prog='neato')
nx.draw(G,pos,node_size=10)
plt.show()
