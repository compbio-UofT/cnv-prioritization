#!/usr/bin/env python

import csv
import pprint
import cPickle as pickle
import networkx as nx
import sys
import argparse
import numpy as np
import time
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

pp = pprint.PrettyPrinter(indent=0)

parser = argparse.ArgumentParser(description='Parse gene mania data to a list of genes mapped to their neighbours sorted by dijkstra distance.')
parser.add_argument('mania_file', help="Gene mania file.")
parser.add_argument('out_file', help="Output.")
parser.add_argument('cutoff', type=float, help="Only consider interactions with a weight >= this cutoff.")
parser.add_argument('start', type=int, help="Line in the input file to start parsing.")
parser.add_argument('length', type=int, help="Number of lines in the input file to parse.")
args = parser.parse_args()

totalstart = time.clock()
start = time.clock()

def timer(s):
    global start
    global totalstart
    totalelapsed = (time.clock() - totalstart)
    elapsed = (time.clock() - start)
    start = time.clock()
    sys.stderr.write('%s...%s...%s\n' % (s, elapsed, totalelapsed))

timer('Loading gene mania file')
genenet=nx.Graph()
f = open(args.mania_file, 'r')
f.readline()
max_weight = 0
min_weight = np.inf
for line in f:
    (a, b, w) = line.strip().split()
    w = float(w)
    if w > max_weight:
        max_weight = w
    if w < min_weight:
        min_weight = w
    if w >= args.cutoff and 'ENSG' in a and 'ENSG' in b:
        genenet.add_edge(a, b, {'weightinv':w})

for a,b in genenet.edges():
    genenet[a][b]['weight'] = max_weight/genenet[a][b]['weightinv']

# HARDCODE
# cutoff to take the top n nodes by dijkstra
top_cutoff = 30

os.system("mkdir -p 1 2 3 4")

# HARDCODE
timer('Loop through nodes: %s' % 'Start')
f = open(args.out_file, 'w')
counter = 0
for nname in genenet.nodes()[args.start - 1 :args.start - 1 + args.length]:
    counter += 1
    
    timer('Loop through nodes: %s' % nname)

    length, path = nx.single_source_dijkstra(genenet, nname)

    import operator
    sorted_length = sorted(length.iteritems(), key=operator.itemgetter(1))

    print "%s_%s" % (counter, nname)
    plt.figure(1)
    plt.hist([b for a,b in sorted_length], 100)
    plt.savefig("1/%s_%s_1.png" % (counter, nname))
    plt.figure(2)
    plt.plot([b for a,b in sorted_length])
    plt.savefig("2/%s_%s_2.png" % (counter, nname))
    plt.figure(3)
    plt.plot([b for a,b in sorted_length[:30]])
    plt.savefig("3/%s_%s_3.png" % (counter, nname))
    plt.figure(4)
    plt.hist([b for a,b in sorted_length[:30]], 10)
    plt.savefig("4/%s_%s_4.png" % (counter, nname))

    from sets import Set
    neighbours_set = Set()
    for name, length in sorted_length[0:top_cutoff]:
        neighbours_set |= Set(path[name])

    genesub = genenet.subgraph(neighbours_set)

    # add neighbours that don't exist in genesub as sinks
    for n in genesub.nodes():
        if n == 'sink':
            continue

        out_neighbours = Set(genenet.neighbors(n)) - Set(genesub.neighbors(n))

        for on in out_neighbours:
            # add current weight of sink to new node if it exists
            w = genenet[n][on]['weight']
            winv = genenet[n][on]['weightinv']
            if 'sink' in genesub[n]:
                w += genesub[n]['sink']['weight']
                winv += genesub[n]['sink']['weightinv']
                # HARDCODE: 50
                # w = 50.0

            genesub.add_edge(n, 'sink', weight=w)
            genesub[n]['sink']['weightinv'] = winv

    import numpy as np

    sink_vector = np.zeros(len(genesub))
    if 'sink' in genesub.nodes():
        sink_index = genesub.nodes().index('sink')
        sink_vector[sink_index] = 1

    pi = nx.to_numpy_matrix(genesub).getA()

    # invert to actual weights
    # change inf -> 0
    invpi = max_weight/pi
    indpi = invpi == np.inf
    invpi[indpi] = 0

    # sum to a vector
    spi = np.sum(invpi, axis=1)
    # put the sum along the diagonal to create self weights at 1/2
    d = (np.diag(spi))*1
    # add the diagonal to the original transition matrix
    dpi = invpi + d
    # calc a new sum
    sdpi = np.sum(dpi, axis=1)
    # normalize over the sum
    ndpi = np.array(dpi / sdpi[:,np.newaxis])
    # set sink column to keep the values
    ndpi[0,:] = sink_vector

    # put full initial probability on the starting gene
    s = np.zeros(len(genesub))
    s[genesub.nodes().index(nname)] = 1

    # walk the graph
    ss = s
    # HARDCODE: number of steps
    for j in range(3):
        ss = ss.dot(ndpi)

    sum(ss)
    sss = ss.argsort()

    ret = []
    sorted_gene_by_length = map(lambda x: x[0], sorted_length[0:top_cutoff])
    sorted_length_dict = dict(sorted_length[0:top_cutoff])

    for i in range(1, len(sss)+1):

        #node name
        nn = genesub.nodes()[sss[-i]]

        if nn == 'sink':
            dijkstra_rank = top_cutoff+1+1
            genesub.node['sink']['dijkstra'] = 1.0
        else:
            if nn == nname:
                genesub.node[nn]['dijkstra'] = 1.0
            else:
                genesub.node[nn]['dijkstra'] = sorted_length_dict[nn]
            dijkstra_rank = sorted_gene_by_length.index(nn)+1

        # walk weight
        walkw = ss[sss[-i]]
        # print walkw
        tup = (nn, walkw, i, dijkstra_rank, genesub.node[nn]['dijkstra'])
        print tup

        if nn != nname and nn != 'sink':
            ret.append('%s:%s' % (nn, max_weight/genesub.node[nn]['dijkstra']))

        genesub.node[nn]['walk_rank'] = float(i)
        genesub.node[nn]['dijstra_rank'] = float(dijkstra_rank)

    print

    print >> f, '%s\t%s' % (nname, ','.join(ret))
    sys.stdout.flush()

    nx.write_graphml(genesub, '%s_%s.xml' % (counter, nname))

f.close()

timer('Loop through nodes: %s' % 'Done')

f = open('pm_%s_%s.txt' % (args.start, args.length), 'w') 
print >> f, len(genenet.nodes())
print >> f, genenet.nodes()
f.close()
