'''
Created on Mar 16, 2013

@author: Prateek
'''
import networkx as nx
import mst
import sys
#===============================================================================
# INITIALIZATIONS
#===============================================================================
fil = sys.argv[1]
G = nx.read_gexf(fil)
nodes = G.nodes()
start = sorted(G.nodes())[0]
tour = [start]
source = start
cost = 0
S = []
fx = 0
hx = 0 
gx = 0
#===============================================================================
# LET THE SALESMAN TRAVEL!
#===============================================================================
while len(tour) != len(nodes):
    for i in G.neighbors(start):
        if i in tour:
            continue
        diff = set(tour) - set(source)
        setdiff = set(nodes)-diff
        #------- Subgraph fed to form MST for heuristic
        subG = G.subgraph(list(setdiff))
        hx = mst.heuristic(subG)
        gx = G[start][i]['weight']
        fx = gx+hx
        S.append((fx,i))
#===============================================================================
# S now contains all the (key,item) pairs thrown into heap
#===============================================================================
    H = mst.makeheap(S)
    S = []
    low = mst.deletemin(H)
    #----------------------------------- Extracting the node with the lowest key
    cur = low.item
    cost = cost+G[start][cur]['weight']
    tour.append(cur)
    start = cur
end = tour[-1]
cost = cost + G[source][end]['weight']
print "Tour: ",' '.join(tour)
print "Cost: ", cost