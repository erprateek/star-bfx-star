'''
Created on Mar 16, 2013

@author: Prateek
'''
#------------------------------------------------- Reference: mst.py from class
#---------- Addition: funtion heuristic to calculate the heuristic cost of a mst
import math
ARITY = 3  # the branching factor of the d-Heaps

#=======================================================================
# d-Heap
#=======================================================================

class HeapItem(object):
    """Represents an item in the heap"""
    def __init__(self, key, item):
        self.key = key
        self.item = item
        self.pos = None

def makeheap(S):
    """Create a heap from set S, which should be a list of pairs (key, item)."""
    heap = list(HeapItem(k,i) for k,i in S)
    for pos in xrange(len(heap)-1, -1, -1):
        siftdown(heap[pos], pos, heap)
    return heap

def findmin(heap):
    """Return element with smallest key, or None if heap is empty"""
    return heap[0] if len(heap) > 0 else None

def deletemin(heap):
    """Delete the smallest item"""
    if len(heap) == 0: return None
    i = heap[0]
    last = heap[-1]
    del heap[-1]
    if len(heap) > 0:
        siftdown(last, 0, heap)
    return i

def heapinsert(key, item, heap):
    """Insert an item into the heap"""
    heap.append(None)
    hi = HeapItem(key,item)
    siftup(hi, len(heap)-1, heap)
    return hi

def heap_decreasekey(hi, newkey, heap):
    """Decrease the key of hi to newkey"""
    hi.key = newkey
    siftup(hi, hi.pos, heap)

def siftup(hi, pos, heap):
    """Move hi up in heap until it's parent is smaller than hi.key"""
    p = parent(pos)
    while p is not None and heap[p].key > hi.key:
        heap[pos] = heap[p]
        heap[pos].pos = pos
        pos = p
        p = parent(p)
    heap[pos] = hi
    hi.pos = pos

def siftdown(hi, pos, heap):
    """Move hi down in heap until its smallest child is bigger than hi's key"""
    c = minchild(pos, heap)
    while c != None and heap[c].key < hi.key:
        heap[pos] = heap[c]
        heap[pos].pos = pos
        pos = c
        c = minchild(c, heap)
    heap[pos] = hi
    hi.pos = pos

def parent(pos):
    """Return the position of the parent of pos"""
    if pos == 0: return None
    return int(math.ceil(pos / ARITY) - 1)

def children(pos, heap):
    """Return a list of children of pos"""
    return xrange(ARITY * pos + 1, min(ARITY * (pos + 1) + 1, len(heap)))

def minchild(pos, heap):
    """Return the child of pos with the smallest key"""
    minpos = minkey = None
    for c in children(pos, heap):
        if minkey == None or heap[c].key < minkey:
            minkey, minpos = heap[c].key, c
    return minpos

#===============================================================================
# UNION FIND FOR KRUSKAL'S
#===============================================================================

class ArrayUnionFind:
    """Holds the three "arrays" for union find"""
    def __init__(self, S):
        self.group = dict((s,s) for s in S) # group[s] = id of its set
        self.size = dict((s,1) for s in S) # size[s] = size of set s
        self.items = dict((s,[s]) for s in S) # item[s] = list of items in set s
        
def make_union_find(S):
    """Create a union-find data structure"""
    return ArrayUnionFind(S)
    
def find(UF, s):
    """Return the id for the group containing s"""
    return UF.group[s]

def union(UF, a,b):
    """Union the two sets a and b"""
    assert a in UF.items and b in UF.items
    # make a be the smaller set
    if UF.size[a] > UF.size[b]:
        a,b = b,a
    # put the items in a into the larger set b
    for s in UF.items[a]:
        UF.group[s] = b
        UF.items[b].append(s)
    # the new size of b is increased by the size of a
    UF.size[b] += UF.size[a]
    # remove the set a (to save memory)
    del UF.size[a]
    del UF.items[a]


def kruskal_mst(G):
    """Return a minimum spanning tree using kruskal's algorithm"""
    # sort the list of edges in G by their length
    Edges = [(u, v, G[u][v]['weight']) for u,v in G.edges()]
    Edges.sort(cmp=lambda x,y: cmp(x[2],y[2]))

    UF = make_union_find(G.nodes())  # union-find data structure

    # for edges in increasing weight
    mst = [] # list of edges in the mst
    for u,v,d in Edges:
        setu = find(UF, u)
        setv = find(UF, v)
        # if u,v are in different components
        if setu != setv:
            mst.append((u,v))
            union(UF, setu, setv)
    return mst


def heuristic(G):
    mst = kruskal_mst(G)
    mst_cost = 0
    for (x,y) in mst:
        wt = G[x][y]['weight']
        mst_cost = mst_cost+wt
    return mst_cost