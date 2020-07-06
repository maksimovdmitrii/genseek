from ase import neighborlist
from ase.build import molecule
from scipy import sparse
import numpy as np

def construct_graph(connectivity_matrix):
    """Constructu graph from connectivity matrix"""
    graph = {}
    for i in connectivity_matrix.keys():
        if i[0] not in graph:
            graph[i[0]] = [i[1]]
        else:
            graph[i[0]].append(i[1])
    for i in connectivity_matrix.keys():
        if i[1] not in graph:
            graph[i[1]] = [i[0]]
        else:
            graph[i[1]].append(i[0])
    return graph


# Cycles detection
# https://algocoding.wordpress.com/2015/04/02/detecting-cycles-in-an-undirected-graph-with-dfs-python/
def cycle_exists(G):                      # - G is an undirected graph.              
    marked = { u : False for u in G }     # - All nodes are initially unmarked.
    found_cycle = [False]                 # - Define found_cycle as a list so we can change
                                          # its value per reference, see:
                                          # http://stackoverflow.com/questions/11222440/python-variable-reference-assignment
 
    for u in G:                           # - Visit all nodes.
        if not marked[u]:
            dfs_visit(G, u, found_cycle, u, marked)     # - u is its own predecessor initially
        if found_cycle[0]:
            break
    return found_cycle[0]
 
def dfs_visit(G, u, found_cycle, pred_node, marked):
    if found_cycle[0]:                                # - Stop dfs if cycle is found.
        return
    marked[u] = True                                  # - Mark node.
    for v in G[u]:                                    # - Check neighbors, where G[u] is the adjacency list of u.
        if marked[v] and v != pred_node:              # - If neighbor is marked and not predecessor,
            found_cycle[0] = True                     # then a cycle exists.
            return
        if not marked[v]:                             # - Call dfs_visit recursively.
            dfs_visit(G, v, found_cycle, u, marked)

def create_torsion_list(bond, graph):
    
    return (
    [i for i in graph[bond[0]] if i!=bond[1]][0],
    bond[0],
    bond[1],
    [i for i in graph[bond[1]] if i!=bond[0]][0])


def set_centre_of_mass(atoms, new_com):
    
    old_positions = atoms.get_positions()
    old_com = atoms.get_center_of_mass() 
    atoms.set_positions(old_positions - old_com + new_com)

def create_connectivity_matrix(atoms):
    cutOff = neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=False)
    neighborList.update(atoms)
    connectivity_matrix = neighborList.get_connectivity_matrix()
    return connectivity_matrix

def detect_rotatble(atoms):
    """Detection of all rotatable bonds
    2. The bonds does not contain terminate atoms
    2. 
    3. 
    """
    connectivity_matrix = create_connectivity_matrix(atoms)
    graph = construct_graph(connectivity_matrix)
    indx_not_terminal = [i for i in graph if len(graph[i]) > 1]
    conn = [i for i in connectivity_matrix.keys() 
            if all(k in indx_not_terminal for k in i)]   
    list_of_torsions = []
    # If no cycles in the molecule
    if not cycle_exists(graph):
        for bond in conn:
            # Check for the index order
            list_of_torsions.append(create_torsion_list(bond, graph))
    return list_of_torsions



def getroots(aNeigh):
    #    source: https://stackoverflow.com/questions/10301000/python-connected-components
    def findroot(aNode,aRoot):
        while aNode != aRoot[aNode][0]:
            aNode = aRoot[aNode][0]
        return aNode, aRoot[aNode][1]
    myRoot = {}
    for myNode in aNeigh.keys():
        myRoot[myNode] = (myNode,0)
    for myI in aNeigh:
        for myJ in aNeigh[myI]:
            (myRoot_myI, myDepthMyI) = findroot(myI, myRoot)
            (myRoot_myJ, myDepthMyJ) = findroot(myJ, myRoot)
            if myRoot_myI != myRoot_myJ:
                myMin = myRoot_myI
                myMax = myRoot_myJ
                if myDepthMyI > myDepthMyJ:
                    myMin = myRoot_myJ
                    myMax = myRoot_myI
                myRoot[myMax] = (myMax,max(myRoot[myMin][1]+1,myRoot[myMax][1]))
                myRoot[myMin] = (myRoot[myMax][0],-1)
    mytoret = {}
    for myI in aNeigh:
        if myRoot[myI][0] == myI:
            mytoret[myI] = []
    for myI in aNeigh:
        mytoret[findroot(myI,myRoot)[0]].append(myI)
    return mytoret

def insertbreak(graph, atom1, atom2):
    graph[atom1].pop(graph[atom1].index(atom2))
    graph[atom2].pop(graph[atom2].index(atom1))
    return graph


def carried_atoms(atoms, positions):
    """ Returns list of carried atoms """
    connectivity_matrix = create_connectivity_matrix(atoms)
    graph = construct_graph(connectivity_matrix)
    graph_with_break = insertbreak(graph, positions[1], positions[2])
    if positions[2] in list(getroots(graph_with_break).values())[0]:
        return list(getroots(graph_with_break).values())[0]
    else:
        return list(getroots(graph_with_break).values())[1]



