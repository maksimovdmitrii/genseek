#    Copyright 2020 Dmitrii Maksimov

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

def detect_rotatble(connectivity_matrix):
    """Detection of all rotatable bonds
    1. The bonds does not contain terminate atoms
    2. 
    3. 
    """
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