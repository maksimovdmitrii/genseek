from ase import neighborlist
from ase.build import molecule
from scipy import sparse
import numpy as np
import operator
import sys

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

def create_connectivity_matrix(atoms, bothways):
    cutOff = neighborlist.natural_cutoffs(atoms)
    neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=bothways)
    neighborList.update(atoms)
    connectivity_matrix = neighborList.get_connectivity_matrix()

    return connectivity_matrix

def detect_rotatble(connectivity_matrix):
    """Detection of all rotatable bonds
    2. The bonds does not contain terminate atoms
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


def carried_atoms(connectivity_matrix_isolated, positions):
    """ Returns list of carried atoms """
    graph = construct_graph(connectivity_matrix_isolated)
    graph_with_break = insertbreak(graph, positions[1], positions[2])
    if positions[2] in list(getroots(graph_with_break).values())[0]:
        return list(getroots(graph_with_break).values())[0]
    else:
        return list(getroots(graph_with_break).values())[1]



def unit_vector(vector):
    """ Returns the unit vector of the vector."""
    if np.linalg.norm(vector) == 0.:
        vector = np.array([0.,0.,0.0000000001])    #       Not to forget to check again this section
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Returns angle between two vectors"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.dot(v1_u, v2_u))*180./np.pi


def translate(point, coord):
    translated = coord[:] - point[:]
    return translated


def translate_back(point, coord):
    translated = coord[:] + point[:]
    return translated


def mult_quats(q_1, q_2):
    Q_q_2 = np.array([[q_2[0], q_2[1], q_2[2], q_2[3]],
                     [-q_2[1], q_2[0], -q_2[3], q_2[2]],
                     [-q_2[2], q_2[3], q_2[0], -q_2[1]],
                     [-q_2[3], -q_2[2], q_2[1], q_2[0]]])
    q_3 = np.dot(q_1, Q_q_2)
    return q_3


def unit_quaternion(q):
    ones = np.ones((1,4))
    ones[:,0] = np.cos(q[0]*np.pi/180/2)
    vec = np.array([q[1], q[2], q[3]])
    vec = unit_vector(vec)
    ones[:,1:] = vec*np.sin(q[0]*np.pi/180/2)
    quaternion = ones[0]
    return quaternion


def rotation_quat(coord, q):
    q = unit_quaternion(q)
    R_q = np.array([[1 - 2*q[2]**2 - 2*q[3]**2, 2*q[1]*q[2] - 2*q[0]*q[3], 2*q[1]*q[3] + 2*q[0]*q[2]],
                    [2*q[2]*q[1] + 2*q[0]*q[3], 1 - 2*q[3]**2 - 2*q[1]**2, 2*q[2]*q[3] - 2*q[0]*q[1]],
                    [2*q[3]*q[1] - 2*q[0]*q[2], 2*q[3]*q[2] + 2*q[0]*q[1], 1 - 2*q[1]**2 - 2*q[2]**2]])
    rotated = np.dot(R_q, coord.transpose())
    return rotated.transpose()


def Rotation(coord, point, quaternion):
    trans = translate(point, coord)
    rotate = rotation_quat(trans, quaternion)
    final = translate_back(point, rotate)
    return np.array(final)


def produce_quaternion(angle, vector):
    ones = np.ones((1, 4))
    ones[:, 0] = angle
    ones[:, 1:] = unit_vector(vector[:])
    quaternion = ones[0]
    return quaternion


def produce_coords_and_masses(coords, masses):
    zeros = np.zeros((len(coords), 4))
    zeros[:, :3] = coords[:]
    zeros[:, 3] = masses[:]
    return zeros

def measure_quaternion(atoms, atom_1_indx, atom_2_indx):

    # To revise and test
    coords = atoms.get_positions()
    orient_vec = unit_vector(coords[atom_2_indx] -
                                      coords[atom_1_indx])
    x_axis = np.array([1, 0, 0])
    z_axis = np.array([0, 0, 1])
    center = atoms.get_center_of_mass()          
    inertia_tensor = atoms.get_moments_of_inertia(vectors=True)    
    eigvals = inertia_tensor[0]
    eigvecs = inertia_tensor[1]     
    z_index = np.argmax(eigvals)                                      
    x_index = np.argmin(eigvals)                                      
    if np.dot(unit_vector(eigvecs[z_index]), orient_vec) < 0:
        eigvecs[z_index] = -eigvecs[z_index]
    ang_1 = angle_between(eigvecs[z_index], z_axis)                   
    vec_1 = np.cross(eigvecs[z_index], z_axis)                        
    quat_1 = produce_quaternion(ang_1, vec_1)                           
    rotated_1 = Rotation(coords, center, quat_1)       
    atoms.set_positions(rotated_1 )
    orient_vec_2 = unit_vector(rotated_1[atom_2_indx] - rotated_1[atom_1_indx])
    eigs_after = atoms.get_moments_of_inertia(vectors=True)[1]
    if np.dot(unit_vector(eigs_after[x_index]), orient_vec_2) < 0:
        eigs_after[x_index] = -eigs_after[x_index]
    angle_x = angle_between(eigs_after[x_index], x_axis)
    if np.dot(np.cross(unit_vector(eigs_after[x_index]), x_axis), z_axis) > 0:
        angle_x = -angle_x
    quaternion_of_the_molecule = np.array([angle_x, eigvecs[z_index, 0], eigvecs[z_index, 1], eigvecs[z_index, 2]])
    return quaternion_of_the_molecule

def align_to_axes(atoms, atom_1_indx, atom_2_indx):

    coords = atoms.get_positions() 
    center = atoms.get_center_of_mass()
    quaternion = measure_quaternion(atoms, atom_1_indx, atom_2_indx)
    vec = np.cross(quaternion[1:], np.array([0, 0, 1]))
    angle = angle_between(quaternion[1:], np.array([0, 0, 1]))
    quat_1 = produce_quaternion(angle, vec)
    rotation_1 = Rotation(coords, center, quat_1)
    angle_2 = -quaternion[0]
    quat_2 = produce_quaternion(angle_2, np.array([0, 0, 1]))
    rotation_2 = Rotation(rotation_1, center, quat_2)
    return atoms.set_positions(rotation_2)


def quaternion_set(atoms, quaternion, atom_1_indx, atom_2_indx):

    coords = atoms.get_positions()
    center = atoms.get_center_of_mass()
    align_to_axes(atoms, atom_1_indx, atom_2_indx)
    first_rot = produce_quaternion(quaternion[0], np.array([0, 0, 1]))
    rotation_1 = Rotation(atoms.get_positions(), center, first_rot)
    angle_2 = angle_between(np.array([0, 0, 1]), quaternion[1:])
    vec_2 = np.cross(np.array([0, 0, 1]), quaternion[1:])
    quat_2 = produce_quaternion(angle_2, vec_2)
    rotation_2 = Rotation(rotation_1, center, quat_2)
    return atoms.set_positions(rotation_2)


def internal_clashes(molecules, connectivity_matrix_full):
    
    clashes = True
    for i in range(len(molecules)):
        a = sorted(create_connectivity_matrix(molecules[i], bothways=True).keys(), key=lambda element: (element[1:]))
        b = sorted(connectivity_matrix_full.keys(), key=lambda element: (element[1:]))
        if operator.eq(set(a), set(b)):
            pass
        else:
            clashes = False

    return clashes


def intramolecular_clashes(molecules, mic):

    all_atoms = molecules[0].copy()
    for molecule in molecules[1:]:
        all_atoms.extend(molecule)
    distances = all_atoms.get_all_distances(mic=mic).reshape(len(all_atoms), len(all_atoms))

    for i in range(len(molecules)):
        values = np.ones(len(molecules[i])**2).reshape(len(molecules[i]), len(molecules[i])) * 100
        distances[len(molecules[i])*i:len(molecules[i])*i + len(molecules[i]) ,
                  len(molecules[i])*i:len(molecules[i])*i + len(molecules[i]) ] = values

    return not all(i >= 1.5 for i in distances.flatten()) 


def clashes_with_fixed_frame(molecules, fixed_frame, mic):

    mols = molecules[0].copy()
    for molecule in molecules[1:]:
        mols.extend(molecule)
    all_atoms = mols + fixed_frame
    distances = all_atoms.get_all_distances(mic=mic).reshape(len(all_atoms), len(all_atoms))
    values_mol = np.ones(len(mols)**2).reshape(len(mols), len(mols)) * 100
    distances[0:len(mols) ,
              0:len(mols) ] = values_mol
    values_fixed = np.ones(len(fixed_frame)**2).reshape(len(fixed_frame), len(fixed_frame)) * 100
    distances[len(mols) :len(mols) + len(fixed_frame)  ,
              len(mols) :len(mols) + len(fixed_frame) ] = values_fixed

    return not all(i >= 2.0 for i in distances.flatten())

def all_right(molecules, fixed_frame, connectivity_matrix_full, mic):

    ready = False
    if internal_clashes(molecules, connectivity_matrix_full):
        if not intramolecular_clashes(molecules, mic=mic):
            if not clashes_with_fixed_frame(molecules, fixed_frame, mic):
                ready = True
    return ready

def create_blacklist(molecules, list_of_torsions):

    if len(molecules) > 1:
        torsions = np.array([0 for i in list_of_torsions])
        quaternion = produce_quaternion(0, np.array([0, 0, 1]))
        value_com = np.array([0, 0, 0])
        blacklist_one = np.hstack((torsions, quaternion, value_com))
        blacklist = np.hstack((torsions, quaternion, value_com)) 
        for i in range(len(molecules) -1):
            blacklist = np.concatenate((blacklist, blacklist_one), axis=0)
    else:
        torsions = np.array([0 for i in list_of_torsions])
        quaternion = produce_quaternion(0, np.array([0, 0, 1]))
        value_com = np.array([0, 0, 0])
        blacklist = np.hstack((torsions, quaternion, value_com))        
    return blacklist

from random import random, randint
def create_internal_vector(molecules, list_of_torsions):

    if len(molecules) > 1:
        torsions = np.array([randint(0, 360) for i in list_of_torsions])
        quaternion = produce_quaternion(randint(0, 360), np.array([random() for i in range(3)]))
        value_com = np.array([randint(0, 13), randint(0, 13), randint(33, 38)])
        vector = np.hstack((torsions, quaternion, value_com))
        for i in range(len(molecules) -1):
            torsions = np.array([randint(0, 360) for i in list_of_torsions])
            quaternion = produce_quaternion(randint(0, 360), np.array([random() for i in range(3)]))
            value_com = np.array([randint(0, 13), randint(0, 13), randint(33, 38)])
            vec = np.hstack((torsions, quaternion, value_com))
            vector = np.hstack((vector, vec))
    else:
        torsions = np.array([randint(0, 360) for i in list_of_torsions])
        quaternion = produce_quaternion(randint(0, 360), np.array([random() for i in range(3)]))
        value_com = np.array([randint(0, 13), randint(0, 13), randint(33, 38)])
        vector = np.hstack((torsions, quaternion, value_com))        
    return vector

def add_to_blacklist(vector, blacklist):
    return np.vstack((blacklist, vector))

def not_in_blacklist(vector, blacklist):
    check = False
    dist = [np.linalg.norm(point - vector) for point in blacklist]
    if len(dist) > 1:
        if all(i > 100 for i in dist):
            check= True
    else:
        if dist[0] > 100:
            check= True
    return check

def apply_to_molecules(molecules, vector, list_of_torsions, connectivity_matrix_isolated):
    for i in range(len(molecules)):
        k=-1
        for torsion in list_of_torsions:
            k+=1
            z=i*(len(list_of_torsions)+4+3)+k
            fixed_indices = carried_atoms(connectivity_matrix_isolated, torsion)
            molecules[i].set_dihedral(angle=vector[z],
                                      a1=torsion[0],
                                      a2=torsion[1],
                                      a3=torsion[2],
                                      a4=torsion[3],
                                      indices=fixed_indices)
        # Set orientation
        quaternion_set(molecules[i], produce_quaternion(vector[z+1], 
                                                        np.array([vector[z+2],
                                                                 vector[z+3],
                                                                 vector[z+4]])),
                                                                0, len(molecules[i])-1)
        # Set center of mass
        set_centre_of_mass(molecules[i], np.array([vector[z+5], vector[z+6], vector[z+7]]))

def measure_molecules(molecules, list_of_torsions):
    vector = []
    for i in range(len(molecules)):
        for torsion in list_of_torsions:
            vector.append(molecules[i].get_dihedral(
                                      a1=torsion[0],
                                      a2=torsion[1],
                                      a3=torsion[2],
                                      a4=torsion[3]))
        # Set orientation
        vector.append(measure_quaternion(molecules[i], 0, len(molecules[i])-1))
        # Set center of mass
        vector.append(molecules[i].get_center_of_mass())
    return vector


def assign_random_state(molecule, list_of_torsions):
    for torsion in list_of_torsions:
        value = randint(0, 360)
        fixed_indices = carried_atoms(connectivity_matrix_isolated, torsion)
        molecules[i].set_dihedral(angle=value,
                                  a1=torsion[0],
                                  a2=torsion[1],
                                  a3=torsion[2],
                                  a4=torsion[3],
                                  indices=fixed_indices)
    quaternion_set(molecules[i], produce_quaternion(value, np.array([value, value, value])), 0, len(atoms)-1)
    value_com = np.array([randint(0, 13),
                  randint(0, 13),
                  randint(32, 36)])
    set_centre_of_mass(molecules[i], value_com)