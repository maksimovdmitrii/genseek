from genseek.modules import *
from ase.io import read, write
from random import random, randint, uniform


class Structure:

    def __init__(self, parameters):
        self.atoms = read(parameters["geometry"][0], 
                        format=parameters["geometry"][1])
        self.connectivity_matrix_full = create_connectivity_matrix(self.atoms, bothways=True) 
        self.connectivity_matrix_isolated = create_connectivity_matrix(self.atoms, bothways=False)
        self.list_of_torsions = detect_rotatble(self.connectivity_matrix_isolated)

        if "mic" in parameters.keys():
            self.pbc = parameters["mic"]
            self.mic = True
            self.atoms.set_cell(self.pbc)
            self.atoms.set_pbc(True)
        self.molecules = [self.atoms.copy() for i in range(parameters["number_of_replicas"])]


    def create_configuration(self, parameters):

        def make_torsion(self, parameters):
            if parameters["configuration"]["torsions"]["values"] == "random":
                torsions = np.array([randint(-180, 180) 
                                    for i in self.list_of_torsions])
            return torsions


        def make_orientation(self, parameters):
            if parameters["configuration"]["orientations"]["values"] == "random":
                quaternion = produce_quaternion(
                    randint(-180, 180), 
                    np.array([random(),
                              random(),
                              random()]))
            else:
                angle = parameters["configuration"]["orientations"]["angle"] 
                x = parameters["configuration"]["orientations"]["x"] 
                y = parameters["configuration"]["orientations"]["y"] 
                z = parameters["configuration"]["orientations"]["z"] 
                quaternion = produce_quaternion(
                    randint(angle[0], angle[1]), 
                    np.array([uniform(x[0], x[1]),
                              uniform(y[0], y[1]),
                              uniform(z[0], z[1])]))             
            return quaternion

        def make_com(self, parameters):
            if parameters["configuration"]["coms"]["values"] == "restricted":
                x = parameters["configuration"]["coms"]["x_axis"] 
                y = parameters["configuration"]["coms"]["y_axis"] 
                z = parameters["configuration"]["coms"]["z_axis"] 
                com = np.array([uniform(x[0], x[1]), 
                                 uniform(y[0], y[1]), 
                                 uniform(z[0], z[1])])
            else:
                 com = np.array([uniform(0, 10), 
                                 uniform(0, 10), 
                                 uniform(0, 10)])
            return com

        torsions = make_torsion(self, parameters)
        quaternion = make_orientation(self, parameters)
        coms = make_com(self, parameters)
        configuration = np.hstack((torsions, quaternion, coms))  


        if len(self.molecules) > 1:
            for i in range(len(self.molecules) -1):
                if parameters["configuration"]["torsions"]["same"]:
                    pass
                else:
                    torsions = make_torsion(self, parameters)

                if parameters["configuration"]["orientations"]["same"]:
                    pass
                else:
                    quaternion = make_orientation(self, parameters)
                
                if parameters["configuration"]["coms"]["same"]:
                    pass
                else:
                    coms = make_com(self, parameters)   
                vec = np.hstack((torsions, quaternion, coms))
                configuration = np.hstack((configuration, vec))
   
        return configuration


    def apply_configuration(self, configuration):
    # molecules, configuration, list_of_torsions, connectivity_matrix_isolated):
        for i in range(len(self.molecules)):
            k=-1
            for torsion in self.list_of_torsions:
                k+=1
                z=i*(len(self.list_of_torsions)+4+3)+k
                fixed_indices = carried_atoms(
                                self.connectivity_matrix_isolated, torsion)

                self.molecules[i].set_dihedral(angle=configuration[z],
                                          a1=torsion[0],
                                          a2=torsion[1],
                                          a3=torsion[2],
                                          a4=torsion[3],
                                          indices=fixed_indices)
            # Set orientation
            quaternion_set(self.molecules[i], 
                            produce_quaternion(configuration[z+1], 
                                                np.array([configuration[z+2],
                                                        configuration[z+3],
                                                        configuration[z+4]])),
                                                0, len(self.molecules[i])-1)
            # Set center of mass
            set_centre_of_mass(self.molecules[i], np.array([configuration[z+5], 
                                                            configuration[z+6], 
                                                            configuration[z+7]]))




class Fixed_frame:

    def __init__(self, parameters):

        # Minimum Image Conventio
        self.mic = False
        if "fixed_frame" in parameters.keys():
            self.fixed_frame = read(parameters["fixed_frame"][0], 
                                format=parameters["fixed_frame"][1])

        if "mic" in parameters.keys():
            self.pbc = parameters["mic"]
            self.mic = True
            self.fixed_frame.set_cell(self.pbc)
            self.fixed_frame.set_pbc(True)

    def get_len(self):
        return len(self.fixed_frame)