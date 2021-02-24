from xml.dom.minidom import parse
import xml.dom.minidom
import numpy as np


AtomMass = {
    "H" : 1.007,
    "C" : 12.010,
    "N" : 14.006,
    "F" : 18.998,
    "Pb" : 207.2,
    "I" : 126.904
}


class OperateXsd(object):
    def __init__(self, file_name):
        DOMTree = xml.dom.minidom.parse(file_name)
        self.XSD = DOMTree.documentElement
        if self.XSD.hasAttribute("WrittenBy"):
            print("\nProgram:", self.XSD.getAttribute("WrittenBy"), '\n')

    def norm_vector(self, vector):
        return np.sqrt(float(vector[0])**2 + float(vector[1])**2 + float(vector[2])**2)

    def get_atom_infos(self):
        identity_mapping = self.XSD.getElementsByTagName("IdentityMapping")[0]
        
        atoms = identity_mapping.getElementsByTagName("Atom3d")
        positions = np.empty([len(atoms),3])
        elements = []
        for i,atom in enumerate(atoms):
            positions[i] = np.array(atom.getAttribute("XYZ").split(','))
            elements.append(atom.getAttribute("Components"))
        
        a_vector = identity_mapping.getElementsByTagName("SpaceGroup")[0].getAttribute("AVector").split(',')
        b_vector = identity_mapping.getElementsByTagName("SpaceGroup")[0].getAttribute("BVector").split(',')
        c_vector = identity_mapping.getElementsByTagName("SpaceGroup")[0].getAttribute("CVector").split(',')
        
        return positions, elements, [self.norm_vector(a_vector), self.norm_vector(b_vector), self.norm_vector(c_vector)]

    def get_geometric_center(self):
        positions, elements, vectors = self.get_atom_infos()
        geometric_center_frac = np.sum(positions, axis=0)/len(elements)
        geometric_center_cart = geometric_center_frac * np.array(vectors)
        return geometric_center_cart

    def get_barycenter(self):
        positions, elements, vectors = self.get_atom_infos()
        mass_sum1 = np.zeros(3)
        mass_sum2 = np.zeros(3)
        for i, elem in enumerate(elements):
            mass_sum1 += positions[i] * AtomMass[elem]
            mass_sum2 += AtomMass[elem]
        barycenter_frac = mass_sum1/mass_sum2
        barycenter_cart = barycenter_frac * np.array(vectors)
        return barycenter_cart
            

if __name__ == '__main__':
    app = OperateXsd("???.xsd")
    print(app.get_geometric_center())
    print(app.get_barycenter())
    
