import numpy as np
import numpy as np
from scipy.spatial.transform import Rotation as R

def angle(vector_1, vector_2): # find angle between vectors, radian
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    return np.arccos(dot_product)

def rotate_system(a: np.array, b: np.array, coord: np.array) -> np.array:
    """
    Parameters
    ----------
    a : rotatable vector
    b : direction vector
    coord : coordinations of our system
    Returns
    -------
    System will be rotate at the same angle as we rotate vector a to make it parallel with vector b
    """
    if np.linalg.norm(np.cross(a, b)) == 0: # calc mult vector, in case of parallel vector we choose any normal vector
        if a[0] == 0:
            rv = [np.pi, 0, 0]
        else:
            rv = np.array([-(a[1]+a[2])/a[0], 1, 1])/np.linalg.norm([-(a[1]+a[2])/a[0], 1, 1]) * np.pi
    else:
        rv = np.cross(a, b)/np.linalg.norm(np.cross(a, b)) * angle(a, b)
    
    r = R.from_rotvec(rv)
    return r.apply(coord)

class MolGro:
    def __init__(self, system=[], coords=[], numFirstAt=1):
        self.system = system
        self.coords = coords
        self.numFirstAt = numFirstAt
        self.index = -1
        
        
    def read(self, file):
        with open(file) as in_file:
            lenght = len(in_file.readlines())
        with open(file) as file:
            for cnt, line in enumerate(file):
                if 1 < cnt < lenght - 1:
                    line = [ int(line[:5]), line[5:10].strip(), line[10:15].strip(), 
                            int(line[15:20]), float(line[20:28]), 
                            float(line[28:36]), float(line[36:44]) ]
                    self.system.append(line[:4])
                    self.coords.append(line[4:])
        self.coords = np.array(self.coords)    
    
    def GC2zero(self):
        self.GC = np.mean(self.coords, axis=0)
        self.coords -= self.GC
    
    def GCback(self):
        self.coords += self.GC
        self.GC = np.array([0, 0, 0])
    
    def rotateSystem(self, a, b):
        self.coords = rotate_system(a, b, self.coords)
        
    def copy(self):
        return MolGro(self.system, self.coords, self.numFirstAt)
    
    def __str__(self): # modify to reses!
        system = [ self.system[i][:3] + [self.numFirstAt + i] + list(self.coords[i]) for i in range(len(self.system))]
        return '\n'.join([ '{:5d}{:<5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(*i) for i in system ]) + '\n'

    def __add__(self, other):
        return MolGro(self.system + other.system, np.vstack((self.coords, other.coords)))
    def __len__(self):
        return len(self.system)
    
    def __getitem__(self, index): # make it better!
        return self.system[index], self.coords[index]
    def __iter__(self):
        return self
    def __next__(self):
        self.index += 1
        return self.system[self.index], self.coords[self.index]
        
        
class ResesGro(MolGro):
    def __init__(self, system=[], coords=[], numFirstAt=1, systemFragments=[]):
        MolGro.__init__(self, system, coords, numFirstAt)
        self.index = -1
        self.systemFragments = systemFragments
        
    def split(self):
        systemReses=[[0, ]]
        for i in range(1, len(self.system)):
            if f'{self.system[i-1][0]}{self.system[i-1][1]}' != f'{self.system[i][0]}{self.system[i][1]}':
                systemReses[-1].append(i)
                systemReses.append([i+1])
        systemReses[-1].append(i+1)
        self.systemFragments = [MolGro(self.system[i[0]: i[1]], self.coords[i[0]: i[1]], self.system[i[0]][3]) for i in systemReses]
    
    def __add__(self, other):
        return ResesGro(systemFragments=self.systemFragments + other.systemFragments)
    def __len__(self):
        return len(self.systemFragments)
    def __getiterm__(self, index):
        return self.systemFragments[index]
    def __iter__(self):
        return self
    def __next__(self):
        self.index += 1
        return self.systemFragments[self.index]

class FragmentsGro(ResesGro):
    def split(self, start=frozenset(('CHT0C4', )), end=frozenset(('CHTNHO1', ))):
        systemFrag=[]
        for i in range(len(self.system)):
            key = f'{self.system[i][1]}{self.system[i][2]}' 
            if  key in start:
                systemFrag.append([i])
            if key in end:
                systemFrag[-1].append(i+1)
        self.systemFragments = [MolGro(self.system[i[0]: i[1]], self.coords[i[0]: i[1]], self.system[i[0]][3]) for i in systemFrag]
              
    
                               
                            
            
        
#    8CHT     C4 1173   3.030   2.242   4.365        
testSys = FragmentsGro()
testSys.read('box.gro')
testSys.split()
# testSys += testSys
