import gRASPA as g
import numpy as np

def CombineFrameworkAdsorbate(FrameworkData, AdsData):
  # Since there are (usually) more framework atoms, just add TrialData to frameworkData
  Data = dict()
  LastFrameworkMolecule = int(FrameworkData['MolID'][-1]) + 1
  #TrialData['Trial_MolID'] += LastFrameworkMolecule
  Data['MolID']  = np.concatenate((FrameworkData['MolID'],  AdsData['MolID'] + LastFrameworkMolecule))
  Data['pos']    = np.concatenate((FrameworkData['pos'],    AdsData['pos']))
  Data['charge'] = np.concatenate((FrameworkData['charge'], AdsData['charge']))
  Data['Type']   = np.concatenate((FrameworkData['Type'],   AdsData['Type']))
  return Data

def calculate_angle(vector1, vector2):
  #Calculate the angle (in degrees) between two vectors using the dot product
  dot_product = np.dot(vector1, vector2)
  magnitude1 = np.linalg.norm(vector1)
  magnitude2 = np.linalg.norm(vector2)

  cos_theta = dot_product / (magnitude1 * magnitude2)
  angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Clip to handle numerical errors
  angle_deg = np.degrees(angle_rad)
  return angle_deg

def WriteDataToPOSCAR(Vars, systemId, DATA, filename = "DATA", exclusion = []):
  Box = g.GetBox(Vars, systemId)
  VectorA = Box["Cell"][0:3]
  VectorB = Box["Cell"][3:6]
  VectorC = Box["Cell"][6:9]
  gamma   = calculate_angle(VectorA, VectorB)
  beta    = calculate_angle(VectorA, VectorC)
  alpha   = calculate_angle(VectorB, VectorC)
  a = Box["Cell"][0]
  b = Box["Cell"][4]
  c = Box["Cell"][8]
  DATA_poscar = np.column_stack((DATA['pos']['x'], DATA['pos']['y'], DATA['pos']['z'], DATA['Type']))
  DATA_poscar = DATA_poscar[DATA_poscar[:, 3].argsort(kind='stable')]
  DATA_poscar_type, DATA_poscar_atomnumber = np.unique(DATA_poscar[:, 3], return_counts=True)
  DATA_poscar_xyz = DATA_poscar[:, :3]
  PseudoAtom_Symbols = Vars.PseudoAtoms.Symbol

  # 写入文件 'POSCAR'
  with open(f'{filename}.vasp', 'w') as file:
    file.write(f'{filename}.vasp\n')
    file.write('   1.00000000000000\n')
    file.write(f'    {VectorA[0]} {VectorA[1]} {VectorA[2]}\n')
    file.write(f'    {VectorB[0]} {VectorB[1]} {VectorB[2]}\n')
    file.write(f'    {VectorC[0]} {VectorC[1]} {VectorC[2]}\n')
    file.write(" ".join([f"{PseudoAtom_Symbols[int(t)]}" for t in DATA_poscar_type]) + "\n")
    file.write(f"{' '.join(map(str, DATA_poscar_atomnumber.astype(int)))}\n")
    file.write('Cartesian\n')
    for row in DATA_poscar_xyz:
      file.write('   '.join(map(str, row)) + '\n')
