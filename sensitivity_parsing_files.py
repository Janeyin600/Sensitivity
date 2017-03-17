import os
import sys

def parsing_mol2(mol2):
  """
  Return the atom type and partial charge of each atom
  """
  # Open the mol2 file of the receptor/ligand
  with open(mol2) as fr:
    lines = fr.read().splitlines()
    lines = list(line for line in lines if line) # Remove blank lines

  for i, line in enumerate(lines):
    if '@<TRIPOS>ATOM' == line.split()[0]:
      startline = i
    if '@<TRIPOS>BOND' == line.split()[0]:
      endline = i

  atom_type_list =  [line.split()[5] for i, line in enumerate(lines) if i < endline and i > startline]
  partial_charge_list = [line.split()[-1] for i, line in enumerate(lines) if i < endline and i > startline]
  
  zipped_list = zip(atom_type_list, partial_charge_list) 

  # Zip two lists together for returning
  return zipped_list
   
def parsing_gaff(atom_types, gaff_file):
  """
  Read the vdw parameters from the gaff data file
  """
  # Generate a small temporary file for faster lookup
  dict_param = {}
  with open(gaff_file) as f:
    lines = f.read().splitlines()
    lines = list(line for line in lines if line) # Remove blank lines

    for line in lines:
      splitline = line.split()
      if len(splitline) != 3:
         print('Aborted.There should be three columns in %s:'%(gaff_file))
         print('  atom type, radius value and epsilon value.\n')
         sys.exit(1)
      if line.split()[0] in atom_types:
         dict_param[splitline[0]] = [float(splitline[1]), float(splitline[2])] 

  return dict_param  

def parsing_water_ion_dat(water_param, ion):
  """
  Read the water/ion parameters from a user defined file
  """   
  ifind = 'no'
  with open(water_param) as f_in:
    lines = (line.strip(' \t\n\r') for line in f_in)
    lines = list(line for line in lines if line)  # Non-blank lines in a list
    if len(lines)!= 3:
      print('Aborted.There should be three columns in %s:'%(water_param))
      print('  atom type, radius value, epsilon value and partial charge.\n')
      sys.exit(1)

  for i in range(len(lines)):
      splitline = lines[i].split()
      if splitline[0].upper() == 'OW':
        ow_param = [float(splitline[1]), float(splitline[2])]
      elif splitline[0].upper() == 'HW':
        hw_param = [float(splitline[1]), float(splitline[2])]
      elif splitline[0] == ion:   
        ion_param = [float(splitline[1]), float(splitline[2])]
        ifind = 'yes'
     
  if ifind == 'no':
    print('Aborted! The ion name you provided in the sensitivity input file is not consistent with the parameter file %s.'%(ion))
    sys.exit(1)
  
  water_ion_params = []
  water_ion_params.append(ow_param)
  water_ion_params.append(hw_param)
  water_ion_params.append(hw_param)
  water_ion_params.append(ion_param)
  
  return water_ion_params  

def analysing_bonded_atoms(mol2, index_shift):
  """
  Return a list that contains all the bonded and one atoms aways (angle) atom pairs
  """
  # Open the mol2 file of the receptor/ligand
  # For the repceptor molecule the shift of the indexes will be 2
  # For the ligand molecule the shift of the indexes will be the number of host atoms plus 2 

  with open(mol2) as fr:
    lines = fr.read().splitlines()
    lines = list(line for line in lines if line) # Remove blank lines

  for i, line in enumerate(lines):
    if '@<TRIPOS>BOND' == line.split()[0]:
      startline = i
    if '@<TRIPOS>SUBSTRUCTURE' == line.split()[0]:
      endline = i
  
  # Find all the bonded atoms   
  bond_atom1 =  [int(line.split()[1])+index_shift for i, line in enumerate(lines) if i < endline and i > startline]
  bond_atom2 = [int(line.split()[2])+index_shift for i, line in enumerate(lines) if i < endline and i > startline]

  bond_atom_pairs = [list(i) for i in zip(bond_atom1, bond_atom2)]
  return bond_atom_pairs
   
