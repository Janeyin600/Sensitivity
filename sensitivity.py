#!/usr/bin/env python2
import os as os
import sys as sys
import glob as glob
import shutil as shutil
import datetime as dt
import signal as signal
import subprocess as sp
import re as re
import math
from def_functions import *
import sensitivity_parsing_files as sensitivity_parsing_files
import sensitivity_calc_coord as sensitivity_calc_coord

class Sensitivity:
    def __init__(self):
        """
        Default user and system parameters specification. These are now overwritten with the values found in sensitivity.in.
        More information on these variables can be found in the input file and the sensitivity paper.
        """
        self.output = 'data.dat'     # default Filename for the output
        self.rmol2 = 'None'          # mol2 file of the receptor 
        self.lmol2 = 'None'          # mol2 file of the ligand
        self.fcoord = 'None'         # coordinate file of the complex
        self.fparam = 'None'         # general force field parameter file   
        self.wparam = 'None'         # Water parameter file      

    def check_arguments(self):
        """
        The location/name of the input file is required; output file is optional  
        :return:
        """
        if len(sys.argv) == 1:
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 2) and ('-h' in sys.argv[1].lower()):
            help_message()
            # Return without an error
            sys.exit(0)
        elif (len(sys.argv) == 15) and ('-r' in sys.argv) and ('-l' in sys.argv) and ('-c' in sys.argv)\
          and ('-p' in sys.argv) and ('-w' in sys.argv):
            print('Running sensitivity analysis ...\n')
        elif (len(sys.argv) == 13) and ('-r' in sys.argv) and ('-l' in sys.argv) and ('-c' in sys.argv)\
          and ('-p' in sys.argv) and ('-w' in sys.argv):                          
            print('Running sensitivity analysis ...') 
        elif (len(sys.argv) < 13):
            help_message()
            print('Too few arguments. You need to provide an input file with the -i flag, the coordinate file of the complex (-c),')
            print ('the parameter file of the general force field (-p) and water/ions (-w), and the mol2 files with the -r and -l flags.\n')
            sys.exit(1)
        else:
            help_message()
            print('Not the right command. You need to provide an input file with the -i flag, the coordinate file of the complex (-c),')
            print ('the parameter file of the general force field (-p) and water/ions (-w), and the mol2 files with the -r and -l flags.\n')
            sys.exit(1)

        # Loop over the list of command line arguments and look for the flags
        for i, argv in enumerate(sys.argv[:-1]):
          if argv.lower() == '-i':           
            self.input_file = sys.argv[i+1]
            if not os.path.isfile(self.input_file):
              print('The input file %s does not exist.' % self.input_file)
              sys.exit(1)
          elif argv.lower() == '-r':
            self.rmol2 = sys.argv[i+1]
            if not os.path.isfile(self.rmol2):
              print('The receptor mol2 file %s does not exist.' % self.rmol2)
              sys.exit(1)
          elif argv.lower() == '-l':
            self.lmol2 = sys.argv[i+1]
            if not os.path.isfile(self.lmol2):
              print('The receptor mol2 file %s does not exist.' % self.lmol2)
              sys.exit(1)
          elif argv.lower() == '-c':
            self.fcoord = sys.argv[i+1]
            if not os.path.isfile(self.fcoord):
              print('The coordinate file for the complex: %s does not exist.' % self.fcoord)
              sys.exit(1)
          elif argv.lower() == '-p':
            self.fparam = sys.argv[i+1]
            if not os.path.isfile(self.fparam):
              print('The parameter file of the general force field: %s does not exist.' % self.fparam)
              sys.exit(1)
          elif argv.lower() == '-w':
            self.wparam = sys.argv[i+1]
            if not os.path.isfile(self.wparam):
              print('The parameter file of water/ions: %s does not exist.' % self.wparam)
              sys.exit(1)
          elif argv.lower() == '-o':
            self.output = sys.argv[i+1]
          
        if self.rmol2 == 'None':
            # if the mol2 file of the receptor is not provided
            help_message()
            print('Please provide the mol2 file of the receptor!')
            sys.exit(1)
        elif self.lmol2 == 'None':
            # if the mol2 file of the ligand is not provided
            help_message()
            print('Please provide the mol2 file of the ligand!')
            sys.exit(1)
        elif self.fcoord == 'None':
            # if the coordinate file of the complex is not provided
            help_message()
            print('Please provide the PDB file of the complex!')
            sys.exit(1)
        elif self.wparam == 'None':
            # if the parameter file of water/ions is not provided
            help_message()
            print('Please provide the PDB file of the complex!')
            sys.exit(1)
        elif self.fparam == 'None':
            # if the parameter file of the general force field is not provided
            help_message()
            print('Please provide the parameter file for the water/ions!')
            sys.exit(1)

    def process_input_file(self):
        """
        Process the sensitivity input file
        """
        with open(self.input_file) as f_in:
            # Delete all spaces/tabs at both ends       
            lines = (line.strip(' \t\n\r') for line in f_in)
            lines = list(line for line in lines if line)  # Non-blank lines in a list

        for i in range(0, len(lines)):
            # Combine the lines that belong to the same entry
            if not lines[i][0] == ';':
                lines[i] = lines[i].split(';')[0].split('=')
                if len(lines[i]) == 1:
                    j = i
                    while True:
                        if lines[j - 1] != ';':
                            lines[j - 1][1] += lines[i][0]
                            lines[i] = ';'
                            break
                        j -= 1
        
        for i in range(0, len(lines)):
            if not lines[i][0] == ';':
                lines[i][0] = lines[i][0].strip().lower()
                lines[i][1] = lines[i][1].strip()
                if lines[i][0] == 'num_water':
                    self.num_water = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'water_sites':
                    self.water_sites = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'ion_name':
                    self.ion_name = ismyinstance('string', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'num_ions':
                    self.num_ions = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'vdw_cutoff':
       	       	    self.vdw_cutoff = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'ee_cutoff':
       	       	    self.ee_cutoff = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                elif lines[i][0] == 'frames':
       	       	    self.frames = ismyinstance('int', lines[i][1], self.input_file, lines[i][0])
                else:
                    print('Wrong entry name in %s! ' % self.input_file)
                    print(lines[i][0] + '\n')
                    print('Please use the same keywords as in the template input file. Aborted.')
                    sys.exit(1)

        """
        print self.num_water
        print self.water_sites
        print self.ion_name
        print self.num_ions
        print self.vdw_cutoff
        print self.ee_cutoff
        print self.frames
        """

    def get_parameters(self):
      """
      Write a big 3D list with all the parameters: [[radius-for-atom1, epsilon-for-atom1, partial-charge-for-atom1],
      [radius-for-atom2, epsilon-for-atom2, partial-charge-for-atom2],...] 
      """
      # Get a list that contains all the atom types in the receptor and the ligand molecules
      receptor_atom_types = [ elem[0] for elem in sensitivity_parsing_files.parsing_mol2(self.rmol2)]
      ligand_atom_types = [ elem[0] for elem in sensitivity_parsing_files.parsing_mol2(self.lmol2)]
      self.num_receptor_atoms = len(receptor_atom_types)
      self.num_ligand_atoms = len(ligand_atom_types)      

      # Get the partial charge list
      charge_list = [ elem[1] for elem in sensitivity_parsing_files.parsing_mol2(self.rmol2)] + \
                    [ elem[1] for elem in sensitivity_parsing_files.parsing_mol2(self.lmol2)] 
      
      # Join two lists
      merged_list = receptor_atom_types + ligand_atom_types
      self.solute_atoms = len(merged_list) + + self.num_ions + 3 # take three dummy atoms into account
      self.total_atoms = self.solute_atoms + self.num_water * self.water_sites
      
      # Remove all the duplicates
      self.all_types = list(set(merged_list))

      vdw_param = sensitivity_parsing_files.parsing_gaff(self.all_types, self.fparam)
      vdw_list = []
      self.all_types.append(self.ion_name)

      # Parameters for the dummy atom hand-coded
      dummy_params = [[0,0],[0,0],[0,0]]

      for i, atom_type in enumerate(merged_list):
        if atom_type in vdw_param.keys():
           vdw_list.append(vdw_param[atom_type])

      # Parameters for the ion parameters
      water_ion_params = sensitivity_parsing_files.parsing_water_ion_dat(self.wparam, self.ion_name)
      
      self.allparam = dummy_params + vdw_list + self.num_ions * [water_ion_params[-1]] + self.num_water * water_ion_params[0:3] 
      self.tip3p_param =  [water_ion_params[0][0],water_ion_params[0][1]]
      self.atom_types =  ['Pb'] * 3 + receptor_atom_types + ligand_atom_types + [self.ion_name] * self.num_ions

    def get_connectivity(self):
      """
      Generate atom lists that contain the connectivity information
      """
      # Build a list that contains all the atom pairs within the solute (receptor and ligand)
      solute_atom_pairs = []

      for i in range(3, self.solute_atoms):
         for j in range(i+1, self.solute_atoms):
             solute_atom_pairs.append([i,j])
      
      receptor_bond_pairs = sensitivity_parsing_files.analysing_bonded_atoms(self.rmol2, 2)
      ligand_bond_pairs = sensitivity_parsing_files.analysing_bonded_atoms(self.lmol2, 2+self.num_receptor_atoms)
      self.bond_pairs = receptor_bond_pairs + ligand_bond_pairs
      
      # Find all the atom pairs that are one atom away (angles)
      angle_pairs = []
      angle_list = []     

      for i in range(len(self.bond_pairs)):
         for j in range(i+1,len(self.bond_pairs)):
           tmp_pair = list(set(self.bond_pairs[i]).symmetric_difference(set(self.bond_pairs[j])))
           if len(tmp_pair) == 2:
              tmp_pair.sort()
              angle_list.append(tmp_pair)
              connecting_atom = list(set(self.bond_pairs[i])&set(self.bond_pairs[j]))[0]
              tmp_pair.insert(1,connecting_atom)

      self.angle_pairs = [[angle[0], angle[2]] for angle in angle_list]
      angle_list =  sorted(angle_list)

      # Find all the atom pairs that are two atoms away (1-4)
      self.pair_list = []   # 1-4 pair list
      for i in range(len(angle_list)):
         for j in range(i+1, len(angle_list)):
             tmp_pair = find_1_4(angle_list[i], angle_list[j])
             if tmp_pair!=[]:
               self.pair_list.append(tmp_pair)

      self.pair_list = remove_duplicates_from_list(self.pair_list)

      # Get a nonbonded list (excluding bond, angle, and 1-4)
      for elem in self.angle_pairs:
           solute_atom_pairs.remove(elem)

      for elem in self.bond_pairs:
           solute_atom_pairs.remove(elem)

      for elem in self.pair_list:
         solute_atom_pairs.remove(elem)

      self.nonbonded_list = solute_atom_pairs
    ################################################################################################
    def calculate_derivative(self):
      """
      Calculate the first-order partial derivatives of the solute
      """
      # Read the coordinate file frame by frame
      num_coords = 3 * len(self.allparam)
      #self.frames = 100 # for testing purpose
      num_lines_per_frame = num_coords/10       

      f = open(self.fcoord,  'r')
      f.next()     
      joined_lines = []    # every element is a single coordinate
      atomcoord = []       # a 3D list that contains x, y, z coordinate of one atom
      allcoord = []        # a list that contains atom coord of all atoms in one frame
   
      if os.path.isfile(self.output):
          shutil.copyfile(self.output, self.output+'.backup')
          os.remove(self.output)   
      fout = open(self.output, 'a')
      for i in range(len(self.all_types)):
         fout.write('%s,%s  %s,%s  '%(self.all_types[i],'radius', self.all_types[i], 'epsilon'))
      fout.write('\n')   

      for frame in range(self.frames):
        lines_in_a_frame = [next(f) for i in range(num_lines_per_frame+1)]
        lines_in_a_frame = [line.strip(' \t\n\r') for line in lines_in_a_frame]

        # Merge all the coords into a big list
        for i in range(len(lines_in_a_frame)):
          splitline =  lines_in_a_frame[i].split()
	  for coord in splitline:
            joined_lines.append(float(coord))

        # Split the list into subsets (each element now contains X Y Z coordinates of one atom) 
        for coord in range(0, num_coords, 3):
           atomcoord = joined_lines[coord:coord+3]  
           allcoord.append(atomcoord)
        # Parsing the box lengths  
        line = f.next()
        splitline = line.split()
        dim_x, dim_y, dim_z = splitline[0], splitline[1], splitline[2]                                  
        dim = [dim_x, dim_y, dim_z]

        ##################             Derivative calculation starts                #######################################
        # Compute the first-order analytical derivative of free energy with respect to the vdw parameters of each atom type
        # Derivatives of binding enthalpies will come next
        self.all_der = []
        for i in range(len(self.all_types)):
            self.all_der.append([0.0,0.0])

        # Compute the derivatives in all 1-4 pairs
        for j in range(len(self.all_types)):
            for i in range(len(self.pair_list)):
              atom1 = self.pair_list[i][0]
              atom2 = self.pair_list[i][1]
              if self.atom_types[atom1] == self.all_types[j] or self.atom_types[atom2] == self.all_types[j]:
                dist = sensitivity_calc_coord.coords_calc(self.pair_list[i], allcoord, dim)
                if dist <= self.vdw_cutoff:
                  if self.atom_types[atom1] == self.all_types[j]:   
                      rad_tmp, eps_tmp = calc_der(dist, 'atom1', self.allparam[atom1][0], self.allparam[atom2][0], self.allparam[atom1][1], self.allparam[atom2][1], 0.5)
                      self.all_der[j][0] += rad_tmp
                      self.all_der[j][1] += eps_tmp 
                  elif self.atom_types[atom2] == self.all_types[j]:
                      rad_tmp, eps_tmp = calc_der(dist, 'atom2', self.allparam[atom1][0], self.allparam[atom2][0], self.allparam[atom1][1], self.allparam[atom2][1], 0.5)
       	       	      self.all_der[j][0] += rad_tmp
       	       	      self.all_der[j][1] += eps_tmp

        # Compute the derivatives in non 1-4 interactions for the solute
        for j in range(len(self.all_types)):
          for i in range(len(self.nonbonded_list)):
            atom1 = self.nonbonded_list[i][0]
            atom2 = self.nonbonded_list[i][1]
            if self.atom_types[atom1] == self.all_types[j] or self.atom_types[atom2] == self.all_types[j]:
              dist = sensitivity_calc_coord.coords_calc(self.nonbonded_list[i], allcoord, dim)
              if dist <= self.vdw_cutoff:
                  if self.atom_types[atom1] == self.all_types[j]:
                      rad_tmp, eps_tmp = calc_der(dist, 'atom1', self.allparam[atom1][0], self.allparam[atom2][0], self.allparam[atom1][1], self.allparam[atom2][1], 1)
       	       	      self.all_der[j][0] += rad_tmp
       	       	      self.all_der[j][1] += eps_tmp
       	       	  elif self.atom_types[atom2] == self.all_types[j]:
                      rad_tmp, eps_tmp = calc_der(dist, 'atom2', self.allparam[atom1][0], self.allparam[atom2][0], self.allparam[atom1][1], self.allparam[atom2][1], 1)
                      self.all_der[j][0] += rad_tmp
                      self.all_der[j][1] += eps_tmp


        # Compute the derivatives in the interactions between the solute and the solvent
        for m in range(len(self.all_types)):
       	  for j in range(3, self.solute_atoms):
            if self.atom_types[j] == self.all_types[m]: 
              for i in range(self.solute_atoms, self.total_atoms, 3 ): # Loop over all oxygen atoms (hw atom type does not have van der Waals)
                dist = sensitivity_calc_coord.coords_calc([i,j], allcoord, dim)
                if dist <= self.vdw_cutoff:
                  rad_tmp, eps_tmp = calc_der(dist, 'atom1', self.allparam[j][0], self.tip3p_param[0], self.allparam[j][1], self.tip3p_param[1], 1)
                  self.all_der[m][0] += rad_tmp
                  self.all_der[m][1] += eps_tmp
        for i in range(len(self.all_der)):
          fout.write('%8.3f%8.3f'%(self.all_der[i][0], self.all_der[i][1]))
        fout.write('\n')
        del joined_lines[:]
        del allcoord[:]
      return 0            

def calc_der(distance, matching_atom, radius_atom1, radius_atom2, epsilon_atom1, epsilon_atom2, scaling_factor):
    radius = radius_atom1 + radius_atom2
    epsilon = math.sqrt(epsilon_atom1 * epsilon_atom2)
    radius_to_dist = radius/distance
    radius_to_dist_pow5 = radius_to_dist**5
    radius_to_dist_pow6 = radius_to_dist_pow5*radius_to_dist
    rad_der  =  scaling_factor*12.0*epsilon*radius_to_dist_pow5*(radius_to_dist_pow6 - 1.0)/distance  # A scaling factor of 0.5 for 1-4 van der Waals
    if matching_atom	== 'atom1': 
       	eps_der  =  scaling_factor*epsilon_atom2*radius_to_dist_pow6*(0.5*radius_to_dist_pow6 - 1.0)/epsilon
    else:
       	eps_der  =  scaling_factor*epsilon_atom1*radius_to_dist_pow6*(0.5*radius_to_dist_pow6 - 1.0)/epsilon
    
    return rad_der, eps_der


def welcome_message():
    print('*************************************************************************************')
    print(' Welcome to Sensitivity Analysis:')
    print(' a program that calculates the first-order partial derivatives from MD trajectories.')
    print('                      Version 1.0 Beta\n')
    print(' Written by: Jian (Jane) Yin')
    print(' Copyright (c) 2016-2017, University of California, San Diego')
    print('**************************************************************************************')
    
def help_message():
    print('Usage:\n')
    print('    -i  Sensitivity input file')
    print('    -r  mol2 file of the receptor')
    print('    -l  mol2 file of the ligand')
    print('    -c  coordinate file of the complex')
    print('    -p  parameter file of the general force field')
    print('    -w  parameter file of water model/counterions')
    print('    -o  Sensitiivty output file (optional; default file name will be provided)\n')
    print('For example:\n')
    print('     python2 sensitivity.py -i sensitivity.in -r receptor.mol2 -l ligand.mol2 -c mdcrd -p gaff.dat -w water.dat -o data.dat')
    print('or')
    print('     python2 sensitivity.py -i sensitivity.in -r oah.mol2 -l oah_g2.mol2 -c 750000frames.inpcrd -p gaff.dat -w tip3p.dat\n')
              
def main():
    welcome_message()
    this = Sensitivity()
    this.check_arguments()
    this.process_input_file()
    this.get_parameters()
    this.get_connectivity()    
    this.calculate_derivative()

if __name__ == "__main__":
    main()
