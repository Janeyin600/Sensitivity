import math
from functools import reduce

def coords_calc(pair_list, coord_list, dim_list):
   """
   Compute the distance between a pairs of atoms provided in pair_list
   """
   atom1 = pair_list[0]  # coordinates of atom one
   atom2 = pair_list[1]  # coordinates of atom two
  
   diff =  [coord_atom1 - coord_atom2 for coord_atom1, coord_atom2 in zip(coord_list[atom1],coord_list[atom2])]
   diff = [boundary_conditions(elem_diff, dim_elem) for elem_diff, dim_elem in zip (diff, dim_list)] 
   
   squared = list(map(lambda x: x**2, diff))
   dist = math.sqrt(sum(squared))
   return dist

def boundary_conditions(distance, edge_length):
    """   
    This code is to deal with the period boundary conditions by putting pairs of atoms near the edges of the primary box
    back within the cutoff distance.
    """
    if (distance <= -float(edge_length)*0.5):
      distance += float(edge_length)
    elif (distance > float(edge_length)*0.5):
      distance -= float(edge_length);
    return distance
