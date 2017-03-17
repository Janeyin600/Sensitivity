def ismyinstance(variable_type, parameter_value, filename, parameter):
    """
    Check the data type of the user inputs
    """ 
    if not parameter_value:
    # If the value of entry is left as blank    
         if variable_type == 'string':
             return 'None'
         elif variable_type == 'list':
             return []
         elif variable_type == 'float':
             return 0.0
         elif variable_type == 'int':
             return 0
    
    # if input was provided
    if variable_type == 'float':
        try:
            float(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a float value for %s.' % parameter)
            sys.exit()

        if float(parameter_value) < 0:
            print('Wrong input in %s:' % (filename))
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return float(parameter_value)

    elif variable_type == 'int':
        try:
            int(parameter_value)
        except ValueError:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter an int value for %s.' % parameter)
            sys.exit()

        if int(parameter_value) < 0:
            print('Wrong input in %s:' % filename)
            print(parameter_value)
            print('Please enter a non-negative value for %s.' % parameter)
            sys.exit()
        else:
            return int(parameter_value)

    elif (variable_type == 'list') or (variable_type == 'string'):
        return parameter_value

def find_1_4(angle1, angle2):

  if angle1[1] == angle2[0] and angle1[2] == angle2[1]:
    tmp_pair =  [angle1[0], angle2[2]]
  elif angle2[1] == angle1[0] and angle2[2] == angle1[1]:
    tmp_pair = [angle1[2], angle2[0]]
  elif angle1[1] == angle2[0] and angle1[0] == angle2[1]:
    tmp_pair =  [angle1[2], angle2[2]]
  elif angle2[1] == angle1[2] and angle2[2] == angle1[1]:
    tmp_pair = [angle1[0], angle2[0]]
  else:
    tmp_pair = []

  return sorted(tmp_pair)

def remove_duplicates_from_list(a):

  dup = []
  for i in range(len(a)):
    for j in range(i+1,len(a)):
      #print a[i],a[j]
      if a[i] == a[j]:
        dup.append(a[i])

  for x in dup:
    if x in a:
     a.remove(x)
    
  return a
