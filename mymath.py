# =============================================================================
# mymath.py - Various mathematical routines.
#
# Freely extensible biomedical record linkage (Febrl) Version 0.2.1
# See http://datamining.anu.edu.au/projects/linkage.html
#
# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.1
#
# The contents of this file are subject to the ANUOS License Version 1.1 (the
# "License"); you may not use this file except in compliance with the License.
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
# the specific language governing rights and limitations under the License.
# The Original Software is "mymath.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module mymath.py - Various mathematical routines.

   See doc strings of individual functions for detailed documentation.
"""

# =============================================================================
# Imports go here

import math

# =============================================================================

def mean(x):
  """Compute the mean (average)  of a list of numbers.
  """

  if (len(x) == 1):  # Only one element in list
    return float(x[0])

  elif (len(x) == 0):  # Empty list
    print 'warning:Empty list given: %s' % (str(x))
    return None

  else:  # Calculate average
    sum = 0.0
    for i in x:
      sum += i

    res = sum / float(len(x))

    return res

# =============================================================================

def stddev(x):
  """Compute the standard deviation of a list of numbers.
  """

  if (len(x) == 1):  # Only one element in list
    return 0.0

  elif (len(x) == 0):  # Empty list
    print 'warning:Empty list given: %s' % (str(x))
    return None

  else:
    sum = 0.0
    for i in x:
      sum += i

    avrg = sum / float(len(x))

    sum = 0.0
    for i in x:
      sum = sum + (i - avrg) * (i - avrg) 

    res = math.sqrt(sum / float(len(x)))

    return res

# =============================================================================

def log2(x):
  """Compute binary logarithm (log2) for a floating-point number.

  USAGE:
    y = log2(x)

  ARGUMENT:
    x  An positive integer or floating-point number

  DESCRIPTION:
    This routine computes and returns the binary logarithm of a positive
    number.
  """

  return math.log(x) / 0.69314718055994529  # = math.log(2.0)    

# =============================================================================

def perm_tag_sequence(in_tag_seq):
  """Create all permuations of a tag sequence.

  USAGE:
    seq_list = perm_tag_sequence(in_tag_seq)

  ARGUMENT:
    in_tag_seq  Input sequence (list) with tags

  DESCRIPTION:
    This routine computes all permutations of the given input sequence. More
    than one permutation is created if at least one element in the input
    sequence contains more than one tag.

    Returns a list containing tag sequences (lists).
  """

  if (not isinstance(in_tag_seq, list)):
    print 'error:Input tag sequence is not a list: %s' % (str(in_tag_seq))

  list_len = len(in_tag_seq)
  out_tag_seq = [[]]  # List of output tag sequences, start with one empty list

  for elem in in_tag_seq:
    if ('/' in elem):  # Element contains more than one tag, covert into a list
      elem = elem.split('/')

    tmp_tag_seq = []

    if (isinstance(elem,str)):  # Append a simple string
      for t in out_tag_seq:
        tmp_tag_seq.append(t + [elem])  # Append string to all tag sequences

    else:  # Process a list (that contains more than one tags)
      for tag in elem:  # Add each tag in the list to the temporary tag list
        for t in out_tag_seq:
          tmp_tag_seq.append(t+[tag])  # Append string to all tag sequences

    out_tag_seq = tmp_tag_seq

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Input tag sequence: %s' % (str(in_tag_seq))
  print '3:  Output permutations:'
  for p in out_tag_seq:
    print '3:    %s' % (str(p))

  return out_tag_seq

# =============================================================================
