# =============================================================================
# mymath.py - Various mathematical routines.
#
# Freely extensible biomedical record linkage (Febrl) Version 0.1
# See http://datamining.anu.edu.au/projects/linkage.html
#
# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.0
#
# The contents of this file are subject to the ANUOS License Version 1.0 (the
# "License"); you may not use this file except in compliance with the License.
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
# the specific language governing rights and limitations under the License.
# The Original Software is "mymath.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module mymath.py - Various mathematical routines.

   See doc strings of individual functions for detailed documentation.
"""

# -----------------------------------------------------------------------------

import math
import types

# import config  # No verbose output or logging needed (so far)
# import inout

# -----------------------------------------------------------------------------

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

  return math.log(x) / 0.69314718055994529 # = math.log(2.0)    

# -----------------------------------------------------------------------------

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

  list_len = len(in_tag_seq)
  out_tag_seq = [[]]  # List of output tag sequences, start with one empty list

  for elem in in_tag_seq:
    if ('/' in elem):  # Element contains more than one tag, covert into a list
      elem = elem.split('/')

    tmp_tag_seq = []

    if (type(elem) == types.StringType):  # Append a simple string
      for t in out_tag_seq:
        tmp_tag_seq.append(t + [elem])  # Append string to all tag sequences

    else:  # Process a list (that contains more than one tags)
      for tag in elem:  # Add each tag in the list to the temporary tag list
        for t in out_tag_seq:
          tmp_tag_seq.append(t+[tag])  # Append string to all tag sequences

    out_tag_seq = tmp_tag_seq

  return out_tag_seq

# -----------------------------------------------------------------------------

def mean(x):
  """Compute the mean of a sequence of numbers.
  """

  sum = 0.0
  for i in x:
    sum += i

  res = sum / float(len(x))

  return(res)

# --------------------------------------------------------------------

def stddev(x):
  """Compute the standard deviation of a sequence of numbers.
  """

  avg = mean(x)

  sum = 0.0
  for i in x:
    sum = sum + (i - avg) * (i - avg) 

  res = math.sqrt(sum / float(len(x)))

  return(res)

# -----------------------------------------------------------------------------

def test():
  """Simple test routine.
  """

  test_tags = [['TI','GM','PR',['GF','SN'],'UN'], \
               [['GM','GF','SN'],'UN','PR',['GF','GM'],['GM','SN']], \
              ]

  for t in test_tags:
    print 'Input tag sequence:', t

    tag_list = perm_tag_sequence(t)

    print 'Output tag sequences:'
    for l in tag_list:
      print ' ', l

# -----------------------------------------------------------------------------
