# =============================================================================
# randomselect.py - Auxiliary program to randomly select records
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
# The Original Software is "randomselect.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module randomselect.py - Auxiliary program to randomly select records.

   USAGE:
     python randomselect.py [in_file] [out_file] -perc [percentage_value]

     or

     python randomselect.py [in_file] [out_file] -num [num_records]

   ARGUMENTS:
     in_file   Name of the input file with the original data records
     out_file  Name of the output file where the randomly selected records are
               written to
     Either '-perc' or '-num' (plus a value) has to be given as third command
     line argument:
     -perc [percentage_value]  Set the percentage value of how many records
                               should be selected randomly as a floating-point
                               number. The percentage value must be between
                               0.0 and 100.0.
                               For example '-perc 0.1' selects 0.1% of all
                               records.
     -num [num_records]        Alternatively, the absolute number of records
                               to be selected randomly can be given as an
                               argument. The value must be positive.

   DESCRIPTION:
     This program reads in a data file and randomly selects record according
     to the selected optional argument. It writes these records unchanged into
     the output file.
"""

# -----------------------------------------------------------------------------

import sys
import random

# -----------------------------------------------------------------------------

def selectrandom():
  """Main routine, open file, read lines, select randomly, write if select.

  USAGE:
    selectrandom()

  ARGUMENTS:
    None

  DESCRIPTION:
    Main routine, see description of module above.
  """

  # Process command line arguments and check for correctness  - - - - - - - - -
  #
  if (len(sys.argv) != 5):
    print '***** Error: %s needs three arguments:' % (sys.argv[0])
    print '*****        - Name of the original input data file'
    print '*****        - Name of the output file with the selected records'
    print '*****        Then either'
    print '*****          -perc [percentage_value]'
    print '*****        or'
    print '*****          -num  [num_records]'
    raise Exception()

  if (sys.argv[1] == sys.argv[2]):
    print '***** Error: Input and output files must differ'
    print '*****        Input file name:  %s' % (sys.argv[1])
    print '*****        Output file name: %s' % (sys.argv[2])
    raise Exception()

  in_file_name = sys.argv[1]
  out_file_name = sys.argv[2]

  if (sys.argv[3][:2] == '-p'):
    select_mode = 'perc'
  elif(sys.argv[3][:2] == '-n'):
    select_mode = 'num'
  else:
    print '***** Error: Illegal random selection argument: %s' % \
          (sys.argv[3])
    print '*****        Possible are:'
    print '*****          -perc [percentage_value]'
    print '*****        or'
    print '*****          -num  [num_records]'
    raise Exception()

  if (select_mode == 'perc'):
    perc_val = float(sys.argv[4])
    if (perc_val <= 0.0) or (perc_val >= 100.0):
      print '***** Error: Illegal value for random percentage: %s' % \
            (sys.argv[4])
      print '*****        Value must be between 0.0 and 100.0'
      raise Exception()

  else:  # Number of records given
    num_rec = int(sys.argv[4])
    if (num_rec <= 0):
      print '***** Error: Illegal value for number of records: %s' % \
            (sys.argv[4])
      print '*****        Value must be positive.'
      raise Exception()

  # Open input file and check number of available records - - - - - - - - - - -
  #
  try:
    f_in = open(in_file_name,'r')
  except:
    print '***** Error: Can not open input file: %s' % (str(in_file_name))
    raise IOError()

  print "Counting lines in source file..."
  line_count = 0
  for line in f_in.xreadlines():
    line_count += 1
  f_in.close()

  if (select_mode == 'num') and (num_rec > line_count):  # Illegal value
    print '***** Error: Illegal values for number of records: %i' % (num_rec)
    print '*****        File only contains %i lines/records' % (line_count)
    raise Exception()

  # Flag records which have been selected randomly  - - - - - - - - - - - - - -
  #
  if (select_mode == 'perc'):  # Compute number of records to be selected
    num_rec = int(perc_val * float(line_count) / 100.0)

  # Create a sequence of bit flags, one for each record in the file
  #
  record_flags = [0] * line_count

  selected_rec_num = 0

  while (selected_rec_num < num_rec):
    random_rec_num = random.randrange(0,line_count,1)

    if (record_flags[random_rec_num] == 0):
      record_flags[random_rec_num] = 1
      selected_rec_num += 1

  print 'Selected %i records' % (selected_rec_num)

  # Open files, read lines and write selected records to output file  - - - - -

  f_in = open(in_file_name,'r')

  try:
    f_out = open(out_file_name,'w')
  except:
    print '***** Error: Can not write to output file: %s' % \
          (str(out_file_name))
    f_in.close()
    raise IOError()

  increment = int(float(line_count) / 100.0)
  
  print 'Processing source file...'   # - - - - - - - - - - - - - - - - - - - -

  line_count = 0
  percent_complete = 0

  for line in f_in.xreadlines():
    if (record_flags[line_count] == 1):  # This line is selected
      f_out.write(line)

    line_count += 1
    if ((line_count % increment) == 0):
      percent_complete += 1
      print '%3d%% completed' %(percent_complete),'\r',

  f_in.close()
  f_out.close()

  print

# ----------------------------------------------------------------------------

# Start main routine
#
selectrandom()

# ----------------------------------------------------------------------------

