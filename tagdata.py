# =============================================================================
# tagdata.py - Main module to tag elements in a training data file
#
# Freely extensible biomedical record linkage (Febrl) Version 0.2.2
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
# The Original Software is "tagdata.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module tagdata.py - Main module to tag elements into a training data file

   DESCRIPTION:
     Ths module can be used to tag records randomly selected from a data set to
     be used to train a hidden Markov model (HMM) using the 'trainhmm.py'
     module. Records are tagged using either lookup-tables for names or for
     address words.

     The 'tag_component' can be either 'name' or 'address', and is used to 
     determine the possible tags to be used. The corresponding list of field
     names has to be given in 'tag_component_fields'.

     The input data set has to be defined as a Febrl data set implementation in
     'read' access mode.

     Records are selected randomly in the range between (and including) the
     'start_rec_number' and 'end_rec_number'. A total of 'num_rec_to_select'
     records will be selected, with their component (name or address) being not
     empty.

     The tagged records will be written to the file given in 'output_file_name'
     with the original component (name or address) being listed as a comment
     line (starting with the Python comment character hash '#')  and the line
     number from the input file (starting with zero) at the beginning of a
     line. After each input data line, one or more lines with tag sequences
     follows.

     Word spilling can be checked if the 'check_word_spilling' flag is set to
     'True' and the 'field_separator' (the character or string which is
     inserted between the original field values) is not set to an empty string
     (e.g. field_separator = ' ').

     If a hidden Markov model (HMM) file is defined using 'hmm_file_name' the
     training records are tagged as well as standardised using this HMM. This
     allows a semi-automatic training process, where the user only has to
     inspect the output training file and change HMM states or tags for cases
     that are standardised incorrectly. This mechanism reduces the time needed
     to create enough records to train a HMM  training.

     If you want to re-process an existing tagged and corrected training file,
     for example after making some adjustments to some tagging look-up tables,
     then set the 'retag_file_name' to the name of the file to be re-tagged.
     Note that the retagged file will be written to the output file name.

     When a hidden Markov model file is defined, it is pssible to save
     probability from the HMM segmentation process into a 'freqs_file_name' if
     defined (not set to None). Frequencies of tag/state patterns are written
     in ascending order into the given file. This is useful for finding
     examples of unusual patterns of names or addresses which might need to be
     added to the training file(s).

     The user has to manually inspect the output file and delete (or comment)
     all lines with tags that are not correct, and insert a (or modify the
     given) HMM state name for each tag in a sequence.

   EXAMPLE:
     The three input data records:

       "dr peter baxter dea"
       "miss monica mitchell meyer"
       "phd tim william jones harris"

     will be processed and written into the output file as:

       # 0: |dr peter baxter dea|
         TI:, GM:, GM:, GF:
         TI:, GM:, SN:, GF:
         TI:, GM:, GM:, SN:
         TI:, GM:, SN:, SN:

       # 1: |miss monica mitchell meyer|
         TI:, UN:, GM:, SN:
         TI:, UN:, SN:, SN:

       # 2: |phd tim william jones harris|
         TI:, GM:, GM:, UN:, SN:

     If the '-hmm' option is given the output will be something like this:

       # 0: |dr peter baxter dea|
         TI:titl, GM:gname1, GM:gname2, GF:sname1
         TI:titl, GM:gname1, SN:sname1, GF:sname2
         TI:titl, GM:gname1, GM:gname2, SN:sname1
         TI:titl, GM:gname1, SN:sname1, SN:sname2

       # 1: |miss monica mitchell meyer|
         TI:titl, UN:gname1, GM:sname1, SN:sname2
         TI:titl, UN:gname1, SN:sname1, SN:sname2

       # 2: |phd tim william jones harris|
         TI:titl, GM:gname1, GM:gname2, UN:sname1, SN:sname2
"""


# =============================================================================
# Imports go here

from febrl import *            # Main Febrl classes
from dataset import *          # Data set routines
from standardisation import *  # Standardisation routines
from lookup import *           # Look-up table routines
from mymath import *           # Mathematical routines
from name import *             # Name tagging functionality
from address import *          # Address tagging functionality
from simplehmm import *        # Hidden Markov model (HMM) functionalities

import os      # Operating system depedent functions
import random  # Python standard module for random number functionality
import time    # Python standard module for time functions

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

tag_febrl = Febrl(description = 'Data tagging Febrl instance',
                   febrl_path = '.')

tag_project = tag_febrl.new_project(name = 'Tag-Data',
                             description = 'Data tagging module',
                               file_name = 'tag.fbr')

# =============================================================================
# Define a project logger

tag_log = ProjectLog(file_name = 'tag_data.log',
                       project = tag_project,
                     log_level = 1,
                 verbose_level = 1,
                     clear_log = True,
                       no_warn = False)

# =============================================================================
# Define settings for data tagging

# Define your original input data set - - - - - - - - - - - - - - - - - - - - -
#
input_data = DataSetCSV(name = 'example1in',
                 description = 'Example data set 1',
                 access_mode = 'read',
                header_lines = 1,
                   file_name = './dbgen/dataset1.csv',
                      fields = {'rec_id':0,
                                'given_name':1,
                                'surname':2,
                                'street_num':3,
                                'address_part_1':4,
                                'address_part_2':5,
                                'suburb':6,
                                'postcode':7,
                                'state':8,
                                'date_of_birth':9,
                                'soc_sec_id':10},
              fields_default = '',
                strip_fields = True,
              missing_values = ['','missing'])

# Define block of records to be used for tagging  - - - - - - - - - - - - - - -
#
start_rec_number = 0
end_rec_number =   input_data.num_records

# Define number of records to be selected randomly  - - - - - - - - - - - - - -
#
num_rec_to_select = (input_data.num_records / 2)  # Use 50% for tagging

# Define name of output data set  - - - - - - - - - - - - - - - - - - - - - - -
#
output_file_name = 'tagged_data.csv'

# Component: Can either be 'name' or 'address'  - - - - - - - - - - - - - - - -
#
tag_component = 'name'

# Define a list with field namess from the input data set in the component  - -
#
tag_component_fields = ['given_name', 'surname']

# Define if word spilling should be checked or not  - - - - - - - - - - - - - -
#
check_word_spilling = True  # Set to True or False

# Define the field separator  - - - - - - - - - - - - - - - - - - - - - - - - -
#
field_separator = ' '

# Use HMM for tagging and segmenting  - - - - - - - - - - - - - - - - - - - - -
#
hmm_file_name = './hmm/name-absdiscount.hmm' # Set to name of a HMM file or
                                             # None

# Retag an existing training file - - - - - - - - - - - - - - - - - - - - - - -
# - Note that retagging is only possible if a HMM file name is given as well
# - If the retag file name is defined, the start and end record numbers as
#   defined above are not used, instead the record numbers in the re tag file
#   are used.
#
retag_file_name = None  # Set to name of an existing training file or None

# Write out frequencies into a file - - - - - - - - - - - - - - - - - - - - - -
#
freqs_file_name = 'freqs.txt'  # Set to a file name or None

# Define and load lookup tables - - - - - - - - - - - - - - - - - - - - - - - -
#
name_lookup_table = TagLookupTable(name = 'Name lookup table',
                                default = '')

name_lookup_table.load(['./data/givenname_f.tbl',
                        './data/givenname_m.tbl',
                        './data/name_prefix.tbl',
                        './data/name_misc.tbl',
                        './data/saints.tbl',
                        './data/surname.tbl',
                        './data/title.tbl'])

name_correction_list = CorrectionList(name = 'Name correction list')

name_correction_list.load('./data/name_corr.lst')

address_lookup_table = TagLookupTable(name = 'Address tagging lookup table',
                                   default = '')

address_lookup_table.load(['./data/country.tbl',
                          './data/address_misc.tbl',
                          './data/address_qual.tbl',
                          './data/institution_type.tbl',
                          './data/locality_name_act.tbl',
                          './data/locality_name_nsw.tbl',
                          './data/post_address.tbl',
                          './data/postcode_act.tbl',
                          './data/postcode_nsw.tbl',
                          './data/saints.tbl',
                          './data/territory.tbl',
                          './data/unit_type.tbl',
                          './data/wayfare_type.tbl'])

address_correction_list = CorrectionList(name = 'Address correction list')

address_correction_list.load('./data/address_corr.lst')

# =============================================================================
# =============================================================================
# Do not change anything below here
# =============================================================================
# =============================================================================

# Test settings - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
if (end_rec_number <=  start_rec_number):
  print 'error:End record number is larger or equal to start record number'
  raise Exception

rec_range = end_rec_number - start_rec_number + 1  # Range of racords

if (rec_range < num_rec_to_select):
  print 'error:Range of input records smaller than records to be selected'
  raise Exception

if (tag_component not in ['name', 'address']):
  print 'error:Illegal tag component: %s' % (str(tag_component))
  raise Exception

if (check_word_spilling not in [True, False]):
  print 'error:Illegal value for "check_word_spilling", must be True or False'
  raise Exception

# Load HMM if file name is given - - - - - - - - - - - - - - - - - - - - - - -

if (hmm_file_name != None):
  if (hmm_file_name == input_data.file_name):
    print 'error:HMM file name is the same as input file name'
    raise Exception
  elif (hmm_file_name == output_file_name):
    print 'error:HMM file name is the same as output file name'
    raise Exception

  the_hmm = hmm('Tagging HMM',[],[])  # Create new empty HMM object
  the_hmm.load_hmm(hmm_file_name)
  the_hmm.print_hmm()  # Print HMM (according to verbose and logging level)

else:
  the_hmm = None  # Define non existing HMM

# If the retag file name is given load an existing training file for re-tagging
#
if (retag_file_name != None):
  if (the_hmm == None):  # Make sure a HMM is given
    print 'error:Re-tagging only possible when a HMM is given'
    raise Exception

  if (retag_file_name == input_data.file_name):
    print 'error:Re-tag file name is the same as input file name'
    raise Exception
  elif (retag_file_name == output_file_name):
    print 'error:Re-tag file name is the same as output file name'
    raise Exception
  elif (retag_file_name == hmm_file_name):
    print 'error:Re-tag file name is the same as HMM file name'
    raise Exception

  try:
    f_in = open(retag_file_name, 'r')  # Test if file is available
  except:
    print 'error:Cannot open file specified for re-tagging: %s' % \
          (str(retag_file_name))
    raise IOError

  # Now gather record numbers and previous tags/states, as well as the original
  # header information. Use a simple state machine to do this.
  #
  tagged_recs  = {}
  cleaned_recs = {}
  original_header_lines = []
  state = -1  # Header lines state
  prevline = ''

  for line in f_in.xreadlines():  # Read training file and process it
    line = line.strip()

    if (state == -1) and (len(line) == 0):  # End of header lines
      state = 0
      prevline = line
      continue

    if (state == -1) and (len(line) > 0) and (line[0] == "#"):
      original_header_lines.append("# " + line)
      prevline = line
      continue
    sline = line.split(' ')

    if (len(sline) > 2) and (len(sline[2]) > 3) and (sline[0] == '#') \
       and (sline[2][0] == '(') and (sline[2][-2:] == '):'):
      try:	
        rec = int(sline[1])  # Original record number 
        tagged_recs[rec]  = None
        cleaned_recs[rec] = None
        state = 1
      except:
        pass
      prevline = line
      continue

    if (state == 1) and (len(line) > 0) and (line[0] != '#'):
      tagged_recs[rec]  = line
      cleaned_recs[rec] = prevline
      state = 0
      prevline = line
      continue

    if (state == 1) and (len(line) > 0):
      prevline = line
      continue

  f_in.close()
  tagged_recs_keys = tagged_recs.keys()

  start_rec_number =  0  # Override specified numbers
  end_rec_number =    input_data.num_records  
  num_rec_to_select = len(tagged_recs_keys)

# If a frequency file name is given, check  - - - - - - - - - - - - - - - - - -
#
if (freqs_file_name != None):
  if (the_hmm == None):  # Make sure a HMM is given
    print 'error:Frequency file output only possible when a HMM is given'
    raise Exception

  if (freqs_file_name == input_data.file_name):
    print 'error:Frequency file name is the same as input file name'
    raise Exception
  elif (freqs_file_name == output_file_name):
    print 'error:Frequency file name is the same as output file name'
    raise Exception
  elif (freqs_file_name == hmm_file_name):
    print 'error:Frequency file name is the same as HMM file name'
    raise Exception

# Print header (and write it to the log file if activated)  - - - - - - - - - -
#
print '1:'
print '1:'+'#'*75
print '1:#'
print '1:# "tagdata.py" - Version 0.2'
print '1:# Process started at: %s' % (time.ctime(time.time()))
print '1:#'
print '1:# Input file name:        %s' % (str(input_data.file_name))
print '1:# Output file name:       %s' % (str(output_file_name))
print '1:# Tag component:          %s' % (str(tag_component))
if (hmm_file_name != None):
  print '1:# Using HMM:              %s' % (str(hmm_file_name))
  if (retag_file_name != None):
    print '1:# Re-tagging file:        %s' % (str(retag_file_name))
  if (freqs_file_name != None):
    print '1:# Writing frequency file: %s' % (str(freqs_file_name))
print '1:#'

# Open output file and write header - - - - - - - - - - - - - - - - - - - - - -
#
try:
  f_out = open(output_file_name,'w')
except:
  print 'error:Cannot write to output file: %s' % (str(output_file_name))
  raise IOError

f_out.write('#'+'#'*70+os.linesep)
f_out.write('# Tagged training data written by "tagdata.py" - Version 0.2'+ \
            os.linesep)
f_out.write('#'+os.linesep)
f_out.write('# Created '+time.ctime(time.time())+os.linesep)
f_out.write('#'+os.linesep)
f_out.write('# Input file name:  '+input_data.file_name+os.linesep)
f_out.write('# Output file name: '+output_file_name+os.linesep)
f_out.write('# Tag component:    '+tag_component+os.linesep)
f_out.write('# Parameters:'+os.linesep)
f_out.write('# - Start of block with training records: '+ \
            str(start_rec_number)+os.linesep)
f_out.write('# - End of block with training records:   '+ \
            str(end_rec_number)+os.linesep)
f_out.write('# - Number of training records:           '+ \
            str(num_rec_to_select)+os.linesep)
if (hmm_file_name != None):
  f_out.write('# Using HMM for standardisation: '+hmm_file_name+os.linesep)
  if (retag_file_name != None):
    f_out.write('# Re-tagging file: '+retag_file_name+os.linesep)
    f_out.write('#   Header lines from original training file follow:'+ \
                os.linesep)
    for header_line in original_header_lines:
      f_out.write(header_line + os.linesep)
  if (freqs_file_name != None):
    f_out.write('# Writing tag/state pattern to frequency file: '+ \
                freqs_file_name+os.linesep)
f_out.write('#'+os.linesep)
f_out.write('#'+'#'*70+os.linesep)
f_out.write(os.linesep)

# Define various variables  - - - - - - - - - - - - - - - - - - - - - - - - - -
#
rec_count    = 0                  # Number of selected records
num_rec_left = num_rec_to_select  # Number of records to be selected left
rec_selected = {}                 # Dictionary of all record numbers that were
                                  # selected
seq_freqs = {}                    # Dictionary to hold examples of tag/state
                                  # patterns

unchanged_loop_cnt = 0  # Counter of how many loops have been done without new
                        # records being selected

prev_num_rec_left = num_rec_to_select  # Number of records left in the previous
                                       # interation

# Due to the random nature of selecting records, and because sometimes  - - -
# a selected component can be empty (and is thus not used for training)
# more than one iteration over the input data set is carried out. In each 
# iteration, records are selected randomly.
#
while (rec_count < num_rec_to_select):  # Loop until 'num_rec_to_select'
                                        # records are selected

  # Read first record in the given range  - - - - - - - - - - - - - - - - - - -
  #
  record = input_data.read_records(start_rec_number,1)[0]

  line_read = start_rec_number  # Number of the line read

  while (rec_count < num_rec_to_select) and (line_read <= end_rec_number):

    if (record != None) and \
       (((retag_file_name != None) and (line_read in tagged_recs_keys)) or \
        ((retag_file_name == None) and \
         (num_rec_left >= random.randrange(0,rec_range, 1)))):

      record_id =  '[RecID: '+str(record['_rec_num_'])+'/'+ \
                   record['_dataset_name_']+'] '
      fields_str = os.linesep+'###            [Fields: '+str(record)+']'

      # Now concatenate field values into a string  - - - - - - - - - - - - - -
      #
      component = ''

      for field in tag_component_fields:

        field_value = record.get(field,'')  # Get the field value

        if (field_value != ''):

          field_sep = field_separator  # Get original field separator value

          # Check field spilling only if field separator is not an empty string
          #
          if (field_separator != '') and (check_word_spilling == True):

            if (tag_component == 'name'):
              spill_flag = check_field_spilling(component, field_value, \
                                                name_lookup_table, record_id,
                                                fields_str)
            else:  # Address component
              spill_flag = check_field_spilling(component, field_value, \
                                                address_lookup_table,
                                                record_id, fields_str)
            if (spill_flag == True):
              field_sep = ''  # A word spills over, so make field separator an
                              # empty string

          # Append separator and field to the component
          #
          if (component == ''):  # Start of component, no need for a separator
            component = field_value
          else:
            component = component+field_sep+field_value

      if (tag_component == 'name'):
        clean_comp = clean_component(component, name_correction_list, \
                                     record_id)
      else:   # Address component
        clean_comp = clean_component(component, address_correction_list, \
                                     record_id)

      # Now tag the component - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      if (clean_comp != '') and (not rec_selected.has_key(line_read)):

        print '3:    Cleaned component: "%s"' % (clean_comp)

        if (tag_component == 'name'):
          [word_list, tag_list] = tag_name_component(clean_comp, \
                                                     name_lookup_table,
                                                     record_id)
        else:   # Address component
          [word_list, tag_list] = tag_address_component(clean_comp, \
                                                        address_lookup_table,
                                                        record_id)

        if (tag_list != []):  # Only process non-empty tag lists

          # Append record number into dictionary of processed records
          #
          rec_selected.update({line_read:line_read})

          # Create all permutation sequences of this tag list - - - - - - - -
          #
          tag_seq = perm_tag_sequence(tag_list)  # From mymath module

          print '3:      Word list: %s' % (str(word_list))
          print '3:      Tag list:  %s' % (str(tag_list))
          print '3:      Tag permutations: %s' % (str(tag_seq))

          # Do HMM processing - - - - - - - - - - - - - - - - - - - - - - - - -
          #
          if (hmm_file_name != None):

            state_seq  = []    # List containing computed HMM state sequences
            max_prob   = -1.0  # maximal probability for a sequence
            max_seq_no = -1    # Number of the seq. with the max. probablity

            # Now give tag sequences to the HMMs to compute state sequences
            #
            i = 0
            for t in tag_seq:
              [obs_seq, prob] = the_hmm.viterbi(t)
              state_seq.append(obs_seq)
              if (prob > max_prob):
                max_prob = prob
                max_seq_no = i
              i += 1

          # Write original component and resulting tag sequences to output - -
          #
          f_out.write('# '+str(line_read)+' ('+str(rec_count)+'): |'+ \
                      component+'|'+os.linesep) # Commented original
          num_len = len(str(line_read))+len(str(rec_count))+6

          f_out.write('#'+num_len*' '+'|'+' '.join(word_list)+'|'+os.linesep)

          for i in range(len(tag_seq)):
            # Convert each tag sequence into a string for file output
            #
            seq_string = '  '

            if (hmm_file_name != None) and (i != max_seq_no):
              seq_string = '# ' # Comment sequences with not max. probability

            for j in range(len(tag_seq[i])):

              if (hmm_file_name != None):
                seq_string = seq_string+' '+tag_seq[i][j]+':'+ \
                             state_seq[i][j]+','
              else:
                seq_string = seq_string+' '+tag_seq[i][j]+':,'

            f_out.write(seq_string[:-1]+os.linesep)  # Write without , at end
            print '3:    %s' % (str(seq_string[:-1]))

          if (hmm_file_name != None):
            f_out.write('# Maximum Viterbi probability: %0.5f'% \
                        (max_prob) + os.linesep)
            print '3:    Maximum Viterbi probability: %0.5f' % (max_prob)

          if (retag_file_name != None) and (tagged_recs[line_read] != None):
            if (tagged_recs[line_read].strip() != seq_string[:-1].strip()):
              f_out.write("# Note: ***** Changed *****" + os.linesep)

              print '3:%s   Note: ***** Changed *****' % (record_id)
              f_out.write('# Was: ' + tagged_recs[line_read]+os.linesep)
                          # Write commented original tag sequence
              print '3:%s    Original tag sequence: %s' % \
                    (record_id, str(tagged_recs[line_read]))

          f_out.write(os.linesep)  # Write an empty line
          print '3:'

          if (hmm_file_name != None):
            seq_key = seq_string[:-1]  # Add sequence to dictionary
            if (seq_freqs.has_key(seq_key)):
              seq_freqs[seq_key].append(['|'+' '.join(word_list)+'|', \
                                         max_prob])
            else:
              seq_freqs[seq_key] = [['|'+' '.join(word_list)+'|',max_prob]]

          rec_count += 1

          # Print process indicator message every 100 records - - - - - - - - -
          #
          if (rec_count % 100 == 0):
            print '1:  Processed %i of %i records' % \
                  (rec_count, num_rec_to_select)

    record = input_data.read_record()  # Read next record
    line_read += 1

  num_rec_left = num_rec_to_select - rec_count

  if (prev_num_rec_left == num_rec_left):  # No new records selected
    unchanged_loop_cnt += 1
  prev_num_rec_left = num_rec_left  # Set to current value

  if (unchanged_loop_cnt > 5):  # Do five loops maximal without selecting
                                # new records
    print 'warning:Can not select more than %i' % (rec_count) + \
          ' records for training. This is probably due to empty input ' + \
          'components. Please reduce the value of "num_rec" or increase ' + \
          'the range between "start_rec_number" and "end_rec_number".'
    break

  if (num_rec_left < 10):  # Only 10 records left to select
    num_rec_left = num_rec_to_select + 1  # Set to more than 100% probablity

  elif (num_rec_left < (num_rec_to_select / 100.0)):  # Less than 1% records
                                                      # left
    num_rec_left = int(num_rec_to_select / 100.0)  # Set to 1%

f_out.close()

# If specified, save Viterbi frequencies to a file  - - - - - - - - - - - - -
#
if (freqs_file_name != None):
  freqs_out = open(freqs_file_name, 'w')  # Open frequency file for writing
  freqs_out.write('# Frequency listing of tag/state patterns written by')
  freqs_out.write('"tagdata.py" - Version 0.2'+os.linesep)
  freqs_out.write('#'+os.linesep)
  freqs_out.write('# Created '+time.ctime(time.time())+os.linesep)
  freqs_out.write('#'+os.linesep)
  freqs_out.write('# Input file name:  '+input_data.file_name+os.linesep)
  freqs_out.write('# Output file name: '+output_file_name+os.linesep)
  freqs_out.write('# Parameters:'+os.linesep)
  freqs_out.write('# - Start of block with training records: '+ \
                  str(start_rec_number)+os.linesep)
  freqs_out.write('# - End of block with training records:   '+ \
                  str(end_rec_number)+os.linesep)
  freqs_out.write('# - Number of training records:           '+ \
                  str(num_rec_to_select)+os.linesep)
  if (hmm_file_name != None):
    freqs_out.write('# Using HMM for standardisation: '+ \
                    hmm_file_name+os.linesep)
    if (retag_file_name != None):
      freqs_out.write('# Re-tagging file: '+retag_file_name+os.linesep)
      freqs_out.write('#   Header lines from original training file follow:'+ \
                      os.linesep)
      for header_line in original_header_lines:
        freqs_out.write(header_line + os.linesep)
  freqs_out.write('#'+os.linesep)
  freqs_out.write('#'+'#'*70+os.linesep)
  freqs_out.write(os.linesep)

  sorted_seq_freqs = []  # Now sort sequences according to their fruequencies
  for key in seq_freqs.keys():
    sorted_seq_freqs.append((len(seq_freqs[key]),key))
  sorted_seq_freqs.sort()

  for skey in sorted_seq_freqs:
    key = skey[1]
    freqs_out.write('# Pattern: '+str(key)+os.linesep)
    freqs_out.write('# Frequency: '+str(skey[0])+os.linesep)
    examples = seq_freqs[key]
    freqs_out.write('# Maximum Viterbi probability: '+ \
                    str(examples[0][1])+os.linesep)
    freqs_out.write('# Examples: '+os.linesep)
    for example in examples:
      freqs_out.write('#    '+str(example[0])+os.linesep)
    freqs_out.write(str(key)+os.linesep)
    freqs_out.write(os.linesep)
  freqs_out.close()

print '1:Read %i lines, processed %i lines' % (line_read, rec_count)

print '1:'
print '1:End.'

# =============================================================================
