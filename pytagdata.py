# =============================================================================
# pyTagData.py - Main module to tag elements in a training data file
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
# The Original Software is "pyTagData.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module pyTagData.py - Main module to tag elements into a training data file

   USAGE:
     pyTagData.py [project_module] [tag_mode] [train_file] [first_rec]
                  [last_rec] [num_rec] [options]
   ARGUMENTS:
     project_module  A Python module containing all project settings.
     tag_mode        Mode for tagging, must be either 'name' or 'locality'
     train_file      Name of the output data file containing the tagged
                     training records (as comma separated values text file)
     first_rec       Start of block with training records
     last_rec        End of block with training records
     num_rec         Number of training records to tag and write to the output

   OPTIONS (after arguments);
     -hmm [file_name]   Load and use a Hidden Markov Model (HMM) to tag the
                        training records.
     -retag [file_name] Re-process an existing training file. Used with the
                        '-hmm' option only. Note that the values given for
                        'first_rec', 'last_rec' and 'num_rec' are overridden
                        when the '-retag' option is used (the number of records
                        in the training file to be re-processed is used
                        instead). However, values for 'first_rec', 'last_rec
                        and 'num_rec' must still be given to act as
                        placeholders on the command line
                        (TO-DO: must fix this).
     -freqs [file_name] Compile and write frequencies of tag/state patterns in
                        ascending order into the given file. This is useful
                        for finding examples of unusual patterns of names or
                        addresses which might need to be added to the training
                        file(s). Used with the '-hmm' option only.
     -l [file_name]     Log warnings and errors into the given file
     -v1                Verbose output level 1 (low output)
     -v2                Verbose output level 2 (high output)
     -nowarn            Don't print warning messages

   DESCRIPTION:
     This program reads in a data file as specified in a configuratiom module
     'project.py' and uses the component as specified by the 'tag_mode'
     argument (either 'name' or 'locality') to create a training file with
     tagged records. Records are tagged using either lookup-tables for names
     or for locality words. These tagged records can then be used to train a
     Hidden Markov Model (HMM) with the 'pyTrainHMM.py' program.

     The program skips the input data file until record number 'first_rec' and
     from then on randomly selects 'num_rec' in the block between 'first_rec'
     and 'last_rec'. For example, if 'first_rec' is set to 0 and 'last_rec' is
     set to the number of records in the data file, 'num_rec' records are
     selected randomly over the whole data set.

     If the option '-hmm' followed by the file name of a HMM file is given, the
     training records are tagged as well as standardised using this HMM. This
     allows a semi-automatic training process, where the user only has to
     inspect the output training file and change HMM states or tags for cases
     that are standardised incorrectly. This mechanism reduces the time needed
     to create enough records to train a HMM  training.

     If you want to re-process an existing tagged and corrected training file,
     for example after making some adjustments to some tagging look-up tables,
     then used the '-retag' option, specifying the name of the existing
     training file.

     For verbose output, add the option '-v1' or '-v2' with the second option
     producing more detailed output than the first one.

     To suppress the printing of warning messages used the '-nowarn' option.
     If logging is activated, warning messages are always logged.

     If you want to save warning and error messages into a file, use the '-l'
     option followed by the file name for the log file. If the log file already
     exists, new information will be appended, otherwise the file will be
     created. Logging is done according to the verbose level. If no verbose
     level is given, only warning and error messages will be saved to the log
     file.

     The selected original input data (name or locality component) is written
     to the output file as comment lines with a hash (#) character and the 
     line number from the input file (starting with zero) at the beginning of
     a line. After each input data line, one or more lines with tag sequences
     follows.

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
         TI:titl, GM:gname1, GM:gname2:, GF:sname1:
         TI:titl, GM:gname1, SN:sname1:, GF:sname2:
         TI:titl, GM:gname1, GM:gname2:, SN:sname1:
         TI:titl, GM:gname1, SN:sname1:, SN:sname2:

       # 1: |miss monica mitchell meyer|
         TI:titl, UN:gname1:, GM:sname1, SN:sname2
         TI:titl, UN:gname1:, SN:sname1, SN:sname2

       # 2: |phd tim william jones harris|
         TI:titl, GM:gname1, GM:gname2, UN:sname1, SN:sname2
"""

# -----------------------------------------------------------------------------

import sys
import os
import types
import random
import time

import config
import mymath
import name
import locality
import inout
import simplehmm

# -----------------------------------------------------------------------------

def tagdata():
  """Main routine, open file, read lines, tag data records, write to out-file.

  USAGE:
    tagdata()

  ARGUMENTS:
    None

  DESCRIPTION:
    Main routine, see description of module above.
  """

  # Process command line arguments and check for correctness  - - - - - - - - -
  #
  if (len(config.options) < 5):
    print '***** Error: %s needs at least six arguments:'% (sys.argv[0])
    print '*****        - Name of the project module'
    print '*****        - Tagging mode: "name" or "locality"'
    print '*****        - Output training file name'
    print '*****        - Start of block with training records'
    print '*****        - End of block with training records'
    print '*****        - Number of training records'
    print '*****          plus options'
    raise Exception()

  if (config.in_file_name == config.options[2]):
    print '***** Error: Input and output files must differ'
    print '*****        Input file name:          ', config.in_file_name
    print '*****        Output training file name:', config.options[2]
    raise Exception()

  first_rec = int(config.options[2])
  last_rec  = int(config.options[3])
  num_rec   = int(config.options[4])
  in_file_name = config.in_file_name
  out_file_name = config.options[1]

  # Check record number values  - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (int(first_rec) >= int(last_rec)) or \
     ((int(num_rec)-1) > (int(last_rec)-int(first_rec))):
    print '***** Error: Illegal values for training records block:'
    print '*****        - Start of block with training records:', first_rec
    print '*****        - End of block with training records:  ', last_rec
    print '*****        - Number of training records:          ', num_rec
    raise Exception()

  rec_range = last_rec-first_rec-1  # Range of records in input file

  # Open input file and check number of available records - - - - - - - - - - -
  #
  try:
    f_in = open(in_file_name,'r')
  except:
    inout.log_message('Cannot open input file: '+in_file_name,'err')
    raise IOError()

  line_count = 0
  for line in f_in.xreadlines():
    line_count += 1
  f_in.close()

  if (last_rec > line_count):  # Illegal value for last record
    print '***** Error: Illegal values for last training records:', last_rec
    print '*****        File only contains',line_count, 'lines/records'
    raise Exception()

  # Get tagging mode/lookup-tables used - - - - - - - - - - - - - - - - - - - -
  #
  tag_mode = config.options[0]
  if (tag_mode in ['name','na','n']):
    tag_mode = 'name'
  elif (tag_mode in ['locality','localty','loc','l']):
    tag_mode = 'loc'
  else:
    print '***** Error: Illegal tagging mode:', tag_mode
    print '*****        Must be either "name" or "locality"'
    raise Exception()

  # Check for optional arguments and process if any - - - - - - - - - - - - - -
  #
  config.verbose = 0     # Default: No verbose output
  config.logging = 0     # Default: No logging into a file
  hmm_file_name  = None  # Default: Do not use HMM to standardise training
                         #          records
  retag_file_name = None # Default: Do not retag an existing training file
  config.nowarn  = 0     # Deactivate no warning flag (print/log warning
                         # messages)
  freqs_file_name = None # Default: Do not write frequencies, no -freqs option

  if (len(config.options) > 5):
    options = config.options[5:]
    while (options != []):  # Do a loop processing all options

      if (options[0] == '-nowarn'):
        config.nowarn = 1  # Activate no warning flag
        options = options[1:]  # Remove processed '-nowarn' option

      elif (options[0] == '-v1'):
        config.verbose = 1  # Set to verbose output level 1
        options = options[1:]  # Remove processed '-v1' option

      elif (options[0] == '-v2'):
        config.verbose = 2  # Set to verbose output level 2
        options = options[1:]  # Remove processed '-v2' option

      elif (options[0] == '-l'):
        config.logging = 1
        if (len(options) > 1):
          if (options[1][0] != '-'):  # Not another option, must be a file name
            config.log_file = options[1]  # Get name of log file
            options = options[1:]  # Remove file_name
        options = options[1:]  # Remove processed -'l' option only

        try:
          f_log = open(config.log_file,'a')  # Test if file is appendable
        except:
          print '***** Error ********************',
          print '***** Cannot write to log file: '+config.log_file
          raise IOError()

        # Write (append) header to log file
        #
        f_log.write(os.linesep)
        f_log.write('##################################################')
        f_log.write('############'+os.linesep)
        f_log.write('#'+os.linesep)
        f_log.write("# 'pyTagData.py - Version 0.1' process started at: ")
        f_log.write(time.ctime(time.time())+os.linesep)
        f_log.write('#'+os.linesep)
        f_log.write("# Input file name:  "+in_file_name+os.linesep)
        f_log.write("# Output file name: "+out_file_name+os.linesep)
        f_log.write("# Tagging mode:     "+tag_mode+os.linesep)
        f_log.write(os.linesep)
        f_log.close()

      elif (options[0] == '-hmm'):
        hmm_file_name = options[1]  # Get file name of the HMM to use
        if (hmm_file_name == out_file_name):
          print '***** Error: HMM file name is the same as output file name!'
          raise Exception()

        try:
          f_in = open(hmm_file_name,'r')  # Test if file is available
        except:
          print '***** Error: Cannot open HMM file specified in "-hmm"',
          print 'option:', hmm_file_name
          raise IOError()
        f_in.close()
        options = options[2:]  # Remove processed '-hmm' option and file name

      elif (options[0] == '-retag'):
        if (hmm_file_name == None) and ('-hmm' not in options):
          print '***** Error: "-retag" option can only be used together with',
          print '"-hmm" option (which is not given).'
          raise Exception()

        retag_file_name = options[1]  # Get file name of the already-tagged
                                      # file to re-process
        if (retag_file_name == out_file_name):
          print '***** Error: Retag file name is the same as output file name!'
          raise Exception()
        elif (retag_file_name == in_file_name):
          print '***** Error: Retag file name is the same as input file name!'
          raise Exception()
        elif (retag_file_name == hmm_file_name):
          print '***** Error: Retag file name is the same as HMM file name!'
          raise Exception()

        try:
          f_in = open(retag_file_name,'r')  # Test if file is available

          # Now gather record numbers and previous tags/states, as well as the
          # original header information. Use a simple state machine to do this.
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

          num_rec = len(tagged_recs_keys)  # Override specified numbers
          first_rec = 0
          last_rec = line_count

        except:
          print '***** Error: Cannot open tagged training file specified',
          print 'in "-retag" option:', retag_file_name
          raise IOError()

        options = options[2:]  # Remove processed '-retag' option and file name

      elif (options[0][:5] == '-freq'):
        if (hmm_file_name == None) and ('-hmm' not in options):
          print '***** Error: "-feqs" option can only be used together with',
          print '"-hmm" option (which is not given).'
          raise Exception()

        freqs_file_name = options[1]  # File name to write the frequencies to
        if (freqs_file_name == out_file_name):
          print '***** Error: Frequency file name is the same as output',
          print 'file name!'
          raise Exception()
        elif (freqs_file_name == in_file_name):
          print '***** Error: Frequency file name is the same as input',
          print 'file name!'
          raise Exception()
        elif (freqs_file_name == hmm_file_name):
          print '***** Error: Frequency file name is the same as HMM',
          print 'file name!'
          raise Exception()

        options = options[2:]  # Remove processed '-freqs' option and file name
        try:  # Check if file writing is possible
          freqs_out = open(freqs_file_name,'w')
	  freqs_out.close()
        except:
          print '***** Error: Cannot write to frequency output file specified',
          print 'in "-freqs" option:', freqs_file_name
          raise IOError()

      else:
        print '***** Error: Illegal option:', options[0]
        raise Exception()

  # If specified initalise and load Hidden Markov Model (HMM) - - - - - - - - -
  #
  if (hmm_file_name != None):
    myhmm = simplehmm.hmm([],[])  # Create new empty HMM object
    myhmm.load_hmm(hmm_file_name)
    myhmm.print_hmm()  # Print HMM (according to verbose and logging level)

  # Open output file and write header - - - - - - - - - - - - - - - - - - - - -
  #
  try:
    f_out = open(out_file_name,'w')
  except:
    inout.log_message('Cannot open output file: '+out_file_name,'err')
    raise IOError()

  f_out.write("# Tagged training data written by 'pyTagData.py -"+ \
              " Version 0.1'"+os.linesep)
  f_out.write('#'+os.linesep)
  f_out.write('# Created '+time.ctime(time.time())+os.linesep)
  f_out.write('#'+os.linesep)
  f_out.write('# Input file name:  '+in_file_name+os.linesep)
  f_out.write('# Output file name: '+out_file_name+os.linesep)
  f_out.write('#'+os.linesep)
  f_out.write('# Parameters:'+os.linesep)
  f_out.write('# - Start of block with training records: '+str(first_rec)+ \
              os.linesep)
  f_out.write('# - End of block with training records:   '+str(last_rec)+ \
              os.linesep)
  f_out.write('# - Number of training records:           '+str(num_rec)+ \
              os.linesep)
  if (hmm_file_name != None):
    f_out.write('#'+os.linesep)
    f_out.write("# - Using HMM file '"+hmm_file_name+"' for standardisation"+ \
                os.linesep)
  if (retag_file_name != None):
    f_out.write('#'+os.linesep)
    f_out.write("# - Reprocessing training file '"+retag_file_name+"'"+ \
                os.linesep)
    f_out.write("#   Header lines from original training file follow:" + \
                os.linesep)
    for header_line in original_header_lines:
	    f_out.write(header_line + os.linesep)
  if (freqs_file_name != None):
    f_out.write('#'+os.linesep)
    f_out.write("# - Tag/state pattern frequencies written to file '" + \
                freqs_file_name + os.linesep)
  f_out.write('#'+'-'*70+os.linesep)
  f_out.write(os.linesep)

  rec_count    = 0        # Number of selected records
  num_rec_left = num_rec  # Number of records to be selected left
  rec_selected = {}       # Dictionary of all record numbers that were selected
  seq_freqs = {}          # Dict to hold examples of tag/state patterns

  unchanged_loop_cnt = 0       # Counter of how many loops have been done
                               # without new records being selected
  prev_num_rec_left = num_rec  # Number of records left in the previous
                               # interation

  # Due to the random nature of selecting records, and because sometimes  - - -
  # a selected component can be empty (and is thus not used for training)
  # more than one iteration over the input data set is carried out. In each 
  # iteration, records are selected randomly.
  #
  while (rec_count < num_rec):  # Loop until 'num_rec' records selected

    # Open input file - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      f_in = open(in_file_name,'r')
    except:
      inout.log_message('Cannot open input file: '+in_file_name,'err')
      raise IOError()

    line_read = 0  # Number of read lines

    # Skip to start of training block - - - - - - - - - - - - - - - - - - - - -
    #
    if (first_rec > 0):
      for i in range(first_rec):
        f_in.readline()

    while (rec_count < num_rec) and (line_read <= (last_rec-first_rec)):
      line = f_in.readline()

      if ((retag_file_name != None) and (line_read in tagged_recs_keys)) or \
         ((retag_file_name == None) and \
          (num_rec_left >= random.randrange(0,rec_range,1))):

        line = line.strip()  # Remove line separators
        config.curr_line = line  # Make a copy of the unprocessed current line

        line = line.lower()  # Make all characters lower case

        inout.log_message(['Record number: '+str(line_read+first_rec)],'v1')
        config.curr_line_no = line_read+first_rec  # Store current line number

        # Process line and extract content into components (name, geocode, etc)
        #
        [name_comp, geocode_comp, locality_comp, date1_comp, date2_comp] = \
           inout.process_line(line)

        # Select component and process it - - - - - - - - - - - - - - - - - - -
        #
        if (tag_mode == 'name'):
          if (type(name_comp) == types.ListType):
            component = name_comp[0].strip()+' '+name_comp[1].strip()
          else:
            component = name_comp.strip()
        else:  # Locality component
          component = geocode_comp.strip()+' '+locality_comp.strip()

        if (component != '') and \
           (not rec_selected.has_key((line_read+first_rec))):

          if (tag_mode == 'name'):
            inout.log_message('  Name component: |'+component+'|','v1')

            component = name.clean_name_component(component)
            [word_list, tag_list] = name.tag_name_component(component)
 
          else:  # Locality component
            inout.log_message('  Locality component: |'+component+'|','v1')

            component = locality.clean_geoloc_component(component)
            [word_list, tag_list] = locality.tag_geoloc_component(component)

          if (tag_list != []):  # Only process non-empty tag lists

            # Append record number into dictionary of processed records
            #
            rec_selected.update({(line_read+first_rec):(line_read+first_rec)})

            # Create all permutation sequences of this tag list - - - - - - - -
            #
            tag_seq = mymath.perm_tag_sequence(tag_list)

            inout.log_message(['  Word list: '+str(word_list), \
                               '  Tag list: '+str(tag_list), \
                               '  Tag sequences:'],'v2')

            # Do HMM processing - - - - - - - - - - - - - - - - - - - - - - - -
            #
            if (hmm_file_name != None):

              state_seq  = []    # List containing computed HMM state sequences
              max_prob   = -1.0  # maximal probability for a sequence
              max_seq_no = -1    # Number of the seq. with the max. probablity

              # Now give tag sequences to the HMMs to compute state sequences
              #
              i = 0
              for t in tag_seq:
                [obs_seq, prob] = myhmm.viterbi(t)
                state_seq.append(obs_seq)
                if (prob > max_prob):
                  max_prob = prob
                  max_seq_no = i
                i += 1

            # Write original component and resulting tag sequences to output
            #
            f_out.write('# '+str(line_read+first_rec)+' ('+str(rec_count)+ \
                        '): |'+component+'|'+os.linesep) # Commented original
            num_len = len(str(line_read+first_rec))+len(str(rec_count))+6

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
              inout.log_message('    '+seq_string[:-1],'v2')

            if (hmm_file_name != None):
              f_out.write('# Maximum Viterbi probability: %0.5f'% \
                          (max_prob) + os.linesep)
              inout.log_message('Maximum Viterbi probability: %0.5f'% \
                                (max_prob), 'v2')

            if (retag_file_name != None) and (tagged_recs[line_read] != None):
              if (tagged_recs[line_read].strip() != seq_string[:-1].strip()):
                f_out.write("# Note: ***** Changed *****" + os.linesep)
                inout.log_message('                      Note:' + \
                                  ' ***** Changed *****','v2')
                f_out.write('# Was: ' + tagged_recs[line_read]+os.linesep)
                            # Write commented original tag sequence
                inout.log_message('Original tag sequence: '+ \
                                  tagged_recs[line_read],'v2')

            f_out.write(os.linesep)  # Write an empty line
            inout.log_message('','v1')  # Print empty lines between records

            if (hmm_file_name != None):
              seq_key = seq_string[:-1]  # Add sequence to dictionary
              if (seq_freqs.has_key(seq_key)):
                seq_freqs[seq_key].append(['|'+' '.join(word_list)+'|', \
                                          max_prob])
              else:
                seq_freqs[seq_key] = [['|'+' '.join(word_list)+'|', \
                                      max_prob]]

            rec_count += 1

            # Print process indicator message
            #
            if (config.proc_ind >= 0) and (rec_count > 0):
              if (rec_count % config.proc_ind == 0):
                print 'Processed line', rec_count, 'of', num_rec

      line_read += 1

    f_in.close()

    num_rec_left = num_rec - rec_count

    if (prev_num_rec_left == num_rec_left):  # No new records selected
      unchanged_loop_cnt += 1
    prev_num_rec_left = num_rec_left  # Set to current value

    if (unchanged_loop_cnt > 5):  # Do five loops maximal without selecting
                                  # new records
      config.curr_line_no = -1  # Set to illegal/empty values, as warning is
      config.curr_line    = ''  # not related to the current input line
      inout.log_message(['Can not select more than '+str(rec_count)+ \
                         ' records for training.', \
                         'This is probably due to empty input components.', \
                         'Please reduce value of "num_rec" or increase ' + \
                         'range','between "first_rec" and "last_rec".'],'warn')
      break

    if (num_rec_left < 10):  # Only 10 records left to select
      num_rec_left = num_rec+1  # Set to more than 100% probablity
    elif (num_rec_left < (num_rec / 100.0)):  # Less than 1% records left
      num_rec_left = int(num_rec / 100.0)  # Set to 1%

  f_out.close()

  # If specified, save Viterbi frequencies to a file  - - - - - - - - - - - - -
  #
  if (freqs_file_name != None):
    freqs_out = open(freqs_file_name,'w')  # Open frequency file for writing
    freqs_out.write('# Frequency listing of tag/state patterns written by')
    freqs_out.write('"pyTagData.py - Version 0.1"'+os.linesep)
    freqs_out.write('#'+os.linesep)
    freqs_out.write('# Created '+time.ctime(time.time())+os.linesep)
    freqs_out.write('#'+os.linesep)
    freqs_out.write("# Input file name:  "+in_file_name+os.linesep)
    freqs_out.write("# Output file name: "+out_file_name+os.linesep)
    freqs_out.write(os.linesep)
    freqs_out.write('# Parameters:'+os.linesep)
    freqs_out.write('# - Start of block with training records: '+ \
                    str(first_rec)+os.linesep)
    freqs_out.write('# - End of block with training records:   '+ \
                    str(last_rec)+os.linesep)
    freqs_out.write('# - Number of training records:           '+ \
                    str(num_rec)+os.linesep)
    if (hmm_file_name != None):
      freqs_out.write('#'+os.linesep)
      freqs_out.write("# - Using HMM file '"+hmm_file_name+ \
                      "' for standardisation"+os.linesep)
    if (retag_file_name != None):
      freqs_out.write('#'+os.linesep)
      freqs_out.write("# - Reprocessing training file '"+retag_file_name+ \
                      "'"+os.linesep)
    freqs_out.write('#'+'-'*70+os.linesep)
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

  inout.log_message(['Read '+str(line_read)+' lines, processed '+ \
                    str(rec_count)+' lines', 'End.'],'v1')

# ----------------------------------------------------------------------------

# Start main data tagging routine
#
tagdata()

# ----------------------------------------------------------------------------
