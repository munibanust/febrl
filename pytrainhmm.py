# =============================================================================
# pyTrainHMM.py - Main module to train a Hidden Markov Model (HMM)
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
# The Original Software is "pyTrainHMM.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module pyTrainHMM.py - Main module to train a Hidden Markov Model (HMM)

   USAGE:
     pyTrainHMM.py  [project_module] [tag_mode] [train_file] [hmm_file]
                    [options]

   ARGUMENTS:
     project_module  A Python module containing all project settings.
     tag_mode        Mode for tagging, must be either 'name' or 'locality'
     train_file      Name of the input data file containing the tagged
                     training records.
     hmm_file        The name of the HMM file to be written. It's a text file
                     containing the state and observation (tag) names; initial,
                     transition and observation probabilities.

   OPTIONS (after arguments);
     -s [method]     Smoothing of HMM parameters. The method can either be set
                     to 'laplace' (for Laplace smoothing) or 'absdiscount'
                     (for the absolute discount smoothing method).
     -l [file_name]  Log warnings and errors into the given file
     -v1             Verbose output level 1 (low output)
     -v2             Verbose output level 2 (high output)
     -nowarn         Don't print warning messages

   DESCRIPTION:
     This module can be used to train a Hidden Markov Model (HMM) using tagged
     data as created by 'pyTagData.py' and then edited manually. The 'tag_mode'
     argument (either 'name' or 'locality') is used to get the possible HMM
     states and observations from the 'config.py' module. They are also given
     in the file: './hmm/hmm-states.txt'

     The format of the input file is as follows:
     - Comment lines must start with a hash character (#).
     - Empty lines are possible.
     - Each non-empty line that is not a comment line must contain one training
       record.
     - Training records must contain a comma separated list of pairs
       'observation:state' (i.e. tag:state). Only known states and observations
       (as listed in 'config.py') are allowed, unknown states or observations
       will result in an error.

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

     The output file is a text file with all parameters of the HMM
     (see 'simplehmm.py' module for more details).

     For more information on the smoothing methods see the 'simplehmm.py'
     module ('train' routine) or e.g.
       V.Borkar et.al., Automatic Segmentation of Text into Structured Records
       Section 2.2

   EXAMPLE:
     The following lines are an example of the input file format:

       # 0: |dr peter baxter dea|
       # TI:, GM:, GM:, GF:
       # TI:, GM:, SN:, GF:
         TI:title, GM:givenname, GM:givenname, SN:surname
       # TI:, GM:, SN:, SN:

       # 1: |miss monica mitchell meyer|
       # TI:, UN:, GM:, SN:
         TI:title, UN:givenname, SN:surname, SN:surname

       # 2: |phd tim william jones harris|
         TI:title, GM:givenname, GM:altgivenname, UN:surname, SN:altsurname
"""

# -----------------------------------------------------------------------------

import sys
import os
import xreadlines
import time

import config
import inout
import simplehmm

# -----------------------------------------------------------------------------

def trainhmm():
  """Main routine, open file, read lines, train HMM and save it to file.

  USAGE:
    trainhmm()

  ARGUMENTS:
    None

  DESCRIPTION:
    Main routine, see description of module above.
  """

  # Process command line arguments and check for correctness  - - - - - - - - -
  #
  if (len(config.options) < 3):
    print '***** Error: %s needs at least four arguments:'% (sys.argv[0])
    print '*****        - Name of the project module'
    print '*****        - Tagging mode: "name" or "locality"'
    print '*****        - Input training file name'
    print '*****        - HMM output file name'
    print '*****          plus options'
    raise Exception()

  if (config.options[1] == config.options[2]):
    print '*** Error: Input and output files must differ'
    print '***        Input training file name:', config.options[1]
    print '***        HMM output file name:    ', config.options[1]
    raise Exception()

  in_file_name  = config.options[1]
  hmm_file_name = config.options[2]

  # Get tagging mode/lookup-tables used - - - - - - - - - - - - - - - - - - - -
  #
  tag_mode = config.options[0]
  if (tag_mode in ['name','na','n']):
    tag_mode = 'name'
  elif (tag_mode in ['locality','lolty','loc','l']):
    tag_mode = 'loc'
  else:
    print '***** Error: Illegal tagging mode:', tag_mode
    print '*****        Must be either "name" or "locality"'
    raise Exception()

  # Check for optional arguments and process if any - - - - - - - - - - - - - -
  #
  config.verbose = 0     # Default: No verbose output
  config.logging = 0     # Default: No logging into a file
  smoothing      = None  # Default: No smoothing
  config.nowarn  = 0     # Deactivate no warning flag (print/log warning
                         # messages)

  if (len(config.options) > 3):
    options =  config.options[3:]
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
        f_log.write("############"+os.linesep)
        f_log.write("#"+os.linesep)
        f_log.write("# 'pyTrainHMM.py - Version 0.1' process started at: ")
        f_log.write(time.ctime(time.time())+os.linesep)
        f_log.write("#"+os.linesep)
        f_log.write("# Input file name: "+in_file_name+os.linesep)
        f_log.write("# HMM file name:   "+hmm_file_name+os.linesep)
        f_log.write(os.linesep)
        f_log.close()

      elif (options[0] == '-s'):
        smoothing = 1  # Set to do a HMM smoothing
        smoothing = options[1]
        if (smoothing in ['l','la','lap','laplac','laplace']):
          smoothing = 'laplace'
        elif (smoothing in ['a','ad','abs','absd','absdis','absdisc',\
               'absdiscount']):
          smoothing = 'absdiscount'
        else:  # Illegal value
          print "*** Error: Illegal value for 'smoothing' argument:", smoothing
          print "***        Possible are: 'laplace' or 'absdiscount'"
          raise Exception()

        options = options[2:]  # Remove processed option

      else:
        print '*** Error: Illegal option:', options[0]
        raise Exception()

  # Get HMM states and observations from configuration module - - - - - - - - -
  #
  if (tag_mode == 'name'): 
    state_list = config.name_hmm_states
    obser_list = config.name_hmm_obser

  else:
    state_list = config.geoloc_hmm_states
    obser_list = config.geoloc_hmm_obser

  # Open input file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  try:
    f_in = open(in_file_name,'r')
  except:
    inout.log_message('Cannot open input file: '+in_file_name,'err')
    raise IOError()

  line_count = 0  # Counter for lines read
  rec_count  = 0  # Counter for training records read

  # Read lines, discard comment lines and process training data lines - - - - -
  #
  training_data = []  # List of training records

  train_list = []  # List of training sequences (dictionaries), extracted from
                   # training data

  for line in xreadlines.xreadlines(f_in):

    if (line[0] != '#') and (line.strip() != ''):
      # Line must contain a training record

      line = line.strip()  # Remove line separators
      config.curr_line = line  # Make a copy of the unprocessed current line

      line_list = line.split(',')  # Split into a list of elements
      line_data = []  # Training data list for one training record

      inout.log_message(['Record number: '+str(rec_count)],'v1')
      config.curr_line_no = line_count  # Store current line number

      for elem in line_list:
        [k,v] = elem.split(':')  # Split into key and value
        tag = k.strip()
        state = v.strip()
        line_data.append((state,tag))

        if (state not in state_list):
          msg = ['Illegal state name in training record: '+state, \
                 'Line: '+str(line_count)+', record: '+str(rec_count), \
                 'Possible values: '+str(state_list)]
          inout.log_message(msg,'err')
          raise Exception()

        if (tag not in obser_list):
          msg = ['Illegal observation (tag) name in training record: '+tag, \
                 'Line: '+str(line_count)+', record: '+str(rec_count), \
                 'Possible values: '+str(obser_list)]
          inout.log_message(msg,'err')
          raise Exception()

      inout.log_message('  Training record '+str(rec_count)+':'+ \
                        str(line_data),'v1')

      train_list.append(line_data)

      rec_count += 1
      inout.log_message('','v1')  # Print empty lines between records

    line_count += 1

  # Close input file  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  f_in.close()

  inout.log_message('','v1')  # Print empty lines between records

  # Initalise HMM and train it with training data - - - - - - - - - - - - - - -
  #
  myhmm = simplehmm.hmm(state_list, obser_list)

  myhmm.train(train_list,smoothing)
  myhmm.print_hmm()

  # Save trained HMM  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  myhmm.save_hmm(hmm_file_name)  

  inout.log_message(['Read '+str(line_count)+' lines, processed '+ \
                    str(rec_count)+' training records', 'End.'],'v1')

# ----------------------------------------------------------------------------

# Start main data HMM training routine
#
trainhmm()

# ----------------------------------------------------------------------------
