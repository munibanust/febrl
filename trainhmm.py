# =============================================================================
# trainhmm.py - Main module to train a Hidden Markov Model (HMM)
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
# The Original Software is "trainhmm.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module trainhmm.py - Main module to train a Hidden Markov Model (HMM)

   DESCRIPTION:
     This module can be used to train a Hidden Markov Model (HMM) using tagged
     data as created by 'tagdata.py' and then edited manually.

     The 'hmm_component' can be either 'name' or 'address', and is used to 
     determine the possible HMM states and observations to be used.

     The format of the input training file is as follows:
     - Comment lines must start with a hash character (#).
     - Empty lines are possible.
     - Each non-empty line that is not a comment line must contain one training
       record.
     - Training records must contain a comma separated list of pairs
       'observation:state' (i.e. tag:state). Only known states and observations
       (as listed in 'config.py') are allowed, unknown states or observations
       will result in an error.

     The output file is a text file with all parameters of the HMM
     (see 'simplehmm.py' module for more details).

     For more information on the smoothing methods see the 'simplehmm.py'
     module ('train' routine) or e.g.
       V.Borkar et.al., Automatic Segmentation of Text into Structured Records
       Section 2.2

   EXAMPLE:
     The following lines are an example of the input training file format:

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

# =============================================================================
# Imports go here

from febrl import *      # Main Febrl classes
from dataset import *    # Data set routines
from simplehmm import *  # Hidden Markov model (HMM) routines

import xreadlines  # Python standard module for efficient text file reading
import time        # Python standard module for time functions

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

hmm_febrl = Febrl(description = 'HMM training Febrl instance',
                   febrl_path = '.')

hmm_project = hmm_febrl.new_project(name = 'HMM-Train',
                             description = 'Training module for HMMs',
                               file_name = 'hmm.fbr')

# =============================================================================
# Define a project logger

hmm_log = ProjectLog(file_name = 'hmm_train.log',
                       project = hmm_project,
                     log_level = 1,
                 verbose_level = 1,
                     clear_log = True,
                       no_warn = False)

# =============================================================================
# Define settings for HMM training

# Name of the file containing training records  - - - - - - - - - - - - - - - -
#
hmm_train_file = 'tagged_data.csv'  # './hmm/name-train.csv'

# Name of the HMM file to be written  - - - - - - - - - - - - - - - - - - - - -
#
hmm_model_file = 'test.hmm'

# Name of the HMM - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
hmm_name = 'Test Name HMM'

# Component: Can either be 'name' or 'address'  - - - - - - - - - - - - - - - -
#
hmm_component = 'name'

# HMM smoothing method, can be either None, 'laplace' or 'absdiscount'  - - - -
#
hmm_smoothing = 'absdiscount'

# =============================================================================
# =============================================================================
# Do not change anything below here
# =============================================================================
# =============================================================================

# State and observation (tag) lists for names and addresses - - - - - - - - - -
#
name_states = ['titl','baby','knwn','andor','gname1','gname2','ghyph',
               'gopbr','gclbr','agname1','agname2','coma','sname1','sname2',
               'shyph','sopbr','sclbr','asname1','asname2','pref1','pref2',
               'rubb']
name_tags = ['NU','AN','TI','PR','GF','GM','SN','ST','SP','HY','CO','NE','II',
             'BO','VB','UN','RU']

address_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty','prna1',
                  'prna2','inna1','inna2','inty','panu','paty','hyph','sla',
                  'coma','opbr','clbr','loc1','loc2','locql','pc','ter1',
                  'ter2','cntr1','cntr2','rubb']
address_tags = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT','LQ','WT',
                'WN','UT','HY','SL','CO','VB','PA','UN','RU']

# Test settings - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
if (not isinstance(hmm_train_file, str)) or (hmm_train_file == ''):
  print 'error:HMM training file is not a valid string: "%s"' % \
        (str(hmm_train_file))
  raise Exception

if (not isinstance(hmm_model_file, str)) or (hmm_model_file == ''):
  print 'error:HMM model file is not a valid string: "%s"' % \
        (str(hmm_model_file))
  raise Exception

if (hmm_train_file == hmm_model_file):
  print 'error:HMM training file and HMM model file must differ'
  raise Exception

if (hmm_smoothing not in [None, 'laplace','absdiscount']):
  print 'error:Illegal HMM smoothing method: %s' % (str(hmm_smoothing))
  raise Exception

if (hmm_component == 'name'):  # Set HMM states and observations for names
  hmm_states = name_states
  hmm_observ = name_tags
elif (hmm_component == 'address'):  # Set HMM states and observations for
  hmm_states = address_states       # addresses
  hmm_observ = address_tags
else:
  print 'error:Illegal HMM component: %s' % (str(hmm_component))
  raise Exception

# Print header (and write it to the log file if activated)  - - - - - - - - - -
#
print '1:'
print '1:'+'#'*75
print '1:#'
print '1:# "hmmtrain.py" - Version 0.2'
print '1:# Process started at: %s' % (time.ctime(time.time()))
print '1:#'
print '1:# Training file:        %s' % (str(hmm_train_file))
print '1:# HMM model file:       %s' % (str(hmm_model_file))
print '1:# HMM component:        %s' % (str(hmm_component))
if (hmm_smoothing != None):
  print '1:# HMM smoothing method: %s' % (str(hmm_smoothing))
print '1:#'

# Open training file  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
try:
  f_in = open(hmm_train_file, 'r')
except:
  print 'error:Cannot open training file: "%s"' % (str(hmm_train_file))
  raise IOError

rec_count  = 0  # Counter for training records read
line_count = 0  # Counter for lines in training file

# Read lines, discard comment lines and process training data lines - - - - - -
#
train_records = []  # List of training records

for line in xreadlines.xreadlines(f_in):

  if (line[0] != '#') and (line.strip() != ''): # Skip comment or empty linea

    line = line.strip()  # Remove whitespaces and line separators

    line_list = line.split(',')  # Split into a list of elements
    line_data = []               # Training data list for one training record

    for elem in line_list:
      [key, values] = elem.split(':')  # Split into key and value
      tag = key.strip()
      state = values.strip()

      if (state not in hmm_states):
        print 'error:Illegal state name "%s" in record ' % (str(state)) + \
              '%s in training file line %i' % (str(rec_count), line_count) + \
              ', possible values: %s' % (str(hmm_states))
        raise Exception

      if (tag not in hmm_observ):
        print 'error:Illegal observation (tag) name "%s"' % (str(state)) + \
              ' in record %i in training file line ' % (rec_count) + \
              '%i, possible values: %s' % (line_count, str(hmm_observ))
        raise Exception

      if (state in hmm_states) and (tag in hmm_observ):
        line_data.append((state, tag))  # Append (state,tag) tuple

    print '3:  Training record %s: %s' % (str(rec_count), str(line_data))

    if (line_data != []):  # Only append training records which are valid
      train_records.append(line_data)
      rec_count += 1

  line_count += 1

# Close input file  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
f_in.close()

# Check if there are training records avaiable  - - - - - - - - - - - - - - - -
#
if (train_records == []):
  print 'error:No training records extracted from training file'
  raise Exception

# Initalise HMM and train it with training data - - - - - - - - - - - - - - - -
#
train_hmm = hmm(hmm_name, hmm_states, hmm_observ)

train_hmm.train(train_records, hmm_smoothing)
train_hmm.print_hmm()

# Save trained HMM  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
train_hmm.save_hmm(hmm_model_file)  

print '1:'
print '1:End.'

# =============================================================================
