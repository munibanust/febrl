# =============================================================================
# febrl.py - Main module and classes for febrl projects.
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
# The Original Software is "febrl.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module febrl.py - Main module and classes for febrl projects.

   TODO
   - Allow saving and loading not only of pickled project files, but also XML,
     text, compressed and uncompressed, etc.
     Then, do not restrict project files to have a '.fbl' extension, but list
     all files found in a project directory (incl. their tyes)
     Changes in Project.save() and Febrl.__init__() needed.
"""

# =============================================================================
# The following flags can be set to True or False
# They are used for testing the parallel functionalities of Febrl and should
# be set to False for normal use.

DO_PARALLEL_TEST =         False  # Perform several parallel tests of
                                  # intermediate results
SAVE_PARALLEL_TEST_FILES = False  # Write intermediate results to file for
                                  # inspection

# =============================================================================

import cPickle
import os
import sys
import time
import traceback
import types
import copy

import parallel
import indexing
import output
import lap

# =============================================================================

class Febrl:
  """Class Febrl - Main class for Febrl projects.
  """

  def __init__(self, **kwargs):
    """Constructor - Set attributes and load list of available projects.
    """

    self.version_major =      '0.2.1'
    self.version_minor =      ''
    self.version =            self.version_major+'.'+self.version_minor
    self.license =            'ANUOS Version 1.1'
    self.copyright =          '(C) 2002, 2003 the Australian National ' + \
                              'University and others'
    self.initial_developers = 'Dr Peter Christen (Department of Computer ' + \
                              'Science, Australian National University) ' + \
                              'and Dr Tim Churches (Centre for Epidemiology'+ \
                              ' and Research, New South Wales Department ' + \
                              'of Health)'
    self.contributors =       ''

    self.description = None
    self.febrl_path =  os.curdir  # Default set to current directory
    self.project_names = []

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'febrl_path'):
        self.febrl_path = value
      elif (keyword == 'description'):
        self.description = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if Febrl projects are available in the 'project_path' directory by
    # scanning the directory for files with '.fbr' file extension
    #
    file_list = os.listdir(self.febrl_path)
    for fn in file_list:
      file_name = fn.strip().lower()
      if (file_name[-4:] == '.fbr'):
        self.project_names.append(file_name)

  # ---------------------------------------------------------------------------

  def __str__(self):
    """Create a string representation of the Febrl object.
    """

    linesep = os.linesep

    rep = 'Febrl (Freely extensible biomedical record linkage)' + linesep
    rep += '---------------------------------------------------' + linesep
    rep += linesep
    rep += '  Version: ' + self.version + linesep
    rep += '  License: ' + self.license + linesep
    rep += '  Copyright: ' + self.copyright + linesep + linesep

    rep += '  Initial developers: ' + self.initial_developers + linesep
    rep += '  Contributors: ' + self.contributors + linesep + linesep

    rep += 'Description: ' + str(self.description) + linesep
    rep += 'Febrl path: ' + str(self.febrl_path) + linesep
    rep += 'Avaliable projects:' + linesep
    if (len(self.project_names) == 0):
      rep += '  ' + str(None) + linesep
    else:
      i = 0
      for pn in self.project_names:
        rep += str(i).rjust(3) + ': ' + pn + linesep
        i += 1
      rep += linesep

    return rep                      

  # ---------------------------------------------------------------------------

  def load_project(self, project, project_path=None):
    """Load a project from file.

       A project can either be the file name as a string or the project
       number (as stored in the list of project names).
    """

    # Check project path  - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (project_path == None):  # No project path given, use default
      file_name = self.febrl_path
    else:
      file_name = project_path

    if (file_name[-1] != os.sep):
      file_name += os.sep  # Make sure path ends with a directory separator

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (type(project) == types.IntType):
      try:
        file_name += self.project_names[project]  # Get project file name
      except:
        print 'error:Illegal project number: %s' % (str(project))
        raise Exception

    elif (type(project) == types.StringType):
      file_name += project
    else:
      print 'error:Illegal type for "project" argument, must be either of ' + \
            'type string or integer'
      raise Exception

    # Open project file and load it - - - - - - - - - - - - - - - - - - - - - -
    #
    f = open(file_name,'r')
    loaded_project = cPickle.loads(f.read())
    f.close()
    loaded_project.febrl = self

    return loaded_project

  # ---------------------------------------------------------------------------

  def new_project(self, **kwargs):
    """Create a new project object and populate it.
    """

    new_project = Project(self, **kwargs)

    return new_project

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a Febrl project.
    """

    print '1:'
    print '1:Febrl stopped.'

    parallel.Barrier()  # Make sure all processes are here

    parallel.Finalize()  # Finalise parallel environment

    sys.exit()

# =============================================================================

class Project:
  """Class for record linkage projects.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, febrl, **kwargs):
    """Constructor - Create a new project object.
    """

    self.febrl =            febrl
    self.name =             ''
    self.description =      ''
    self.file_name =        None
    self.project_path =     febrl.febrl_path  # Inherit path
    self.block_size =       10000  # File blocking size (in number of records)
    self.parallel_write =   'host' # Set either to 'host' (default) or 'all'

    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'file_name'):
        self.file_name = value
      elif (keyword == 'project_path'):
        self.project_path = value

      elif (keyword == 'block_size'):
        if (not isinstance(value, int)) and (value > 0):
          print 'error:Argument "block_size" is not a positive integer'
          raise Esception
        self.block_size = value

      elif (keyword == 'parallel_write'):
        if (value not in ['host','all']):
          print 'error:Argument "parallel_write" must be set to "host"'+ \
                ' or "all"'
          raise Exception
        self.parallel_write = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Set the parallel writing/saving mode  - - - - - - - - - - - - - - - - - -
    #
    parallel.writemode = self.parallel_write

  # ---------------------------------------------------------------------------

  def __str__(self):
    """Create a string representation for this project.
    """

    linesep = os.linesep

    rep = linesep + 'Febrl project: "'+self.name + '"' + linesep
    rep += '  Description:  ' + self.description + linesep
    rep += '  Filename:     ' + self.file_name + linesep
    rep += '  Project path: ' + self.project_path + linesep

    return rep

  # ---------------------------------------------------------------------------

  def save(self, path=None):
    """Save the project into a file (currently a pickled file).
    """

    if (path is None):
      path = self.project_path  # Take the project's path

    if (path[-1] != os.sep):
      path += os.sep  # Make sure path ends with a directory separator

    # Unset the current febrl object
    #
    save_febrl = copy.copy(self.febrl)  # Make a deep copy first
    self.febrl = None

    file_name = path + self.file_name
    f = open(file_name, 'w+')
    f.write(cPickle.dumps(self, 1))
    f.close()

    # Restore febrl object
    #
    self.febrl = save_febrl

  # ---------------------------------------------------------------------------

  def standardise(self, **kwargs):
    """Clean and standardise the given data set using the defined record
       standardiser.

       Records are loaded block wise from the input data set, then standardised
       and written into the output data set.

       If the argument 'first_record' is not given, it will automatically be
       set to the first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it will be set
       to the total number of records in the input data set.

       The output data set can be any data set type except a memory based data
       set (as all standardised records would lost once the program finishes).
       This output data set has to be initialised in 'write', 'append' or
       'readwrite' access mode.
    """

    self.input_dataset =    None   # A reference ot the (raw) input data set
                                   # data set
    self.output_dataset =   None   # A reference to the output data set

    self.rec_standardiser = None   # Reference to a record standardiser
    self.first_record =     None   # Number of the first record to process
    self.number_records =   None   # Number of records to process

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'output_dataset'):
        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value

      elif (keyword == 'first_record'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record" is not a valid integer number'
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records" is not a positive integer '+ \
                'number'
          raise Exception
        self.number_records = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      print 'error:Input data set is not defined'
      raise Exception

    if (self.output_dataset == None):
      print 'error:Output data set is not defined'
      raise Exception
    elif (self.output_dataset.dataset_type == 'MEMORY'):
      print 'error:Output data set can not be a memory based data set'
      raise Exception
    if (self.output_dataset.access_mode not in ['write','append','readwrite']):
      print 'error:Output dataset must be initialised in one of the access' + \
            ' modes: "write", "append", or "readwrite"'
      raise Exception

    if (self.rec_standardiser == None):
      print 'error:Record standardiser is not defined'
      raise Exception

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):  # Process all records
      self.number_records = self.input_dataset.num_records

    print '1:'
    print '1:***** Standardise data set: "%s" (type "%s")' % \
          (self.input_dataset.name, self.input_dataset.dataset_type)
    print '1:*****        into data set: "%s" (type "%s")' % \
          (self.output_dataset.name, self.output_dataset.dataset_type)
    print '1:'

    # Call the main standardisation routine - - - - - - - - - - - - - - - - - -
    # (no indexing will be done)
    #
    [stand_time, comm_time] = do_load_standard_indexing(self.input_dataset,
                                                        self.output_dataset,
                                                        self.rec_standardiser,
                                                        None,
                                                        self.first_record,
                                                        self.number_records,
                                                        self.block_size)

  # ---------------------------------------------------------------------------

  def deduplicate(self, **kwargs):
    """Deduplicate the given data set using the defined record standardiser,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data set, then standardised
       (if the record standardiser is defined, otherwise the input data set is
       directly deduplicated), linked and the results are printed and/or saved
       into the result file(s).

       If the argument 'first_record' is not given, it will automatically be
       set to the first record in the data set (i.e. record number 0).
       Similarly, if the argument 'number_records' is not given, it will be set
       to the total number of records in the input data set.

       The temporary data set must be a random acces data set implementation,
       i.e. either a Shelve or a Memory data set. For large data set it is
       recommended to use a Shelve data set. This temporary data set has to be
       initialised in access mode 'readwrite'.

       Currently, the output can be a printed or saved list of record pairs in
       both a detailed and condensed form (if the arguments
       'output_rec_pair_details' and 'output_rec_pair_weights' are set to
       'True' or to a file name (a string). The output can be filtered by
       setting the 'output_threshold' (meaning all record pairs with a weight
       less then this threshold are not printed or saved).

       In future versions, it will be possible to compile an output data set.

       A histogram can be saved or printed by setting the argument
       'output_histogram' to 'True' or to a file name.

       It is also possible to apply a one-to-one assignment procedure by
       setting the argument 'output_assignment' to 'one2one'.
    """

    self.input_dataset =    None   # A reference ot the (raw) input data set
    self.tmp_dataset =      None   # A reference to a temporary (random access)
                                   # data set
#    self.output_dataset =   None   # A reference to the output data set

    self.rec_standardiser = None   # Reference to a record standardiser
    self.rec_comparator =   None   # Reference to a record comparator
    self.blocking_index =   None   # Reference to a blocking index
    self.classifier =       None   # Reference to a weight vector classifier

    self.first_record =     None   # Number of the first record to process
    self.number_records =   None   # Number of records to process

    self.output_histogram = False         # Set to True, a file name or False
                                          # (default) if a histogram of weights
                                          # should be printed or saved
    self.output_rec_pair_details = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved in details
    self.output_rec_pair_weights = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved with weights
    self.output_threshold = None          # Set to a weight threshold (only
                                          # record pairs with weights equal to
                                          # or above will be saved and or
                                          # printed)
    self.output_assignment = None         # Set to 'one2one' if one-to-one
                                          # assignment should be forced
                                          # (default: None)

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'tmp_dataset'):
        self.tmp_dataset = value
#      elif (keyword == 'output_dataset'):
#        self.output_dataset = value

      elif (keyword == 'rec_standardiser'):
        self.rec_standardiser = value
      elif (keyword == 'rec_comparator'):
        self.rec_comparator = value
      elif (keyword == 'blocking_index'):
        self.blocking_index = value
      elif (keyword == 'classifier'):
        self.classifier = value

      elif (keyword == 'first_record'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record" is not a valid integer number'
          raise Exception
        self.first_record = value
      elif (keyword == 'number_records'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records" is not a positive integer '+ \
                'number'
          raise Exception
        self.number_records = value

      elif (keyword == 'output_rec_pair_details'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_rec_pair_details" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_rec_pair_details = value
      elif (keyword == 'output_rec_pair_weights'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_rec_pair_weights" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_rec_pair_weights = value
      elif (keyword == 'output_histogram'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_histogram" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "output_threshold" is not a number: %s' % \
                (str(value))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          print 'error:Illegal value for argument "output_assignment": %s' % \
                (str(value))
          raise Exception
        else:
          self.output_assignment = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Do some checks on the input arguments - - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      print 'error:Input data set is not defined'
      raise Exception

    if (self.tmp_dataset == None):
      print 'error:Temporary data set is not defined'
      raise Exception
    elif (self.tmp_dataset.dataset_type not in ['SHELVE', 'MEMORY']):
      print 'error:Temporary data set must be a random access data set' + \
            ' (either Shelve or Memory)'
      raise Exception
    if (self.tmp_dataset.access_mode not in ['write','append','readwrite']):
      print 'error:Temporary data set must be initialised in one of the ' + \
            'access  modes: "write", "append", or "readwrite"'
      raise Exception

#    if (self.output_dataset == None):
#      print 'error:Output data set is not defined'
#      raise Exception

    # Make sure at least one output is defined
    #
    if (self.output_rec_pair_weights == False) and \
       (self.output_rec_pair_details == False) and \
       (self.output_histogram == False):
      print 'error:No ouput of results is defined.'
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record == None):
      self.first_record = 0  # Take default first record in data set

    if (self.number_records == None):  # Process all records
      self.number_records = self.input_dataset.num_records

    if (self.rec_comparator == None):
      print 'error:No record comparator defined'
      raise Exception
    if (self.rec_comparator.dataset_a != self.tmp_dataset) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset):
      print 'error:Illegal data set definition in record comparator'
      raise Exception

    if (self.blocking_index == None):
      print 'error:No blocking index defined'
      raise Exception

    if (self.classifier == None):
      print 'error:No classifier defined'
      raise Exception
    if (self.classifier.dataset_a != self.tmp_dataset) or \
       (self.classifier.dataset_b != self.tmp_dataset):
      print 'error:Illegal data set definition in classifier'
      raise Exception

    total_time = time.time()  # Get current time

    print '1:'
    print '1:***** Deduplicate data set: "%s" (type "%s")' % \
          (self.input_dataset.name, self.input_dataset.dataset_type)
    print '1:*****   Temporary data set: "%s" (type "%s")' % \
          (self.tmp_dataset.name, self.tmp_dataset.dataset_type)
    print '1:'
    print '1:Step 1: Loading, standardisation and indexing'
    print '1:-------'
    print '1:'

    step_1_time = time.time()  # Get current time

    # Call the main standardisation routine - - - - - - - - - - - - - - - - - -
    #
    [p, step_1_comm_time] = do_load_standard_indexing(self.input_dataset,
                                                      self.tmp_dataset,
                                                      self.rec_standardiser,
                                                      self.blocking_index,
                                                      self.first_record,
                                                      self.number_records,
                                                      self.block_size)

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index from process %i' % (p)

        self.blocking_index.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:    Sent index to process 0'

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent index to process %i' % (p)

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index from process 0'

        self.blocking_index.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index.compact()

    step_1_time = time.time() - step_1_time  # Calculate time for step 1
    step_1_time_string = output.time_string(step_1_time)

    print '1:'
    print '1:Step 1 finished in %s' % (step_1_time_string)

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.Barrier()  # Make sure all processes are here

    # Now re-initialise the temporary data set in read access mode only - - - -
    #
    self.tmp_dataset.re_initialise('read')

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      #f = open('tmp_data_set-dedup-'+str(parallel.rank())+'-'+ \
      #         str(parallel.size()),'w')
      #tmp_list = self.tmp_dataset.dict.keys()
      #tmp_list.sort()

      #for r in tmp_list:
      #  rec = self.tmp_dataset.dict[r]
      #  rec_items = rec.items()
      #  rec_items.sort()
      #  rec = str(r)+': '+str(rec_items)
      #  f.write(rec+os.linesep)
      #f.close()

      f = open('indexes-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index.num_indexes):
        tmp_index = self.blocking_index.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    print '1:'
    print '1:Step 2: Perform deduplication within blocks'
    print '1:-------'
    print '1:'

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    [rec_pair_dict, rec_pair_cnt] = \
           indexing.deduplication_rec_pairs(self.blocking_index)
    rec_pair_time = time.time() - tmp_time
    rec_pair_time_string = output.time_string(rec_pair_time)

    print '1:'
    print '1:  Built record pair dictionary with %i entries in %s' % \
          (rec_pair_cnt, rec_pair_time_string)

    # And do the comparisons of record pairs into classifer - - - - - - - - - -
    #
    [p] = do_comparison(self.tmp_dataset, self.tmp_dataset,
                        self.rec_comparator, self.classifier,
                        rec_pair_dict, rec_pair_cnt, self.block_size)

    # Now gather classifier results on process 0 and merge  - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier_results = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          print '1:    Received classifier from process %i and merged it' % (p)

          self.classifier.merge(tmp_classifier_results)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier.results, 0) # Send classifier results to
                                                  # process 0
        step_2_comm_time += (time.time() - tmp_time)
        print '1:    Sent classifier to process 0'

    # If run in parallel, broadcast the classifier results from process 0 - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.classifier.results, p)
          step_2_comm_time += (time.time() - tmp_time)
          print '1:    Sent classifier to process %i' % (p)

      else:
        tmp_time = time.time()
        self.classifier.results = parallel.receive(0)
        step_2_comm_time += (time.time() - tmp_time)
        print '1:    Received classifier from process 0'

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      tmp_list = rec_pair_dict.keys()
      tmp_list.sort()
      f = open('rec-pair-dict-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for rp in tmp_list:
        rec_list = rec_pair_dict[rp].items()
        rec_list.sort()
        r = str(rp)+': '+str(rec_list)
        f.write(r+os.linesep)
      f.close()

      tmp_list = self.classifier.results.keys()
      tmp_list.sort()
      f = open('classifier_results_dict-dedup-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for c in tmp_list:
        res = self.classifier.results[c].items()
        res.sort()
        ce = str(c)+': '+str(res)
        f.write(ce+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    step_2_time = time.time() - step_2_time  # Calculate time for step 2
    step_2_time_string = output.time_string(step_2_time)

    print '1:'
    print '1:Step 2 (deduplication) finished in %s' % (step_2_time_string)
    print '1:  Totally %i record pair comparisons' % (rec_pair_cnt)

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print '1:'
    print '1:Step 3: Output and assignment procedures'
    print '1:-------'
    print '1:'

    step_3_time = time.time()  # Get current time

    # Get the results dictionary with all the record pairs and their weights
    #
    results_dict = self.classifier.results

    if (results_dict == {}):
      print 'warning:Results dictionary empty'

    else:  # There are results

      # Do assignment restrictions if they are defined  - - - - - - - - - - - -
      #
      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('auction', results_dict, \
                                        'deduplication', self.output_threshold)
      else:  # No one-to-one assignment, set one2one result to None
        o2o_results_dict = None

      #################### START PARALLEL TEST CODE ###########################

      if (SAVE_PARALLEL_TEST_FILES == True):
        tmp_list = o2o_results_dict.items()
        tmp_list.sort()
        f = open('one2one-dedup-'+str(parallel.rank())+'-'+ \
                 str(parallel.size()),'w')
        for c in tmp_list:
          f.write(str(c)+os.linesep)
        f.close()

      #################### END PARALLEL TEST CODE #############################

      if (parallel.rank() == 0):  # Only processor 0 prints results

        # Print or save weights histogram - - - - - - - - - - - - - - - - - - -
        #
        if (self.output_histogram == True):
          output.histogram(results_dict)
        elif (self.output_histogram != False):
          output.histogram(results_dict, self.output_histogram)

        # Print or save detailed record pairs - - - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_details == True):
          output.rec_pair_details(self.tmp_dataset, self.tmp_dataset,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold)
        elif (self.output_rec_pair_details != False):
          output.rec_pair_details(self.tmp_dataset, self.tmp_dataset,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold,
                                  self.output_rec_pair_details)

        # Print or save record pairs with weights - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_weights == True):
          output.rec_pair_weights(self.tmp_dataset.name,
                                  self.tmp_dataset.name,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold)
        elif (self.output_rec_pair_weights != False):
          output.rec_pair_weights(self.tmp_dataset.name,
                                  self.tmp_dataset.name,
                                  results_dict, o2o_results_dict,
                                  self.output_threshold,
                                  self.output_rec_pair_weights)

    step_3_time = time.time() - step_3_time  # Calculate time for step 3
    step_3_time_string = output.time_string(step_3_time)

    print '1:'
    print '1:Step 3 (output and assignments) finished in %s' % \
          (step_3_time_string)
    print '1:'

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    total_time = time.time() - total_time  # Calculate total time

    total_time_string =       output.time_string(total_time)
    step_1_comm_time_string = output.time_string(step_1_comm_time)
    step_2_comm_time_string = output.time_string(step_2_comm_time)

    print '1:Total time needed for deduplication of %i records: %s' % \
          (self.number_records, total_time_string)
    print '1:  Time for step 1 (standardisation):       %s' % \
          (step_1_time_string)
    print '1:  Time for step 2 (deduplication):         %s' % \
          (step_2_time_string)
    print '1:  Time for step 3 (assignment and output): %s' % \
          (step_3_time_string)
    print '1:  Time for communication in step 1: %s' % \
          (step_1_comm_time_string)
    print '1:  Time for communication in step 2: %s' % \
          (step_2_comm_time_string)

  # ---------------------------------------------------------------------------

  def link(self, **kwargs):
    """Link the given two data set using the defined record standardisers,
       record comparators, blocking indexes and classifiers.

       Records are loaded block wise from the input data sets, then
       standardised (if the record standardisers are defined, otherwise an
       input data set is directly taken for the linkage process), linked and
       the results are printed and/or saved into the result file(s).

       If the arguments 'first_record_a' and/or 'first_record_b' are not given,
       they will automatically be set to the first record in the data sets
       (i.e. record number 0). Similarly, if the arguments 'number_records_a'
       and/or 'number_records_b' are not given, they will be set to the total
       number of records in the input data sets.

       The temporary data sets must be random acces data set implementations,
       i.e. either Shelve or a Memory data sets. For large data set it is
       recommended to use Shelve data sets. This temporary data sets have to be
       initialised in access mode 'readwrite'.

       Currently, the output can be a printed or saved list of record pairs in
       both a detailed and condensed form (if the arguments
       'output_rec_pair_details' and 'output_rec_pair_weights' are set to
       'True' or to a file name (a string). The output can be filtered by
       setting the 'output_threshold' (meaning all record pairs with a weight
       less then this threshold are not printed or saved).

       In future versions, it will be possible to compile an output data set.

       A histogram can be saved or printed by setting the argument
       'output_histogram' to 'True' or to a file name.

       It is also possible to apply a one-to-one assignment procedure by
       setting the argument 'output_assignment' to 'one2one'.
    """

    self.input_dataset_a =    None   # A reference to the first input data set
    self.tmp_dataset_a =      None   # A reference to the first temporary
                                     # (random access) data set
    self.input_dataset_b =    None   # A reference ot the second input data set
    self.tmp_dataset_b =      None   # A reference to the second temporary
                                     # (random access) data set
#    self.output_dataset =     None   # A reference to the output data set

    self.rec_standardiser_a = None   # Reference to a record standardiser for
                                     # the first data set (A)
    self.rec_standardiser_b = None   # Reference to a record standardiser for
                                     # the second data set (B)
    self.blocking_index_a =   None   # Reference to a blocking index for data
                                     # the first set (A)
    self.blocking_index_b =   None   # Reference to a blocking index for data
                                     # the second set (B)
    self.rec_comparator =     None   # Reference to a record comparator
    self.classifier =         None   # Reference to a weight vector classifier

    self.first_record_a =     None   # Number of the first record to process
                                     # in the first data set (A)
    self.number_records_a =   None   # Number of records to process in the
                                     # first data set (A)
    self.first_record_b =     None   # Number of the first record to process
                                     # in the second data set (B)
    self.number_records_b =   None   # Number of records to process in the
                                     # second data set (B)

    self.output_histogram = False         # Set to True, a file name or False
                                          # (default) if a histogram of weights
                                          # should be printed or saved
    self.output_rec_pair_details = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved in details
    self.output_rec_pair_weights = False  # Set to True, a file name or False
                                          # (default) if record pairs should
                                          # be printed or saved with weights
    self.output_threshold = None          # Set to a weight threshold (only
                                          # record pairs with weights equal to
                                          # or above will be saved and or
                                          # printed)
    self.output_assignment = None         # Set to 'one2one' if one-to-one
                                          # assignment should be forced
                                          # (default: None)

    for (keyword, value) in kwargs.items():
      if (keyword == 'input_dataset_a'):
        self.input_dataset_a = value
      elif (keyword == 'input_dataset_b'):
        self.input_dataset_b = value
      elif (keyword == 'tmp_dataset_a'):
        self.tmp_dataset_a = value
      elif (keyword == 'tmp_dataset_b'):
        self.tmp_dataset_b = value
#      elif (keyword == 'output_dataset'):
#        self.output_dataset = value

      elif (keyword == 'rec_standardiser_a'):
        self.rec_standardiser_a = value
      elif (keyword == 'rec_standardiser_b'):
        self.rec_standardiser_b = value
      elif (keyword == 'rec_comparator'):
        self.rec_comparator = value
      elif (keyword == 'blocking_index_a'):
        self.blocking_index_a = value
      elif (keyword == 'blocking_index_b'):
        self.blocking_index_b = value
      elif (keyword == 'classifier'):
        self.classifier = value

      elif (keyword == 'first_record_a'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record_a" is not a valid integer number'
          raise Exception
        self.first_record_a = value
      elif (keyword == 'first_record_b'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "first_record_b" is not a valid integer number'
          raise Exception
        self.first_record_b = value
      elif (keyword == 'number_records_a'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records_a" is not a positive '+ \
                'integer number'
          raise Exception
        self.number_records_a = value
      elif (keyword == 'number_records_b'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "number_records_b" is not a positive '+ \
                'integer number'
          raise Exception
        self.number_records_b = value

      elif (keyword == 'output_rec_pair_details'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_rec_pair_details" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_rec_pair_details = value
      elif (keyword == 'output_rec_pair_weights'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_rec_pair_weights" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_rec_pair_weights = value
      elif (keyword == 'output_histogram'):
        if (not isinstance(value, str)) and (value not in [True, False]):
          print 'error:Argument "output_histogram" must be ' + \
                'a file name or "True" or "False"'
          raise Exception
        self.output_histogram = value
      elif (keyword == 'output_threshold'):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "output_threshold" is not a number: %s' % \
                (str(value))
        self.output_threshold = value
      elif (keyword == 'output_assignment'):
        if (value not in ['one2one', None]):
          print 'error:Illegal value for argument "output_assignment": %s' % \
                (str(value))
          raise Exception
        else:
          self.output_assignment = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset_a == None):
      print 'error:Input data set A is not defined'
      raise Exception
    if (self.input_dataset_b == None):
      print 'error:Input data set B is not defined'
      raise Exception

    if (self.tmp_dataset_a == None):
      print 'error:Temporary data set A is not defined'
      raise Exception
    elif (self.tmp_dataset_a.dataset_type not in ['SHELVE', 'MEMORY']):
      print 'error:Temporary data set A must be a random access data set' + \
            ' (either Shelve or Memory)'
      raise Exception
    if (self.tmp_dataset_a.access_mode not in ['write','append','readwrite']):
      print 'error:Temporary data set A must be initialised in one of the ' + \
            'access  modes: "write", "append", or "readwrite"'
      raise Exception

    if (self.tmp_dataset_b == None):
      print 'error:Temporary data set B is not defined'
      raise Exception
    elif (self.tmp_dataset_b.dataset_type not in ['SHELVE', 'MEMORY']):
      print 'error:Temporary data set B must be a random access data set' + \
            ' (either Shelve or Memory)'
      raise Exception
    if (self.tmp_dataset_b.access_mode not in ['write','append','readwrite']):
      print 'error:Temporary data set B must be initialised in one of the ' + \
            'access  modes: "write", "append", or "readwrite"'
      raise Exception

    # Check if there are file names for the temporary data sets and - - - - - -
    # if they differ
    #
    tmp_file_name_a = getattr(self.tmp_dataset_a, 'file_name', None)
    tmp_file_name_b = getattr(self.tmp_dataset_b, 'file_name', None)
    if (tmp_file_name_a != None) and (tmp_file_name_b != None):
      if (tmp_file_name_a == tmp_file_name_b):
        print 'error:The same file names for both temporary data sets'
        raise Exception

#    if (self.output_dataset == None):
#      print 'error:Output data set is not defined'
#      raise Exception

    # Make sure at least one output is defined
    #
    if (self.output_rec_pair_weights == False) and \
       (self.output_rec_pair_details == False) and \
       (self.output_histogram == False):
      print 'error:No ouput of results is defined.'
      raise Exception
    #
    # Code above to be removed once output data set functionality implemented

    if (self.first_record_a == None):
      self.first_record_a = 0  # Take default first record in data set
    if (self.first_record_b == None):
      self.first_record_b = 0  # Take default first record in data set

    if (self.number_records_a == None):  # Take all records
      self.number_records_a = self.input_dataset_a.num_records
    if (self.number_records_b == None):  # Take all records
      self.number_records_b = self.input_dataset_b.num_records

    if (self.rec_comparator == None):
      print 'error:No record comparator defined'
      raise Exception
    if (self.rec_comparator.dataset_a != self.tmp_dataset_a) or \
       (self.rec_comparator.dataset_b != self.tmp_dataset_b):
      print 'error:Illegal data set definition in record comparator'
      raise Exception

    if (self.blocking_index_a == None):
      print 'error:No blocking index for data set A defined'
      raise Exception
    if (self.blocking_index_b == None):
      print 'error:No blocking index for data set B defined'
      raise Exception

    if (self.classifier == None):
      print 'error:No classifier defined'
      raise Exception
    if (self.classifier.dataset_a != self.tmp_dataset_a) or \
       (self.classifier.dataset_b != self.tmp_dataset_b):
      print 'error:Illegal data set definition in classifier'
      raise Exception

#    if (self.rec_standardiser_a != None):
#      if (self.rec_standardiser_a.input_dataset != self.input_dataset_a):
#        print 'error:Illegal input data set definition in record '+ \
#              'standardiser A: %s (should be: %s)' % \
#              (str(self.rec_standardiser_a.input_dataset.name), \
#               str(self.input_dataset_a.name))
#        raise Exception
#      if (self.rec_standardiser_a.output_dataset != self.tmp_dataset_a):
#        print 'error:Illegal output data set definition in record '+ \
#              'standardiser A: %s (should be: %s)' % \
#              (str(self.rec_standardiser_a.output_dataset.name), \
#               str(self.tmp_dataset_a.name))
#        raise Exception

#    else:  # No standardiser for data set A defined, so field names in input
#           # and temporary data sets must be the same
#      input_field_name_list = self.input_dataset_a.fields.keys()
#      input_field_name_list.sort()
#      tmp_field_name_list = self.tmp_dataset_a.fields.keys()
#      tmp_field_name_list.sort()
#
#      if (input_field_name_list != tmp_field_name_list):
#        print 'error:Field names differ in input and temporary data sets ' + \
#              '(with no record standardiser for data set A defined)'
#
#    if (self.rec_standardiser_b != None):
#      if (self.rec_standardiser_b.input_dataset != self.input_dataset_b):
#        print 'error:Illegal input data set definition in record '+ \
#              'standardiser B: %s (should be: %s)' % \
#              (str(self.rec_standardiser_b.input_dataset.name), \
#               str(self.input_dataset_b.name))
#        raise Exception
#      if (self.rec_standardiser_b.output_dataset != self.tmp_dataset_b):
#        print 'error:Illegal output data set definition in record '+ \
#              'standardiser B: %s (should be: %s)' % \
#              (str(self.rec_standardiser_b.output_dataset.name), \
#               str(self.tmp_dataset_b.name))
#        raise Exception
#
#    else:  # No standardiser for data set B defined, so field names in input
#           # and temporary data sets must be the same
#      input_field_name_list = self.input_dataset_b.fields.keys()
#      input_field_name_list.sort()
#      tmp_field_name_list = self.tmp_dataset_b.fields.keys()
#      tmp_field_name_list.sort()
#
#      if (input_field_name_list != tmp_field_name_list):
#        print 'error:Field names differ in input and temporary data sets ' + \
#              '(with no record standardiser for data set B defined)'
#
#    if (self.blocking_index_a.dataset != self.tmp_dataset_a):
#      print 'error:Illegal data set definition in blocking index A'
#      raise Exception
#    if (self.blocking_index_b.dataset != self.tmp_dataset_b):
#      print 'error:Illegal data set definition in blocking index B'
#      raise Exception

    total_time = time.time()  # Get current time

    print '1:'
    print '1:***** Link data set: %s (type "%s") with data set: %s ' % \
           (self.input_dataset_a.name, self.input_dataset_a.dataset_type, \
           self.input_dataset_b.name) + '(type "%s")' % \
           (self.input_dataset_b.dataset_type)
    print '1:*****   Temporary data set A: "%s" (type "%s")' % \
          (self.tmp_dataset_a.name, self.tmp_dataset_a.dataset_type)
    print '1:*****   Temporary data set B: "%s" (type "%s")' % \
          (self.tmp_dataset_b.name, self.tmp_dataset_b.dataset_type)
    print '1:'
    print '1:Step 1: Loading, standardisation and indexing'
    print '1:-------'
    print '1:'

    step_1_time = time.time()  # Get current time

    # Call the main standardisation routine for data set A  - - - - - - - - - -
    #
    [p, step_1_comm_time] = do_load_standard_indexing(self.input_dataset_a,
                                                      self.tmp_dataset_a,
                                                      self.rec_standardiser_a,
                                                      self.blocking_index_a,
                                                      self.first_record_a,
                                                      self.number_records_a,
                                                      self.block_size)

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index A from process %i' % (p)

        self.blocking_index_a.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index_a.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:    Sent index A to process 0'

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_a.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent index A to process %i' % (p)

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index A from process 0'

        self.blocking_index_a.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index_a.compact()

    print '1:'

    # Call the main standardisation routine for data set B  - - - - - - - - - -
    #
    [p, tmp_time] = do_load_standard_indexing(self.input_dataset_b,
                                              self.tmp_dataset_b,
                                              self.rec_standardiser_b,
                                              self.blocking_index_b,
                                              self.first_record_b,
                                              self.number_records_b,
                                              self.block_size)
    step_1_comm_time += tmp_time

    # If Febrl is run in parallel, collect blocking index in process 0  - - - -
    #
    if (parallel.rank() == 0):
      for p in range(1, parallel.size()):
        tmp_time = time.time()
        tmp_indexes = parallel.receive(p)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index B from process %i' % (p)

        self.blocking_index_b.merge(tmp_indexes)

    else:  # Send index to process 0
      tmp_time = time.time()
      parallel.send(self.blocking_index_b.index, 0) # Send indexes to process 0
      step_1_comm_time += (time.time() - tmp_time)
      print '1:    Sent index B to process 0'

    # If run in parallel, broadcast the blocking index from process 0 - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.blocking_index_b.index, p)
          step_1_comm_time += (time.time() - tmp_time)
          print '1:    Sent index B to process %i' % (p)

      else:
        tmp_time = time.time()
        tmp_indexes = parallel.receive(0)
        step_1_comm_time += (time.time() - tmp_time)
        print '1:    Received index B from process 0'

        self.blocking_index_b.merge(tmp_indexes)

    # Compact the blocking index  - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.blocking_index_b.compact()

    step_1_time = time.time() - step_1_time  # Calculate time for step 1
    step_1_time_string = output.time_string(step_1_time)

    print '1:'
    print '1:Step 1 finished in %s' % (step_1_time_string)

    # End of step 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.Barrier()  # Make sure all processes are here

    # Now re-initialise the temporary data sets in read access mode only  - - -
    #
    self.tmp_dataset_a.re_initialise('read')
    self.tmp_dataset_b.re_initialise('read')

    #################### START PARALLEL TEST CODE #############################
    # Save temporary data sets and indexes to files (on all processes)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      #f = open('tmp_data_set_a-link-'+str(parallel.rank())+'-'+ \
      #         str(parallel.size()),'w')
      #tmp_list = self.tmp_dataset_a.dict.keys()
      #tmp_list.sort()

      #for r in tmp_list:
      #  rec = self.tmp_dataset_a.dict[r]
      #  rec_items = rec.items()
      #  rec_items.sort()
      #  rec = str(r)+': '+str(rec_items)
      #  f.write(rec+os.linesep)
      #f.close()

      f = open('indexes_a-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index_a.num_indexes):
        tmp_index = self.blocking_index_a.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index_a.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

      f = open('indexes_b-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')

      for i in range (self.blocking_index_b.num_indexes):
        tmp_index = self.blocking_index_b.index[i].keys()
        tmp_index.sort()

        for bi in tmp_index:
          ind = self.blocking_index_b.index[i][bi]
          ind_list = ind.items()
          ind_list.sort()

          ii = str(i)+'_'+str(bi)+': '+str(ind_list)
          f.write(ii+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    print '1:'
    print '1:Step 2: Perform linkage within blocks'
    print '1:-------'
    print '1:'

    step_2_time = time.time()  # Get current time
    step_2_comm_time = 0.0

    # Get the record pairs which have to be compared  - - - - - - - - - - - - -
    #
    tmp_time = time.time()
    [rec_pair_dict, rec_pair_cnt] = \
           indexing.linkage_rec_pairs(self.blocking_index_a,
                                      self.blocking_index_b)
    rec_pair_time = time.time() - tmp_time
    rec_pair_time_string = output.time_string(rec_pair_time)

    print '1:'
    print '1:  Built record pair dictionary with %i entries in %s' % \
          (rec_pair_cnt, rec_pair_time_string)

    # And do the comparisons of record pairs into classifer - - - - - - - - - -
    #
    [p] = do_comparison(self.tmp_dataset_a, self.tmp_dataset_b,
                        self.rec_comparator, self.classifier,
                        rec_pair_dict, rec_pair_cnt, self.block_size)

    # Now gather classifier results on process 0 and merge  - - - - - - - - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          tmp_classifier_results = parallel.receive(p)
          step_2_comm_time += (time.time() - tmp_time)
          print '1:    Received classifier from process %i and merged it' % (p)

          self.classifier.merge(tmp_classifier_results)

      else:
        tmp_time = time.time()
        parallel.send(self.classifier.results, 0) # Send classifier results to
                                                  # process 0
        step_2_comm_time += (time.time() - tmp_time)
        print '1:    Sent classifier to process 0'

    # If run in parallel, broadcast the classifier results from process 0 - - -
    #
    if (parallel.size() > 1):
      if (parallel.rank() == 0):
        for p in range(1, parallel.size()):
          tmp_time = time.time()
          parallel.send(self.classifier.results, p)
          step_2_comm_time += (time.time() - tmp_time)
          print '1:    Sent classifier to process %i' % (p)

      else:
        tmp_time = time.time()
        self.classifier.results = parallel.receive(0)
        step_2_comm_time += (time.time() - tmp_time)
        print '1:    Received classifier from process 0'

    #################### START PARALLEL TEST CODE #############################
    # Save classifiers and weight vectors to files (only process 0)
    #
    if (SAVE_PARALLEL_TEST_FILES == True):
      tmp_list = rec_pair_dict.keys()
      tmp_list.sort()
      f = open('rec-pair-dict-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for rp in tmp_list:
        rec_list = rec_pair_dict[rp].items()
        rec_list.sort()
        r = str(rp)+': '+str(rec_list)
        f.write(r+os.linesep)
      f.close()

      tmp_list = self.classifier.results.keys()
      tmp_list.sort()
      f = open('classifier_results_dict-link-'+str(parallel.rank())+'-'+ \
               str(parallel.size()),'w')
      for c in tmp_list:
        res = self.classifier.results[c].items()
        res.sort()
        ce = str(c)+': '+str(res)
        f.write(ce+os.linesep)
      f.close()

    #################### END PARALLEL TEST CODE ###############################

    step_2_time = time.time() - step_2_time  # Calculate time for step 2
    step_2_time_string = output.time_string(step_2_time)

    print '1:'
    print '1:Step 2 (linkage) finished in %s' % (step_2_time_string)
    print '1:  Totally %i record pair comparisons' % (rec_pair_cnt)

    # Output the results  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print '1:'
    print '1:Step 3: Output and assignment procedures'
    print '1:-------'
    print '1:'

    step_3_time = time.time()  # Get current time

    # Get the results dictionary with all the record pairs and their weights
    #
    results_dict = self.classifier.results

    if (results_dict == {}):
      print 'warning:Results dictionary empty'

    else:  # There are results

      # Do assignment restrictions if they are defined  - - - - - - - - - - - -
      #
      if (self.output_assignment != None):  # An output assignment is defined

        if (self.output_assignment == 'one2one'):

          # Do a one-to-one assignment on the classifier results dict
          #
          o2o_results_dict = lap.do_lap('auction', results_dict, \
                                        'linkage', self.output_threshold)
      else:  # No one-to-one assignment, set one2one result to None
        o2o_results_dict = None

      #################### START PARALLEL TEST CODE ###########################

      if (SAVE_PARALLEL_TEST_FILES == True):
        tmp_list = o2o_results_dict.items()
        tmp_list.sort()
        f = open('one2one-link-'+str(parallel.rank())+'-'+ \
                 str(parallel.size()),'w')
        for c in tmp_list:
          f.write(str(c)+os.linesep)
        f.close()

      #################### END PARALLEL TEST CODE #############################

      if (parallel.rank() == 0):  # Only processor 0 prints results

        # Print or save weights histogram - - - - - - - - - - - - - - - - - - -
        #
        if (self.output_histogram == True):
          output.histogram(results_dict)
        elif (self.output_histogram != False):
          output.histogram(results_dict, self.output_histogram)

        # Print or save detailed record pairs - - - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_details == True):
          output.rec_pair_details(self.tmp_dataset_a, self.tmp_dataset_b, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)
        elif (self.output_rec_pair_details != False):
          output.rec_pair_details(self.tmp_dataset_a, self.tmp_dataset_b, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold, \
                                  self.output_rec_pair_details)

        # Print or save record pairs with weights - - - - - - - - - - - - - - -
        #
        if (self.output_rec_pair_weights == True):
          output.rec_pair_weights(self.tmp_dataset_a.name, \
                                  self.tmp_dataset_b.name, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold)
        elif (self.output_rec_pair_weights != False):
          output.rec_pair_weights(self.tmp_dataset_a.name, \
                                  self.tmp_dataset_b.name, \
                                  results_dict, o2o_results_dict, \
                                  self.output_threshold, \
                                  self.output_rec_pair_weights)

    step_3_time = time.time() - step_3_time  # Calculate time for step 3
    step_3_time_string = output.time_string(step_3_time)

    print '1:'
    print '1:Step 3 (output and assignments) finished in %s' % \
          (step_3_time_string)
    print '1:'

    parallel.Barrier()  # Wait here for all processes - - - - - - - - - - - - -

    total_time = time.time() - total_time  # Calculate total time

    total_time_string =       output.time_string(total_time)
    step_1_comm_time_string = output.time_string(step_1_comm_time)
    step_2_comm_time_string = output.time_string(step_2_comm_time)

    print '1:Total time needed for linkage of %i records with %i records: %s' \
          % (self.number_records_a, self.number_records_b, total_time_string)
    print '1:  Time for step 1 (standardisation):       %s' % \
          (step_1_time_string)
    print '1:  Time for step 2 (linkage):               %s' % \
          (step_2_time_string)
    print '1:  Time for step 3 (assignment and output): %s' % \
          (step_3_time_string)
    print '1:  Time for communication in step 1: %s' % \
          (step_1_comm_time_string)
    print '1:  Time for communication in step 2: %s' % \
          (step_2_comm_time_string)

# =============================================================================

class ProjectLog:
  """Class for Febrl project logs.

     Project logs capture print statements with a special form as well as
     Python exceptions and warnings, and writes them to a log file and prints
     them to the standard terminal output.

     A log file name can be given with the attribute 'file_name' when a project
     log is initialised. If no file name is given, no logging information is
     written to a file (but still printed to standard output).

     When a project log is specified, the user needs to set both a verbose and
     a log level in the range 0 to 3 using the attributes 'log_level' and
     'verbose_level'. A level of 0 means nothing no log messages are logged or
     printed to standard output, level one means only important messages are
     logged/printed, level 2 means medium output and level 3 finally means a
     high volume output.

     Error messages are looged and printed in any case.

     Logging and printing of warning massages can be suppressed by setting the
     flag 'no_warnings' to 'True' when the project log file is initialised.
     The default is 'False' meaning that warning messages are printed.

     With the flag 'clear_log' the user can clear the contents of a log file
     when it is initialised by setting the value of this flag to 'True'. The
     default is 'False', in which case logging messages are appended to the
     given log file.

     Within Febrl modules, log messages are normal Python statements that must
     start with either of the following substrings followed by the message:

       '1:'        (a high priority log message)
       '2:'        (a medium priority log message)
       '3:'        (a low priority log message)
       'warning:'  (a warning message)
       'error:'    (an error message)

     All other print statements are handled like normal prints and sent to
     standard output 'sys.stdout'.

     The 'parallel_print' defines how printing and logging is done when Febrl
     is run in parallel. Possible are either 'host' (the default) in which case
     only the host node (where Febrl was started on) is printing and logging
     messages, or 'all' in which case all processes are printing and logging.
     Note that currently both error and warning messages are printed by all
     processes.
     The 'parallel_print' is used to set 'printmode' in the parallel.py module.

     The 'message' is the string that is printed to standard output and written
     to the log file if 'level' is equal to smaller than the 'verbose_level' or
     the 'log_level' respectively.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor - Create a new project log object.
    """

    self.febrl =          None    # A reference to the febrl object.
    self.project =        None    # A reference to the project object.
    self.log_level =      0       # Level of looging to file.
    self.verbose_level =  0       # Level of verbose printing.
    self.no_warnings =    False   # Supress warning messages
    self.clear_log =      False   # A flag (True/False), if set to True all
                                  # content of the log will be cleared first.
                                  # Default value is False, i.e. messages will
                                  # be appended to the existing log file.
    self.file_name =      None    # The log file name.
    self.log_file =       None    # The log file pointer
    self.parallel_print = 'host'  # Either set to 'host' (default) or 'all'

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'project'):
        self.project = value
        if (self.project.febrl != None):
          self.febrl = self.project.febrl  # Reference to the febrl object
        else:
          print 'error:Febrl object not defined'
          raise Exception

      elif (keyword == 'file_name'):
        if (not isinstance(value,str)):
          print 'error:Argument "file_name" is not a string'
          raise Exception
        self.file_name = value

      elif (keyword == 'clear_log'):
        if (value not in [True, False]):
          print 'error:Argument "clear_log" must be "True" or "False"'
          raise Exception
        self.clear_log = value

      elif (keyword in ['no_warnings','no_warn']):
        if (value not in [True, False]):
          print 'error:Argument "no_warnings" must be "True" or "False"'
          raise Exception
        self.no_warnings = value

      elif (keyword == 'parallel_print'):
        if (value not in ['host','all']):
          print 'error:Argument "parallel_print" must be set to "host"'+ \
                ' or "all"'
          raise Exception
        self.parallel_print = value

      elif (keyword == 'log_level'):
        if (not isinstance(value,int)) or (value < 0) or (value > 3):
          print 'error:Argument "log_level" must zero or a positive integer'
          raise Exception
        self.log_level = value

      elif (keyword == 'verbose_level'):
        if (not isinstance(value,int)) or (value < 0) or (value > 3):
          print 'error:Argument "verbose_level" must zero or a positive '+ \
                'integer'
          raise Exception
        self.verbose_level = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.project == None):
      print 'error:Project not defined'
      raise Exception

    # Set the parallel printing mode  - - - - - - - - - - - - - - - - - - - - -
    #
    parallel.printmode = self.parallel_print

    # Check if file can be opened - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.file_name != None):  # File logging activated

      if (self.clear_log == True):  # Open log file in 'write' mode
        try:
          self.log_file = open(self.file_name,'w')
        except:
          print 'error:Can not open log file "%s" for writing' % \
                (self.file_name)
          raise IOError

      else:  # Open log file in 'append' mode (don't delete old content)
        try:
          self.log_file = open(self.file_name,'a')
        except:
          print 'error:Can not open log file "%s" for appending' % \
                (self.file_name)
          raise IOError

      # Write a header to the log file - - - - - - - - - - - - - - - - - - - -
      #
      self.append('#'*75)
      self.append('# Febrl project log file')
      self.append('#')
      self.append('# Project name: '+self.project.name)
      self.append('# Project log file name: '+self.file_name)
      self.append('#')
      self.append('# Time and date: '+time.strftime('%d %b %Y %H:%M:%S', \
                                      time.localtime(time.time())))
      self.append('#')
      self.append('# Febrl version:   '+self.febrl.version)
      self.append('# Febrl license:   '+self.febrl.license)
      self.append('# Febrl copyright: '+self.febrl.copyright)
      self.append('#')
      self.append('# Febrl description: '+self.febrl.description)
      self.append('# Febrl path:        '+self.febrl.febrl_path)
      self.append('#'*75)
      self.append('')
      self.flush()

      # Redirect the system stdout so prints can be logged - - - - - - - - - -
      #
      sys.stdout = LogPrinter(self)

      # Set the system exception hook so exceptions and warnings are logged - -
      #
      sys.excepthook = self.except_hook

  # ---------------------------------------------------------------------------

  def except_hook(self, exc_type, value, trace_back):
    """A routine to catch exceptions and warnings and process them.
    """

    bug_report1 = 'Please submit an error report by sending an e-mail'+ \
                  ' to the Febrl authors'
    bug_report2 = 'and attach this error message.'

    time_stamp = time.strftime('%d %b %Y %H:%M:%S',time.localtime(time.time()))

    # Get complete trace stack
    #
    trace_list = traceback.extract_tb(trace_back)
    trace_stack_size = len(trace_list)

    # Create a message list (one element is one line) - - - - - - - - - - - - -
    #
    except_msg = []
    except_msg.append(parallel.prompt+'#'*75)
    except_msg.append(parallel.prompt+'### Exception: '+str(exc_type))
    except_msg.append(parallel.prompt+'###   Time:       '+time_stamp)
    except_msg.append(parallel.prompt+'###   Message:    '+str(value))
    except_msg.append(parallel.prompt+'###   Trace stack:')
    for lev in range(trace_stack_size):
      spc = '  '*lev
      except_msg.append(parallel.prompt+'###     '+spc+'-'*(67-lev*2))
      except_msg.append(parallel.prompt+'###     '+spc+'Module:   '+ \
                  str(trace_list[lev][0]))
      except_msg.append(parallel.prompt+'###     '+spc+'Function: '+ \
                  str(trace_list[lev][2]))
      except_msg.append(parallel.prompt+'###     '+spc+'Line:     '+ \
                  str(trace_list[lev][1]))
      except_msg.append(parallel.prompt+'###     '+spc+'Text:     '+ \
                  str(trace_list[lev][3]))
    except_msg.append(parallel.prompt+'#'*75)
    except_msg.append(parallel.prompt+'### '+bug_report1)
    except_msg.append(parallel.prompt+'### '+bug_report2)
    except_msg.append(parallel.prompt+'#'*75)

    # Print and log the message - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for msg_line in except_msg:
      self.append(msg_line)
      sys.__stdout__.write(msg_line+os.linesep)

    # And flush the standard out and the log file - - - - - - - - - - - - - - -
    #
    self.flush()
    sys.__stdout__.flush()

  # ---------------------------------------------------------------------------

  def close(self):
    """Close the log file if it was opened.
    """

    if (self.log_file != None):
      self.log_file.close()
      self.log_file = None

  # ---------------------------------------------------------------------------

  def append(self, message):
    """Write the given message to the log file if logging is activated (i.e. if
       the log file is opened).

       The message can be a string or a list of strings, in which case each
       list element will be written as one line. A line separator is appended
       to the end of each line.
    """

    if (self.log_file != None):  # Only if the log file is open
      if (isinstance(message, str)):
        if (len(message) > 0) and (message[-1] == os.linesep):
          self.log_file.write(message)
        else:
          self.log_file.write(message+os.linesep)

      elif (isinstance(message, list)):
        for m in message:
          if (isinstance(m, str)):
            if (len(m) > 0) and (message[-1] == os.linesep):
              self.log_file.write(m)
            else:
              self.log_file.write(m+os.linesep)

          else:
            print 'error:Element in log file message list is not a ' + \
                  'string: %s' % (str(m))
            raise Exception
      else:
        print 'error:Log file message is not a string or a list: %s' % \
              (str(message))
        raise Exception

  # ---------------------------------------------------------------------------

  def flush(self):
    """Flush the log file if it is opened to make sure it is written out.
    """

    if (self.log_file != None):
      self.log_file.flush()

  # ---------------------------------------------------------------------------

  def print_log(self, print_type):
    """Print the contents of the log file in various types:
       Text, HTML, LaTex, XML
    """

    pass

    # Check if file is open

    # re-open in read mode

    # print all lines in format

    # close in read mode

# =============================================================================

class LogPrinter:
  """Class that replaces sys.stdout

     Each normal print statement is analysed in the 'write' method, if it
     starts with a 'error:', 'warning:' or '1:', '2:', '3:' it will also be
     passed to the project logger.
     Otherwise, it will simply be given to the original standard out.

     For parallel runs, error and warning message will be printed (and logged)
     from all processes, but normal verbose message are printed acording to the
     value of 'printmode' ('all' or 'host') as defined in parallel.py

     Messages which are neither error nor warning message, nor have a level at
     the beginning are also printed on all processes (if Febrl is run in
     parallel).
  """

  def __init__(self, project_log):
    self.project_log = project_log

    # Now check if printing of level messages is done - - - - - - - - - - - - -
    #
    if ((parallel.printmode == 'host') and (parallel.rank() == 0)) or \
       (parallel.printmode == 'all'):
      self.level_print = True  # Do print level messages
    else:
      self.level_print = False  # Don't print level messages

  def write(self, msg): # - - - - - - - - - - - - - - - - - - - - - - - - - - -

    #if (len(msg) > 0) and (msg[-1] != os.linesep):  # Append a line separator
    #  msg += os.linesep

    call_function =   sys._getframe(1).f_code.co_name
    call_linenumber = str(sys._getframe(1).f_lineno)
    call_module =     sys._getframe(1).f_code.co_filename
    time_stamp = time.strftime('%d %b %Y %H:%M:%S',time.localtime(time.time()))

    # Get message type and split message at line separators
    #
    if (msg[:6].lower() == 'error:'):  # An error message - - - - - - - - - - -

      msg = msg[6:]
      msg_list = msg.split(os.linesep)

      bug_report1 = 'Please submit an error report by sending an e-mail'+ \
                    ' to the Febrl authors'
      bug_report2 = 'and attach this error message.'

      error_msg = []  # Create a message list (one element is one line)

      error_msg.append(parallel.prompt+'#'*75)
      error_msg.append(parallel.prompt+'### Error')
      error_msg.append(parallel.prompt+'###   Module:  '+ call_module + \
                       ', function: ' + call_function + ', line number: ' + \
                       call_linenumber)
      error_msg.append(parallel.prompt+'###   Time:    '+time_stamp)
      error_msg.append(parallel.prompt+'###   Message: '+msg_list[0])
      for msg_line in msg_list[1:]:
        error_msg.append(parallel.prompt+'###            '+msg_line)
      error_msg.append(parallel.prompt+'#'*75)
      error_msg.append(parallel.prompt+'### '+bug_report1)
      error_msg.append(parallel.prompt+'### '+bug_report2)
      error_msg.append(parallel.prompt+'#'*75)

      for msg_line in error_msg:  # Print and log the message
        self.project_log.append(msg_line)
        sys.__stdout__.write(msg_line+os.linesep)

      self.project_log.flush()  # And flush the standard out and the log file
      sys.__stdout__.flush()

    elif (msg[:8].lower() == 'warning:'):  # A warning message  - - - - - - - -

      if (self.project_log.no_warnings == False):

        msg = msg[8:]
        msg_list = msg.split(os.linesep)

        warn_msg = []  # Create a message list (one element is one line)

        warn_msg.append(parallel.prompt)
        warn_msg.append(parallel.prompt+'### Warning')
        warn_msg.append(parallel.prompt+'###   Module:  '+  call_module + \
                        ', function: ' + call_function + ', line number: ' + \
                        call_linenumber)
        warn_msg.append(parallel.prompt+'###   Time:    '+time_stamp)
        warn_msg.append(parallel.prompt+'###   Message: '+msg_list[0])
        for msg_line in msg_list[1:]:
          warn_msg.append(parallel.prompt+'###            '+msg_line)

      for msg_line in warn_msg:  # Print and log the message
        self.project_log.append(msg_line)
        sys.__stdout__.write(msg_line+os.linesep)

      self.project_log.flush()  # And flush the standard out and the log file
      sys.__stdout__.flush()

    # A verbose level message - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    elif (msg[0] in '123') and (len(msg) >= 2) and (msg[1] == ':'): 

      if (self.level_print == True):

        msg_level = int(msg[0])
        msg = msg[2:]
        msg_list = msg.split(os.linesep)

        # Check if the level of the message is good for printing and/or logging
        #
        if (msg_level <= self.project_log.log_level):
          for msg_line in msg_list:
            self.project_log.append(parallel.prompt+msg_line)
          self.project_log.flush()

        if (msg_level <= self.project_log.verbose_level):
          for msg_line in msg_list:
            sys.__stdout__.write(parallel.prompt+msg_line+os.linesep)
          sys.__stdout__.flush()

    else:  # Print 'normal' print commands  - - - - - - - - - - - - - - - - - -

      msg_list = msg.split(os.linesep)

      for msg_line in msg_list:

        if ((len(msg_line) == 1) and (ord(msg_line) != 10)) or \
            (len(msg_line) > 1):
          if (sys.platform[:3] == 'win'): # No line separator needed on Windows
            sys.__stdout__.write(parallel.prompt+msg_line)
          else:
            sys.__stdout__.write(parallel.prompt+msg_line+os.linesep)
      sys.__stdout__.flush()

# =============================================================================
# The following are the main routines for standardisation and record linkage as
# used withing the project methods 'standardise', 'deduplicate' and 'link'

def do_load_standard_indexing(input_dataset, output_dataset,
                              record_standardiser, blocking_index,
                              first_record, number_records,
                              febrl_block_size):

  """The main routine that does the loading of records from the input data set
     into the output data set, if a record standardsier is given these records
     are cleaned and standardised, and if a blocking index is given such a
     index will be built and returned as well.

     If no standardisation is needed the argument 'record_standardiser' has to
     be set to None, and if no indexing is needed the 'blocking_index' argument
     has to be set to None.

     The output dataset must not be a Memory based data set. And of course the
     output data set must be initialised in access mode "write", "append" or
     "readwrite".

     If run in parallel, all processes will open the input data set and read
     records from it (and clean and standardised them), and the processed
     records will then be sent to process 0 (the host process) and saved into
     the output data set.
  """

  # Check if the given record numbers are valid - - - - - - - - - - - - - - - -
  #
  last_record = first_record + number_records

  if (first_record < 0) or (last_record > input_dataset.num_records):
    print 'error:Record range too large: (%i,%i)' % \
          (first_record, last_record)
    raise Exception

  # Check if the data sets are set correctly within the record standardiser - -
  #
  if (record_standardiser != None):
    if (record_standardiser.input_dataset != input_dataset):
      print 'error:Illegal input data set definition in record '+ \
            'standardiser: %s (should be: %s)' % \
            (str(record_standardiser.input_dataset.name), \
            str(input_dataset.name))
      raise Exception
    if (record_standardiser.output_dataset != output_dataset):
      print 'error:Illegal output data set definition in record '+ \
            'standardiser: %s (should be: %s)' % \
            (str(record_standardiser.output_dataset.name), \
             str(output_dataset.name))
      raise Exception

  else:
    # If no record standardiser for the data set is defined, the field names in
    # the input and the output data sets must be the same
    #
    input_field_name_list = input_dataset.fields.keys()
    input_field_name_list.sort()
    output_field_name_list = output_dataset.fields.keys()
    output_field_name_list.sort()

    if (input_field_name_list != output_field_name_list):
      print 'error:Field names differ in input and output data sets ' + \
            '(with no record standardiser defined)'
      raise Exception

  # Check if the data set defined in the blocking index (if defined) is correct
  #
  if (blocking_index != None):
    if (blocking_index.dataset != output_dataset):
      print 'error:Illegal data set definition in blocking index'
      raise Exception

  if (record_standardiser != None) and (blocking_index != None):
    do_string = 'Load, standardise and index records'
  elif (record_standardiser != None):
    do_string = 'Load and standardise records'
  elif (blocking_index != None):
    do_string = 'Load and index records'
  else:
    do_string = 'Load records'
  print '1:'
  print '1:  %s, write them into output data set' % (do_string)
  print '1:    First record: %i' % (first_record)
  print '1:    Last record:  %i' % (last_record-1)
  print '1:'

  start_time = time.time()  # Get current time
  comm_time  = 0.0          # Communication time

  input_rec_counter = first_record  # Current record pointer

  block_cnt = 0  # A round robin block counter, used for parallelism

  # Load records in a blocked fashion - - - - - - - - - - - - - - - - - - - - -

  while (input_rec_counter < last_record):

    block_size = min(febrl_block_size, (last_record - input_rec_counter))

    # Distribute blocks equally to all processors
    #
    if ((block_cnt % parallel.size()) == parallel.rank()):

      # Load original records from input data set
      #
      input_recs = input_dataset.read_records(input_rec_counter, block_size)
      print '1:    Loaded records %i to %i' % \
              (input_rec_counter, input_rec_counter+block_size)

      # Standardise them if a standardiser is defined - - - - - - - - - - - - -
      #
      if (record_standardiser != None):
        clean_recs = record_standardiser.standardise_block(input_recs)
        print '1:      Standardised records %i to %i' % \
              (input_rec_counter, input_rec_counter+block_size)
      else:
        clean_recs = input_recs  # Take the original records directly

      # Insert records into the blocking index if blocking index is defined - -
      #
      if (blocking_index != None):
        blocking_index.build(clean_recs)

      # If Febrl is run in parallel, send cleaned records to process 0
      #
      if (parallel.rank() > 0):
        tmp_time = time.time()
        parallel.send(clean_recs, 0)
        comm_time += (time.time() - tmp_time)

      else:  # Process 0, store standardised records

        # Store records in output data set
        #
        output_dataset.write_records(clean_recs)

    # If Febrl is run in parallel, process 0 receives cleaned records
    #
    if (parallel.rank() == 0) and (block_cnt % parallel.size() != 0):

      p = (block_cnt % parallel.size())  # Process number to receive from
      tmp_time = time.time()
      tmp_recs = parallel.receive(p)
      comm_time += (time.time() - tmp_time)

      # Store records in output data set
      #
      output_dataset.write_records(tmp_recs)

    input_rec_counter += block_size  # Increment current record pointer
    block_cnt += 1

    # Now determine timing and print progress report  - - - - - - - - - - - - -
    #
    if ((block_cnt % parallel.size()) == 0):
      used_time = time.time() - start_time
      recs_done = input_rec_counter - first_record
      perc_done = 100.0 * recs_done / number_records
      rec_time  = used_time / recs_done
      todo_time = (number_records - recs_done) * rec_time

      used_time_string = output.time_string(used_time)
      todo_time_string = output.time_string(todo_time)
      rec_time_string  = output.time_string(rec_time)

      print '1:      Processed %.1f%% of records in %s (%s per record)' % \
            (perc_done, used_time_string, rec_time_string)
      print '1:        Estimated %s until finished' % (todo_time_string)

  # End of standardisation  - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  parallel.Barrier()  # Make sure all processes are here

  total_time = time.time() - start_time  # Calculate total time
  total_time_string = output.time_string(total_time)
  print '1:  Total time needed for standardisation of %i records: %s' % \
        (number_records, total_time_string)

  if (parallel.size() > 1):
    comm_time_string = output.time_string(comm_time)
    print '1:    Time for communication: %s' % (comm_time_string)

  return [total_time, comm_time]

# =============================================================================

def do_comparison(dataset_a, dataset_b, record_comparator, classifier,
                  record_pair_dict, num_rec_pairs, febrl_block_size):
  """The main routine that does the comparison of record pairs given in the
     record pair list and the two data sets using the given record comparator.
     The resulting weight vectors are then inserted into the given classifier.
  """

  start_time =    time.time()  # Get current time
  compare_time =  0.0          # Comparison time
  classify_time = 0.0          # Classification time

  rec_pair_cnt = 0  # Loop counter

  # Compare records, and distribute comparisons equally to all processes  - - -
  #
  for rec_num_a in record_pair_dict:

    rec_num_a_dict = record_pair_dict[rec_num_a]

    rec_a = dataset_a.read_record(int(rec_num_a))  # Read the first record

    for rec_num_b in rec_num_a_dict:

      print '2:          Compare records %i with %i' % (rec_num_a, rec_num_b)

      # Read the records from the data set
      #
      rec_b = dataset_b.read_record(rec_num_b)

      # Compare the two records
      #
      tmp_time = time.time()
      weight_vector = record_comparator.compare(rec_a, rec_b)
      compare_time += (time.time() - tmp_time)

      # And insert the weight vector into the classifier
      #
      tmp_time = time.time()
      classifier.classify(weight_vector)  # Classify the weight vector
      classify_time += (time.time() - tmp_time)

      rec_pair_cnt += 1

      # Now determine timing and print progress report  - - - - - - - - - - - -
      #
      if ((rec_pair_cnt % febrl_block_size) == 0):
        used_time =       time.time() - start_time
        perc_done =       100.0 * rec_pair_cnt / num_rec_pairs
        rec_pair_time =   used_time / rec_pair_cnt
        todo_time =       (num_rec_pairs - rec_pair_cnt) * rec_pair_time
        avrg_comp_time =  (compare_time / rec_pair_cnt)
        avrg_class_time = (classify_time / rec_pair_cnt)

        used_time_string =       output.time_string(used_time)
        todo_time_string =       output.time_string(todo_time)
        rec_pair_time_string =   output.time_string(rec_pair_time)
        avrg_comp_time_string =  output.time_string(avrg_comp_time)
        avrg_class_time_string = output.time_string(avrg_class_time)

        print '1:      Processed %.1f%% (%i/%i) of record pairs in %s' % \
                (perc_done, rec_pair_cnt, num_rec_pairs, used_time_string) + \
                ' (%s per record pair)' % (rec_pair_time_string)
        print '1:        Average comparison time:     %s' % \
              (avrg_comp_time_string)
        print '1:        Average classification time: %s' % \
              (avrg_class_time_string)
        print '1:        Estimated %s until finished' % (todo_time_string)

  # Print final time for record pair comparison
  #
  total_time =    time.time() - start_time  # Calculate total time
  rec_pair_time = (total_time / rec_pair_cnt)

  total_time_string =    output.time_string(total_time)
  rec_pair_time_string = output.time_string(rec_pair_time)

  print '1:  Total time needed for comparison and classification of ' + \
        '%i record pairs: %s' % (rec_pair_cnt, total_time_string)
  print '1:    (%s per record pair)' % (rec_pair_time_string)

  return [total_time]

# =============================================================================