# =============================================================================
# indexing.py - Classes for indexing and blocking.
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
# The Original Software is "indexing.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module indexing.py - Classes for indexing and blocking.

   This module provides classes for building blocks and indices that will be
   used for the linkage process. Various derived classes are provided that
   implement different indexing techniques:
   - blocking index: The 'standard' blocking index used for record linkage.
   - windowing index: A sliding window over the sorted values of the index.
   - fuzzy bigram index: Allows fuzzy blocking with 'overlapping' blocks, like
     clustering (not fully implemented yet)

   The initialisation of the different index implementations and the building
   (i.e. filling the data structures) are all the same and commonly handled by
   the base class.

   For parallel runs, it is possible to merge several indexes into one (can be
   used to build one global index using several local indexes from different
   processes).

   After an index has been build (populated), it has to be 'compacted' into a
   data structure which allows efficient retrieval of the record numbers in the
   blocks. This compacting method depends on the index implementation.

   An iterator (based on a Python generator) can then be used to access all the
   blocks in an index.

   Two additional routines can be used to convert the blockss in one index (for
   deduplication) or two indexes (for linkage) into a list of record pairs that
   need to be compared using comparions record comparator functions.

   TODO:
   - implement routine 'linkage_rec_pairs'

"""

# =============================================================================

# Imports go here

from __future__ import generators
import encode
import parallel

# =============================================================================

def get_sublists(alist, length):
  """Routine to recursively compute all combinations of sub-lists of the list
     'alist' of length 'len'.

     Based on Gagan Saksena's routines from Activestate Python cookbook, see

     http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66465

     and modified by Ole Nielsen, MSI ANU, November 2002
  """

  sub_lists = []

  if (length == 0):
    sub_lists.append([])

  else:
    for i in range(len(alist)):
      sub = get_sublists(alist[i+1:], length-1)

      for l in sub:
        l.insert(0, alist[i])

      sub_lists += sub

  return sub_lists

# =============================================================================

def linkage_rec_pairs(index_a, index_b):
  """Routine to extract all record number pairs from the given indexing data
     structures that need to be compared for a record linkage process where two
     data sets are involved.

     This routine works for sequential and parallel runs, in that only a local
     part of all the record pairs are kept for the record pair dictionary.

     Returns a dictionary with record numbers as keys, and for each a
     dictionary with record numbers as keys (and values 1).
  """

  rec_pair_dict = {}
  duplicate_rec_pairs = 0  # Number of duplicate record pairs

  par_size = parallel.size()  # Get number of processes and my process number
  par_rank = parallel.rank()

  # Initialise the blocking iteration generator - - - - - - - - - - - - - - - -
  #
  block_iter_a = index_a.block_iterator()

  num_blocks_a = index_a.num_blocks  # Total number of blocks in index A
  num_blocks_b = index_b.num_blocks  # Total number of blocks in index B
  block_cnt =      0                      # Processed block counter 
  rec_pair_cnt =   0                  # Number of record pairs

  print '1:  Create record pair dictionary for linkage'

  # Loop over all blocks in the compacted index - - - - - - - - - - - - - - - -
  #
  blocking_value_a = block_iter_a.next()  # Get first block

  while blocking_value_a != None:  # Loop over all blocks

    block_var_a =            blocking_value_a[0]
    block_record_numbers_a = blocking_value_a[1]

    # Now get record numbers in this block from data set B
    #
    block_record_numbers_b = index_b.get_block_records(block_var_a)

    if (block_record_numbers_b != None):  # There are records in this block

      # Form cartesian product of all record pairs
      #
      for rec_num_a in block_record_numbers_a:

        if (rec_num_a % par_size == par_rank):  # A 'local' record number

          rec_num_a_dict = rec_pair_dict.get(rec_num_a, {})

          for rec_num_b in block_record_numbers_b:

            if (rec_num_a_dict.has_key(rec_num_b)):
              duplicate_rec_pairs += 1  # Record pair already in dictionary
            else:
              rec_num_a_dict[rec_num_b] = True
              rec_pair_cnt += 1  # One more record pairs

          if (rec_num_a_dict != {}):  # Only insert non-empty dictionary
            rec_pair_dict[rec_num_a] = rec_num_a_dict

    block_cnt += 1

    # Report progress every 10% (only if more than 10000 blocks)  - - - - - - -
    #
    if (num_blocks_a >= 10000) and (block_cnt % int(num_blocks_a / 10) == 0):
      print '1:      %i/%i blocks processed, number of record pairs: %i' % \
            (block_cnt, num_blocks_a, rec_pair_cnt)

    blocking_value_a = block_iter_a.next()  # Get next block

  print '1:    Number of record pairs:   %i' % (rec_pair_cnt)
  print '1:      Duplicate record pairs: %i (deleted)' % (duplicate_rec_pairs)

  return [rec_pair_dict, rec_pair_cnt]

# =============================================================================

def deduplication_rec_pairs(index):
  """Routine to extract all record number pairs from the given indexing data
     structure that need to be compared for a deduplication process where one
     data set is involved.

     This routine works for sequential and parallel runs, in that only a local
     part of all the record pairs are kept for the record pair dictionary.

     Returns a dictionary with record numbers as keys, and for each a
     dictionary with record numbers as keys (and values 1).
  """

  rec_pair_dict = {}
  duplicate_rec_pairs = 0  # Number of duplicate record pairs

  par_size = parallel.size()  # Get number of processes and my process number
  par_rank = parallel.rank()

  # Initialise the blocking iteration generator - - - - - - - - - - - - - - - -
  #
  block_iter = index.block_iterator()

  num_blocks =     index.num_blocks   # Total number of blocks in the index
  block_cnt =      0                  # Processed block counter 
  rec_pair_cnt =   0                  # Number of record pairs

  print '1:  Create record pair dictionary for deduplication'

  # Loop over all blocks in the compacted index - - - - - - - - - - - - - - - -
  #
  blocking_value = block_iter.next()  # Get first block

  while blocking_value != None:  # Loop over all blocks

    block_var =            blocking_value[0]
    block_record_numbers = blocking_value[1]

    # Do not process blocks with one record only
    #
    if (len(block_record_numbers) > 1):

      block_record_numbers.sort()  # Sort the record numbers in the block

      rec_cnt = 1  # Counter of record numbers within the block

      for rec_num_a in block_record_numbers:

        if (rec_num_a % par_size == par_rank):  # A 'local' record number

          rec_num_a_dict = rec_pair_dict.get(rec_num_a, {})

          for rec_num_b in block_record_numbers[rec_cnt:]:

            if (rec_num_a_dict.has_key(rec_num_b)):
              duplicate_rec_pairs += 1  # Record pair already in dictionary
            else:
              rec_num_a_dict[rec_num_b] = True
              rec_pair_cnt += 1  # One more record pairs

          if (rec_num_a_dict != {}):  # Only insert non-empty dictionary
            rec_pair_dict[rec_num_a] = rec_num_a_dict

        rec_cnt += 1

    block_cnt += 1

    # Report progress every 10% (only if more than 10000 blocks)  - - - - - - -
    #
    if (num_blocks >= 10000) and (block_cnt % int(num_blocks / 10) == 0):
      print '1:      %i/%i blocks processed, number of record pairs: %i' % \
            (block_cnt, num_blocks, rec_pair_cnt)

    blocking_value = block_iter.next()  # Get next block

  print '1:    Number of record pairs:   %i' % (rec_pair_cnt)
  print '1:      Duplicate record pairs: %i (deleted)' % (duplicate_rec_pairs)

  return [rec_pair_dict, rec_pair_cnt]

# =============================================================================

class Indexing:
  """Base class for indexing. Handles initialisation, building, and merging of
     indexes.

     Indexing is used to reduce the huge number of comparisons by grouping
     together records into 'blocks' that are similar according to some
     criterium. Indices are defined using a list of lists each containing
     tuples of the form (field_name, method, parameters). Each tuple will
     process the values in the given field (e.g. calculate the Soundex code,
     take the first three characters, etc.). If an index (i.e .a sub-list)
     contains more than one such tuples, these calculated values (strings) are
     concatenated to form a final blocking variable value.
     Each of the lists in an 'index_definition' defines one indexing criteria.

     For example,

       index_definition = [[('sname','soundex', 4, 'reverse')],
                           [('gname','truncate',2),('pcode', 'direct')],
                           [('pcode','truncate',2),('sname','nysiis')],
                          ]

     defines three indexes. The first will use the Soundex code of the
     reversed values in the field 'sname' with a maximal length of 4. The
     second index is the concatenation of the first two characters of the
     values in the field 'gname' plus the values in the field 'pcode' (taken
     directly from the field without any encoding or truncation). The third
     index is made of the first two characters from the values in the field
     'pcode' concatenated with the NYSIIS code of the values in the field
     'sname' (with the default code length).

     For each index, a dictionary will be built, with the values of the
     blocking variables as keys, and the corresponding record numbers stored
     in dictionaries.

     The possible methods for indexing are the encodings 'soundex',
     'mod_soundex', 'phonex', 'nysiis', and 'dmetaphone', plus the 'truncate'
     and the 'direct' methods.

     For all encoding methods a first parameter is the maximal length of the
     code (a positive integer number). If no length is given the default as
     set in module 'encode.py' is used. A second parameter can be the string
     'reverse' in which case the reversed string is encoded.

     For the 'truncate' method the additional parameter is an positive integer
     number, which is the position (starting with zero) where strings are
     truncated.

     No additional parameter is needed for the 'direct' method.

     The following arguments can be given to the constructor of the base class
     and all derived classes:
       name
       description
       dataset           The data set the index is built on
       index_definition  Definitions of the indexes as described above
       skip_missing      A flag, if set to 'True' records which have empty
                         indexing values will be skipped over, if set to
                         'False' a block with an empty indexing value will be
                         included in the index as well.
                         Default value is 'True'.

     Both the data set and the index must be defined at initialisation time.

     The following methods are available in all derived classes:
       __init__           Initialise a indexing object
       build              Build an index with the given records (updates an
                          index when it is not empty).
       merge              Merge an index data structure with another index of
                          the same derived class. This can only be done before
                          an index is compacted.
       compact            After an index has been built, this method can be
                          used to compact an index into a more efficient form
                          for faster retrieval of the blocks it contains.
       block_iterator     An iterator over all blocks in an index. Returns the
                          value of the blocking variable as a string plus a
                          list with the record numbers in this block. If no
                          more blocks are available, None is returned. Uses a
                          Python generator.
       get_block_records  For a given blocking variable value returns a list
                          with record numbers in the corresponding block, or
                          None if this block is not available in the index.

     Before 'block_iterator' method can be used, an index must be compacted
     using the 'compact' method.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all index implementations
    #
    self.name =             ''
    self.description =      ''
    self.dataset =          None
    self.index_definition = None
    self.num_indexes = 0          # Number of indexes in index definition
    self.field_names = []         # A list with field names for each index
    self.skip_missing =     True
    self.index = []               # A list for the block index dictionaries
    self.compact_index = None     # A compacted version of the index for
                                  # faster access

    self.num_blocks = 0           # Total number of blocks in the index
    self.compacted = False        # A flag, set to True once 'compact' method
                                  # has been called

    # Define a list with all possible indexing methods
    #
    self.methods = ['soundex','mod_soundex','phonex','nysiis','dmetaphone', \
                    'truncate','direct']

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'dataset'):
        self.dataset = value

      elif (keyword in ['skip','skip_missing']):
        if (value not in [True, False]):
          print 'error:Argument "skip_missing" must be "True" or "False"'
          raise Exception
        self.skip_missing = value

      elif (keyword in ['index_def', 'index_definition']):
        if (not isinstance(value, list)):
          print'error:Argument "index_definition" is not a list'
          raise Exception
        self.index_definition = value

      else:
        print 'error:Illegal constructor argument keyword: '+keyword
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.dataset == None):
      print 'error:Data set is not defined'
      raise Exception

    if (self.index_definition == None):
      print 'error:Index is not defined'
      raise Exception

    # Check if definition of indexes is correct and fields are valid in - - - -
    # the given data set
    #
    for index in self.index_definition:
      if (not isinstance(index, list)):
        print 'error:Index "%s" in "index_definition" is not a list' % \
              (str(index))
        raise Exception
      for field_tuple in index:
        if (not isinstance(field_tuple, tuple)):
          print 'error:Tuple "%s" in index "%s" is not a tuple' % \
                (str(field_tuple), str(index))
          raise Exception

        if (len(field_tuple) < 2):
          print 'error:Length of tuple "%s" in index "%s" is less than two' % \
                (str(field_tuple), str(index))
          raise Exception

        if (not self.dataset.fields.has_key(field_tuple[0])):
          print 'error:Field "%s" not in data set %s' % \
                (str(field_tuple[0]), str(self.dataset.name))
          raise Exception

    # Process index definitions - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for index in self.index_definition:
      self.index.append({}) # A new empty dictionary for this block index

      field_list = []

      for field_tuple in index:
        field_list.append(field_tuple[0])  # Append the field name

        if (field_tuple[1] not in self.methods):
          print 'error:Illegal method for field "%s" in index: "%s": %s' % \
                (str(field_tuple[0]), str(index), str(field_tuple[1]))
          raise Exception

        if (field_tuple[1] == 'truncate') and (len(field_tuple) != 3):
          print 'error:Missing "length" with method "truncate" in ' + \
                'index: "%s"' % (str(index))
          raise Exception

        elif (field_tuple[1] in ['soundex','mod_soundex','phonex','nysiis', \
                                 'dmetaphone']) and (len(field_tuple) >= 3):
          if (not isinstance(field_tuple[2], int)) or (field_tuple[2] <= 0):
            print 'error:Length given to encoding method is not a positive '+ \
                  'integer'
            raise Exception

        if (field_tuple[1] in ['soundex','mod_soundex','phonex','nysiis', \
                               'dmetaphone']) and (len(field_tuple) == 4):
          if (field_tuple[3] not in ['reverse','rev']):
            print 'error:Second parameter for tuple "%s" is not "reverse"' % \
                  (str(field_tuple))
            raise Exception

        if (len(field_tuple) > 4):
          print 'error:Field tuple has to many elements: %s' % \
                (str(field_tuple))
          raise Exception

      self.field_names.append(field_list)

    if (len(self.index) != len(self.field_names)):
      print 'error:Illegal number of fields or indexes'
      raise Exception

    self.num_indexes = len(self.index)

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised blocking index:'
    print '1:  Name:                %s' % (str(self.name))
    print '1:  Data set:            %s' % (str(self.dataset.name))
    print '1:  Skip missing values: %s' % (str(self.skip_missing))
    print '1:  Used field names:    %s' % (str(self.field_names))
    print '1:  Number of indexes:   %i' % (self.num_indexes)
    print '2:  Indexes:'
    for i in self.index_definition:
      print '2:      %s ' %(str(i))

  # ---------------------------------------------------------------------------

  def build(self, record_list):
    """Build one or more blocking indexes.
       The indexes can already be populated from an earlier call to 'build'.

       The given records are processed, with the indexing variables values and
       record numbers being inserted into the index.
    """

    # Check that index has not been compacted yet - - - - - - - - - - - - - - -
    #
    if (self.compacted == True):
      print 'error:Index is already compacted, building not possible'
      raise Exception

    if (not isinstance(record_list, list)):
      print 'error:Argument "record_list" is not a list'
      raise Exception

    num_empty_block_var = 0  # Number of records with empty block variable
    num_rec_inserted =    0  # Number of records inserted into the index

    # Loop over all records in the list - - - - - - - - - - - - - - - - - - - -
    #
    for record in record_list:

      # Make sure record is a dictionary
      #
      if (not isinstance(record, dict)):
        print 'error:Illegal record type: "%s", must be a dictionary' % \
              (str(type(record)))
        raise Exception

      record_id = '[RecID: %s/%s]' % \
                  (str(record['_rec_num_']),record['_dataset_name_'])
      print '3:%s  Record: %s' % (record_id,str(record))

      rec_num = record['_rec_num_']  # Get number of the current record

      # Extract needed fields and process them for each block index
      #
      for i in range(self.num_indexes):

        j = 0  # Index into current fields and functions
        block_var = ''  # A new value for the block variable

        for field_name in self.field_names[i]:  # Get the name of the field

          # Get the value from the record, or an empty string if not available
          #
          field_value = str(record.get(field_name,''))

          if (field_value != ''):

            # Now process methods and set up fields and functions lists - - - -
            #
            if (self.index_definition[i][j][1] == 'direct'):
              block_var += field_value

            elif (self.index_definition[i][j][1] == 'truncate'):
              block_var += field_value[:self.index_definition[i][j][2]]

            else:  # Must be an encoding method

              if (len(self.index_definition[i][j]) == 4) and \
                 (self.index_definition[i][j][3] in ['rev','reverse']):

                reverse_list = list(field_value)  # Reverse value
                reverse_list.reverse()
                field_value = ''.join(reverse_list)

              if (self.index_definition[i][j][1] == 'soundex'):
                if (len(self.index_definition[i][j]) == 2):
                  block_var += encode.soundex(field_value)
                else:  # A length must have been given
                  block_var += encode.soundex(field_value, 
                                              self.index_definition[i][j][2])

              elif (self.index_definition[i][j][1] == 'mod_soundex'):
                if (len(self.index_definition[i][j]) == 2):
                  block_var += encode.mod_soundex(field_value)
                else:
                  block_var += encode.mod_soundex(field_value,
                                               self.index_definition[i][j][2])

              elif (self.index_definition[i][j][1] == 'phonex'):
                if (len(self.index_definition[i][j]) == 2):
                  block_var += encode.phonex(field_value)
                else:
                  block_var += encode.phonex(field_value,
                                             self.index_definition[i][j][2])

              elif (self.index_definition[i][j][1] == 'nysiis'):
                if (len(self.index_definition[i][j]) == 2):
                  block_var += encode.nysiis(field_value)
                else:
                  block_var += encode.nysiis(field_value,
                                             self.index_definition[i][j][2])

              elif (self.index_definition[i][j][1] == 'dmetaphone'):
                if (len(self.index_definition[i][j]) == 2):
                  block_var += encode.dmetaphone(field_value)
                else:
                  block_var += encode.dmetaphone(field_value,
                                               self.index_definition[i][j][2])
              else:
                print 'error:Illegal encoding method: %s' % \
                      (str(self.index_definition[i][j][1]))
                raise Exception

          j += 1

        if (block_var == ''):
          num_empty_block_var += 1

        # Insert record number into current index - - - - - - - - - - - - - - -
        # (if 'skip_missing' flag is set to True don't insert if block variable
        # is empty)
        if (block_var != '') or (self.skip_missing == False):

          block_var = str(i)+':'+block_var  # Insert index number at beginning
          index_val = self.index[i].get(block_var,{})
          index_val[rec_num] = 1
          self.index[i][block_var] = index_val
          num_rec_inserted += 1

          print '3:%s    Inserted record number %s with index value ' % \
                (record_id, str(rec_num))+ \
                '"%s" into index %i' % (str(block_var), i)

    print '1:      Inserted %i record numbers into index' % (num_rec_inserted)
    if (self.skip_missing == False):
      print '1:        Number of empty blocking values: %i' % \
            (num_empty_block_var) + ' (which were inserted)'
    else:
      print '1:        Number of empty blocking values: %i' % \
            (num_empty_block_var) + ' (which were not inserted)'

    # Count and print the total number of blocks  - - - - - - - - - - - - - - -
    #
    self.num_blocks = 0
    for i in range(self.num_indexes):
      self.num_blocks += len(self.index[i])

    print '1:    Total number of blocks in index: %i' % (self.num_blocks)

  # ---------------------------------------------------------------------------

  def merge(self, other_index):
    """Merge an indexing data structure with another indexing data structure.

       The 'other_index' must only contain the list of indexes from another
       index, but not the complete index object.

       The number of indexes (indexes are a list of index dictionaries) must
       be the same in self.index and other_index.

       Each of the index dictionaries are then merged.
    """

    # Check if this index is not compacted (unfortunatey we can not check this
    # for the other index)
    #
    if (self.compacted == True):
      print 'error:Cannot merge compacted indexes'
      raise Exception

    # Check if both indexes have the same numbe rof index dictionaries
    #
    if (len(other_index) != self.num_indexes):
      print 'error:Different number of index dictionaries in merge routine' + \
            ': self: %i, other: %i' % (self.num_indexes, len(other_index))

    for i in range(self.num_indexes):

      other_index_dict = other_index[i]

      # Get keys (block values) from other index
      #
      for block_value in other_index_dict:

        # Get record numbers with this block value
        #
        rec_dict = other_index_dict[block_value]

        index_val = self.index[i].get(block_value,{})
        index_val.update(rec_dict)  # Merge dictionaries for this value

        self.index[i][block_value] = index_val  # Store back into this index

    # Count and print the total number of blocks  - - - - - - - - - - - - - - -
    #
    self.num_blocks = 0
    for i in range(self.num_indexes):
      self.num_blocks += len(self.index[i])

    print '1:    Merged indexes (new total number of blocks in index: %i)' % \
          (self.num_blocks)

  # ---------------------------------------------------------------------------

  def compact(self):
    """Build a compacted version of the index.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

  # ---------------------------------------------------------------------------

  def block_iterator(self):
    """An iterator over all blocks in an index. Returns the numbers of the
       records in a block as a list. If no more block is available, None is
       returned. Uses a Python generator.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

  # ---------------------------------------------------------------------------

  def get_block_records(self, block_var):
    """For a given blocking variable value returns a list with record numbers
       in the corresponding block, or None if this block is not available in
       the index.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

# =============================================================================

class BlockingIndex(Indexing):
  """Class that implements the 'classical' blocking used for record linkage.

     When compacting the index, nothing needs to be done, as index is already
     stored in dictionaries and ready to be retrieved block wise.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Simply call the base class constructor.
    """

    Indexing.__init__(self, kwargs)  # Initialise base class

    print '1:  Indexing type:       Classical blocking index'

  # ---------------------------------------------------------------------------

  def build(self, record_list):
    """Build the index, insert the given records. Is done by the base class.
    """

    Indexing.build(self, record_list)

  # ---------------------------------------------------------------------------

  def merge(self, other_index):
    """Merge an indexing data structure with another indexing data structure.
       Is done by the base class.
    """

    Indexing.merge(self, other_index)

  # ---------------------------------------------------------------------------

  def compact(self):
    """Build a compacted version of the blocking indexes.

       Nothing needs to be done, as index is already stored in dictionaries and
       ready to be retrieved block wise.
    """

    self.compact_index = self.index
    self.compacted = True

    print '1:    Index compacted (nothing to be done)'
    print '1:      Total number of blocks in index: %i' % (self.num_blocks)

  # ---------------------------------------------------------------------------

  def block_iterator(self):
    """An iterator over all blocks in the index.
    """

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.compacted == False):
      print 'error:Index has not been compacted, block iteration not possible'
      raise Exception

    # Outer loop is over the index dictionaries - - - - - - - - - - - - - - - -
    #
    for i in range(self.num_indexes):

      for (block_var, rec_dict) in self.compact_index[i].items():

        rec_list = rec_dict.keys()  # Get the record numbers for in this block

        yield [block_var, rec_list]  # Return blocking variable and record list

    yield None  # End of iteration, no more blocks available

  # ---------------------------------------------------------------------------

  def get_block_records(self, block_var):
    """For a given blocking variable value returns a list with record numbers
       in the corresponding block, or None if this block is not available in
       the index.
    """

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.compacted == False):
      print 'error:Index has not been compacted, block iteration not possible'
      raise Exception

    index_num = int(block_var[0])  # Get the number of the index

    if self.compact_index[index_num].has_key(block_var):
      rec_dict= self.compact_index[index_num][block_var]  # Found

      return rec_dict.keys()

    return None

# =============================================================================

class SortingIndex(Indexing):
  """Class that implements a sorting index and uses a moving window over the
     sorted indexes.

     When a sorting index is initialised, one argument (besides the base class
     arguments) that needs to be given is:

     window_size  A positive integer that gives the size of the sliding window
                  in records
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the window size argument, then simply call the
       base class constructor.
    """

    self.window_size = None  # Set the window size to not defined

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['window','window_size']):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Window size is not a positive integer'
          raise Exception
        self.window_size = value

      else:
        base_kwargs[keyword] = value

    if (self.window_size == None):
      print 'error:Window size is not defined'
      raise Exception

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    print '1:  Window size:         %i' % (self.window_size)
    print '1:  Indexing type:       Sliding window index'

  # ---------------------------------------------------------------------------

  def build(self, record_list):
    """Build the index, insert the given records. Is done by the base class.
    """

    Indexing.build(self, record_list)

  # ---------------------------------------------------------------------------

  def merge(self, other_index):
    """Merge an indexing data structure with another indexing data structure.
       Is done by the base class.
    """

    Indexing.merge(self, other_index)

  # ---------------------------------------------------------------------------

  def compact(self):
    """Build a compacted version of the sorting indexes.

       Convert each of the indexes into a sorted (according to the key values)
       list.
    """

    self.compact_index = []  # Make a list of lists for sorted indexes

    # Loop over all indexes
    #
    for i in range(self.num_indexes):

      # Sort according to keys (field values)
      #
      block_keys = self.index[i].keys()
      block_vals = self.index[i].values()
      block_list = map(None, block_keys, block_vals)
      block_list.sort()

      self.compact_index.append(block_list)

    self.compacted = True

    print '1:    Index compacted (sorted indexes)'
    print '1:      Total number of blocks in index: %i' % (self.num_blocks)

  # ---------------------------------------------------------------------------

  def block_iterator(self):
    """An iterator over all blocks in the index.
    """

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.compacted == False):
      print 'error:Index has not been compacted, block iteration not possible'
      raise Exception

    # Loop over all block indexes
    #
    for i in range(self.num_indexes):
      index_len =   len(self.index[i])
      block_start = 0
      block_end =   min(block_start+self.window_size, index_len)

      while (block_start < index_len):

        block = self.compact_index[i][block_start:block_end]

        rec_nums = {}
        for (key, rec_list) in block:
          rec_nums.update(rec_list)

        block_var =  block[0][0]+'...'+block[-1][0]  # Make a 'name' for the
                                                     # blocking variable
        block_recs = rec_nums.keys()

        yield [block_var, block_recs]

        block_start += 1
        block_end =    min(block_start+self.window_size, index_len)

    yield None

  # ---------------------------------------------------------------------------

  def get_block_records(self, block_var):
    """For a given blocking variable value returns a list with record numbers
       in the corresponding block, or None if this block is not available in
       the index.
       Returns a 'block' where the given value is the first ...

    """

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.compacted == False):
      print 'error:Index has not been compacted, block querying not possible'
      raise Exception

    index_num = int(block_var[0])  # Get the number of the index

    [block_var1, block_var2] = block_var.split('...')

    block_var1 = block_var1[:-1]

    # Find the entry in the index that starts with either 'block_var1' or is
    # is the one before (alphabetically)
    # Use binary search as the index is a sorted list
    #
    var_len = len(block_var1)
    start = 0
    end = len(self.compact_index[index_num]) - 1
    found = -1

    while (found == -1):
      if (end < start):
        found = middle
      middle = (start + end) / 2
      middle_index_val = self.compact_index[index_num][middle]

      if (middle_index_val[0][:var_len] < block_var1):
        start = middle+1
      elif (middle_index_val[0][:var_len] > block_var1):
        end = middle-1
      else:
        found = middle

    #if self.index[index_num].has_key(block_var):
    #  rec_dict= self.index[index_num][block_var]  # Found
    #
    #  return rec_dict.keys()

    return None

  # ===========================================================================

class BigramIndex(Indexing):
  """Class that implements an indexing structure based on Bigrams and allows
     for fuzzy 'blocking'.

     When the index is compacted, a value in a blocking variable will be
     converted into a list of bigrams and permutations of sub-lists will be
     built using the given threshold (a number between 0.0 and 1.0) of all
     possible permutations. The resulting bigram lists will be inserted into
     an inverted index, i.e. the record number will be inserted into
     dictionaries for each bigram. This inverted index is then the compacted
     index which will be used to retrieve the blocks.

     When a bigram index is initialised, one argument (besides the base class
     arguments) that need to be given is:

     threshold  A number between 0.0 (not included) and 1.0

     For example, assume a block definition contains the tuple:

       block_definition = [[('sname','direct')], ...]

     and the bigram threshold is set to 0.8. If a value 'peter' is given in an
     'sname' field, the corresponding bigram list is then ['pe','et','te','er']
     with four elements, so using the 0.8 threshold results in 4*0.8 = 3.2
     rounded to 3, which means all permutations of length 3 are calculated. For
     the given example they are ['pe','et','te'], ['pe','et','er'],
     ['pe','te','er'] and ['et','te','er']. So, the record number of this
     example will be inserted into the inverted index blocks with keys
     'peette', 'peeter', 'peteer, and 'etteer'.

     The lower the threshold, the shorter the sub-lists, but also the more
     sub-lists there will be per field value, resulting in more (smaller
     blocks) in the inverted index.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the threshold argument, then simply call the
       base class constructor.
    """

    self.threshold = None  # Set the threshold to not defined

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['thres','threshold']):
        if (not isinstance(value, float)) or (value <= 0.0) or \
           (value > 1.0):
          print 'error:Threshold is not a valid floating-point number '+ \
                '(must be between 0.0 and 1.0): %s' % (str(value))
          raise Exception
        self.threshold = value

      else:
        base_kwargs[keyword] = value

    if (self.threshold == None):
      print 'error:Threshold is not defined'
      raise Exception

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    print '1:  Threshold:           %f' % (self.threshold)
    print '1:  Indexing type:       Bigram index'

  # ---------------------------------------------------------------------------

  def build(self, record_list):
    """Build the index, insert the given records. Is done by the base class.
    """

    Indexing.build(self, record_list)

  # ---------------------------------------------------------------------------

  def merge(self, other_index):
    """Merge an indexing data structure with another indexing data structure.
       Is done by the base class.
    """

    Indexing.merge(self, other_index)

  # ---------------------------------------------------------------------------

  def compact(self):
    """Build a compacted version of the bigram indexes.

       For each of the blocks in the original index, the bigram permutations
       are constructed and the record numbers in the block are inserted into
       the corresponding blocks in the inverted index.
    """

    self.compact_index = []  # Make a list of lists for sorted indexes
    self.num_blocks = 0

    # Loop over all indexes
    #
    for i in range(self.num_indexes):

      inv_index = {}  # Start a new inverted index

      # Loop over all blocks in this index
      #
      for (block_var, rec_dict) in self.index[i].items():

        block_var_ind = block_var[:2]  # Get the index number and ':'
        block_var_val = block_var[2:]

        # Create the bigram list and all sublists of this blocking value
        #
        if (block_var_val != ''):

          field_bigram = []
          for j in range(len(block_var_val)-1):
            # Don't insert duplicates (discussion with Tim Churches, 21/2/2003)
            #
            if (block_var_val[j:j+2] not in field_bigram):
              field_bigram.append(block_var_val[j:j+2])

          field_bigram.sort()  # Sort the bigram list (discussion with Tim)

          # Compute length of sub-lists needed from given 'threshold'
          #
          num_value_bigram = len(field_bigram)

          bigrams_needed = int(round(float(num_value_bigram) * self.threshold))
          if (bigrams_needed < 1):
            bigrams_needed = 1

          # Form all sub-sets (combinations) with bigrams
          #
          bigram_sublists = get_sublists(field_bigram, bigrams_needed)

          # Convert bigram sublists into strings and insert into compact index
          #
          for bg in bigram_sublists:
            bg_string = block_var_ind+''.join(bg)

            block_dict = inv_index.get(bg_string,{})
            block_dict.update(rec_dict)
            inv_index[bg_string] = block_dict

      self.compact_index.append(inv_index)
      self.num_blocks += len(inv_index)

    self.compacted = True

    print '1:    Index compacted (bigram indexes)'
    print '1:      Total number of blocks in index: %i' % (self.num_blocks)

  # ---------------------------------------------------------------------------

  def block_iterator(self):
    """An iterator over all blocks in the index.
    """

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.compacted == False):
      print 'error:Index has not been compacted, block iteration not possible'
      raise Exception

    # Outer loop is over the index dictionaries - - - - - - - - - - - - - - - -
    #
    for i in range(self.num_indexes):

      for (block_var, rec_dict) in self.compact_index[i].items():

        rec_list = rec_dict.keys()  # Get the record numbers for in this block

        yield [block_var, rec_list]  # Return blocking variable and record list

    yield None  # End of iteration, no more blocks available

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

#### PC, 21/2/2003 (old code below)

  def num_blocks(self):
    """Returns the number of blocks in the index.
    """

    if (self.index == {}):
      return None

    else:
      return len(self.index)

# =============================================================================
# =============================================================================

  def get_similar(self, value_dict, perc_common):
    """Get record numbers with similar field values as the given argument.
       'value_dict' is a dictionary where the keys are field names and the
       corresponding values are used to search for similar values in the
       index.

       'perc_common' gives the minimum percentage that values have to be
       similar in order to be included into the similar record list.

       This routine returns a list of record numbers.
    """

    assert (isinstance(perc_common,float)) and \
           (perc_common > 0.0) and (perc_common < 1.0), \
           'Argument "perc_common" must be a number between 0.0 and 1.0'

    start_done = False  # Set to 'True' once the first fields is processed

    # Loop over all field/value pairs - - - - - - - - - - - - - - - - - - - - -
    #
    for (field,value) in value_dict.items():

      assert self.dataset.fields.has_key(field), \
               'Field "+str(field)+" not in given data set'

      # Make a bigram list for the current value
      #
      value_bigram = []
      for i in range(len(value)-1):
        value_bigram.append(value[i:i+2])

      # Compute length of sub-lists needed from given 'perc_common'
      #
      num_value_bigram = len(value_bigram)
      bigrams_needed = int(float(num_value_bigram) * perc_common)

      if (bigrams_needed < 1):
        bigrams_needed = 1

      # Form all sub-sets (combinations) with bigrams
      #
      bigram_sublists = self.get_sublists(value_bigram, bigrams_needed)

      field_rec_dict = {}  # A dictionary that will contain the union of all
                           # record numbers in the intersection lists

      # For each sub-list with bigrams get corresponding dictionaries
      #
      for bigram_list in bigram_sublists:

        bigram_dict_list = []  # A list with dictionaries of bigrams in a list

        for b in bigram_list:  # For each bigram in this sub-list

          if (self.indices[field].has_key(b)):
            bigram_dict_list.append(self.indices[field][b])
          else:
            bigram_dict_list.append({}) # Empty dictionary if a bigram is not
                                        # in the dictionary for this field

        if ({} not in bigram_dict_list):  # If all bigrams were in dictionary

          intersect_list = bigram_dict_list[0].keys()
          for d in bigram_dict_list[1:]:
            intersect_list = filter(d.has_key,intersect_list)

        else:
          intersect_list = []  # No records for this bigram sub-list

        # Now insert the intersection list into record dictionary (union) - - -
        #
        for rec_num in intersect_list:
          field_rec_dict[rec_num] = 1

      a = field_rec_dict.keys()
      #a.sort()
      #print field, a
      #print

      # Do intersection of all field record lists into final list - - - - - - -
      #
      if (not start_done):
        start_done = True
        similar_recs = field_rec_dict.keys()

      else:  # Make intersection
        similar_recs = filter(field_rec_dict.has_key, similar_recs)

    return similar_recs

  # ---------------------------------------------------------------------------

# =============================================================================
# =============================================================================
