# =============================================================================
# comparison.py - Classes for field and record comparions.
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
# The Original Software is "comparison.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
# Daniel R. Sabath M.Sc. (Research Consultant, Harborview Injury Prevention and
#                         Research Center, Seattle)
# =============================================================================

"""Module comparison.py - Classes for field and record comparions.

   This module provides classes for field and record comparisons which are
   used for the linkage process.

   The FieldComparator base class contains different field comparison derived
   classes:
     FieldComparatorExactString     Exact string comparison
     FieldComparatorTruncateString  Exact string comparison with truncated
                                    strings
     FieldComparatorApproxString    Approximate string comparison
     FieldComparatorEncodeString    Phonetic string encoding comparison
     FieldComparatorKeyDiff         Keying error comparison
     FieldComparatorNumericPerc     Comparison of numeric values with
                                    percentage tolerance
     FieldComparatorNumericAbs      Comparison of numeric values with
                                    absolute difference tolerance
     FieldComparatorDate            Date comparison with day difference
                                    tolerance
     FieldComparatorAge             Age comparison with percentage tolerance
     FieldComparatorTime            Time comparison with minute tolerance
     FieldComparatorDistance        Distance comparison with kilometer
                                    tolerance

   TODO
   - The age comparator might alternatively take two fields, check if they are
     integer numbers (ages).
   - What do we do if we want to link two data sets, and we know that we have
     different distributions for e.g. surnames, but only a lookup table for the
     surnames from one data set?
   - Add frequency table capabilities for numeric, data and age field
     comparators.
   - Do a proper date/time comparison function (maybe using mxDateTime module)
     that takes both date and time into consideration.
"""

# =============================================================================
# Imports go here

import math
import time

import mymath     # Module with mathematical routines (e.g. log2)
import stringcmp  # Module with algorithms for approximate string comparison
import encode     # Module with string encodings (Soundex, NYSIIS, etc.)
import date       # Module with routines for dates

# =============================================================================

class RecordComparator:
  """Class that implements record comparator to compare two records and compute
     a weight vector.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, dataset_a, dataset_b, field_comparators):
    """Constructor.
    """

    self.dataset_a = dataset_a
    self.dataset_b = dataset_b

    if (not isinstance(field_comparators, list)):
      print 'error:Argument "field_comparators" is not a list'
      raise Excpetion
    self.field_comparators = field_comparators

    # Set the data set attributes in the field comparators and check  - - - - -
    # if the fields are proper fields of the data sets
    #
    for fc in self.field_comparators:
      fc.dataset_a = self.dataset_a
      fc.dataset_b = self.dataset_b

      for field in fc.fields_a:
        if (not self.dataset_a.fields.has_key(field)):
          print 'error:Data set A does not have field: "%s"' % (str(field))
          raise Exception
      for field in fc.fields_b:
        if (not self.dataset_b.fields.has_key(field)):
          print 'error:Data set B does not have field: "%s"' % (str(field))
          raise Exception

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised record comparator'
    print '1:  Data set A: %s' % (str(self.dataset_a.name))
    print '1:  Data set B: %s' % (str(self.dataset_b.name))
    print '1:  Field comparators:'
    for fc in self.field_comparators:
      print '1:    Name: %s' % (str(fc.name))
      print '1:      Fields A: %s' % (str(fc.fields_a))
      print '1:      Fields B: %s' % (str(fc.fields_b))

  # ---------------------------------------------------------------------------

  def compare(self, rec_a, rec_b):
    """Compare two records (dictionaries for fields) and return vector with
       weight values (floating-point numbers)
    """

    records_id = '[RecID A: %s/%s, RecID B: %s/%s]' % \
                 (str(rec_a['_rec_num_']), str(rec_a['_dataset_name_']), \
                  str(rec_b['_rec_num_']), str(rec_b['_dataset_name_']))

    # Initialise weight vector with record identifier values at beginning - - -
    #
    weight_vector = [rec_a['_dataset_name_'], rec_a['_rec_num_'], \
                     rec_b['_dataset_name_'], rec_b['_rec_num_']]

    # Compute a weight for each field comparator
    #
    for fc in self.field_comparators:

      # Extract field values from data sets
      #
      fields_a = []
      fields_b = []
      for field in fc.fields_a:
        fields_a.append(rec_a.get(field, fc.dataset_a.fields_default))
      for field in fc.fields_b:
        fields_b.append(rec_b.get(field, fc.dataset_b.fields_default))

      w = fc.compare(fields_a, fields_b, records_id)
      weight_vector.append(w)

    # A log message for medium volume log output (level 2)  - - - - - - - - - -
    #
    print '2:  Record comparison:'
    print '2:    %s, comparison vector: %s' % (records_id, str(weight_vector))

    return weight_vector

  # ---------------------------------------------------------------------------

  def compare_block(self, rec_list_a, rec_list_b):
    """Method to compare each record in list A with all records in list B.

       If the records in list A are the same as the records in list B (e.g. for
       a deduplication task), then the number of comparisons is reduced by half
       by only comparing records in the list B that are 'after' a record.
    """

    if (not isinstance(rec_list_a, list)) and \
       (not isinstance(rec_list_a, list)):
      print 'error:Both record blocks must be of type list'
      raise Exception

    if (rec_list_a == rec_list_b):
      same_lists = True
    else:
      same_lists = False

    weights = []
    rec_cnt = 1
    for rec_a in rec_list_a:

      if (same_lists == False):  # Compare all records
        for rec_b in rec_list_b:
          weights.append(self.compare(rec_a, rec_b))

      else:  # Only do comparisons not done yet
        for rec_b in rec_list_b[rec_cnt:]:
          weights.append(self.compare(rec_a, rec_b))

      rec_cnt += 1

    return weights

# =============================================================================

class FieldComparator:
  """Base class for field comparators.

     The following arguments must be given to the constructor of the base class
     and all derived classes:
       fields_a
       fields_b
       missing_weight
       m_probability *
       u_probability *

     Agreement and disagreement weights are computed using the m- and
     u-probabilities.

     Two data set attributes will be set by the record comparator so that field
     comparators have access to the missing value attribute of the data sets. 

     * For the "FieldComparatorDate" separate m- and u-probabilities need to be
       given for day, month and year.

     The following optional arguments can be used with most field comparators:
       name
       description
       frequency_table
       freq_table_max_weight
       freq_table_min_weight

     Default value for 'missing_weight' is zero.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all data sets
    #
    self.name =                  ''
    self.description =           ''
    self.dataset_a =             None  # Will be set by record comparator
    self.dataset_b =             None  # Will be set by record comparator
    self.fields_a =              None
    self.fields_b =              None
    self.m_probability =         None
    self.u_probability =         None
    self.agree_weight =          None
    self.disagree_weight =       None
    self.missing_weight =        0.0
    self.frequency_table =       None
    self.freq_table_max_weight = None
    self.freq_table_min_weight = None

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'fields_a'):
        if (isinstance(value, str)):
          self.fields_a = [value]
        elif (isinstance(value, list)):
          self.fields_a = value
        else:
          print 'error:Argument "fields_a" must be of type string or list'
          raise Exception

      elif (keyword == 'fields_b'):
        if (isinstance(value, str)):
          self.fields_b = [value]
        elif (isinstance(value, list)):
          self.fields_b = value
        else:
          print 'error:Argument "fields_b" must be of type string or list'
          raise Exception

      elif (keyword in ['m', 'm_prob', 'm_probability']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "m_probability" is not a number'
          raise Exception

        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability" value: %s' % (str(value))
          raise Exception

        self.m_probability = value

      elif (keyword in ['u', 'u_prob', 'u_probability']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "u_probability" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability" value: %s' % (str(value))
          raise Exception
        self.u_probability = value

      elif (keyword in ['miss', 'miss_w', 'missing_w', 'missing_weight']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "missing_weight" is not a number'
          raise Exception
        self.missing_weight = value

      elif (keyword == 'frequency_table'):
        self.frequency_table = value

      elif (keyword in ['freq_table_max_w', 'freq_table_max_weight']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "freq_table_max_weight" is not a number'
          raise Exception
        self.freq_table_max_weight = value

      elif (keyword in ['freq_table_min_w', 'freq_table_min_weight']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "freq_table_min_weight" is not a number'
          raise Exception
        self.freq_table_min_weight = value

      else:
        print 'error:Illegal constructor argument keyword: %s' % (str(keyword))
        raise Exception

    # Check if fields are not empty - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.fields_a == None) or (self.fields_b == None):
      print 'error:Fields are not defined'
      raise Exception

    if (len(self.fields_a) == 0) or (len(self.fields_b) == 0):
      print 'error:Fields must not be empty lists'
      raise Exception

    # Compute weights if probabilities are given  - - - - - - - - - - - - - - -
    #
    if (self.m_probability != None) and (self.u_probability != None):
      self.agree_weight = mymath.log2(self.m_probability / self.u_probability)
      self.disagree_weight = mymath.log2((1.0-self.m_probability) /
                                         (1.0-self.u_probability))

    # Check values of weights - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.agree_weight != None) and (self.disagree_weight != None):
      if (self.agree_weight < self.disagree_weight):
        print 'error:Disagreement weight is larger than agreement weight'
        raise Exception

      if (self.agree_weight < self.missing_weight):
        print 'error:Agreement weight is smaller than missing weight'
        raise Exception

  # ---------------------------------------------------------------------------

  def set_probabilities(**kwargs):
    """Re-set m- and/or u-probabilities to new values. Agreement and
       disagrement weights are re-calculated unsing the new values.
    """

    for (keyword, value) in kwargs.items():
      if (keyword in ['m', 'm_prob', 'm_probability']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "m_probability" is not a number'
          raise Exception

        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability" value: %s' % (str(value))
          raise Exception
        self.m_probability = value

      elif (keyword in ['u', 'u_prob', 'u_probability']):
        if (not isinstance(value, int)) and (not isinstance(value, float)):
          print 'error:Argument "u_probability" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability" value: %s' % (str(value))
          raise Exception
        self.u_probability = value

      else:
        print 'error:Illegal constructor argument keyword: %s' % (str(keyword))
        raise Exception

    self.agree_weight =    mymath.log2(self.m_probability / self.u_probability)
    self.disagree_weight = mymath.log2((1.0-self.m_probability) /
                                       (1.0-self.u_probability))

    # Check values of weights - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.agree_weight != None) and (self.disagree_weight != None):
      if (self.agree_weight < self.disagree_weight):
        print 'error:Disagreement weight is larger than agreement weight'
        raise Exception

      if (self.agree_weight < self.missing_weight):
        print 'error:Agreement weight is smaller than missing weight'
        raise Exception

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b):
    """Compare two field lists to compute a weight value.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

  # ---------------------------------------------------------------------------

  def compute_basic_weights(self, fields_a, fields_b):
    """Check if field values are missing, are exactly the same, or differ.
       If a frequency table is given, compute weights according to the
       frequency of the given values, otherwise take generic weights based
       on M- and U- probability.

       This routine returns a dictionary with the keyword being either
         'missing'  and the missing weight as value
         'agree'    and the computed agreement weight as value
         'disagree' and two strings, agreement and disagreement weights as list
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        return {'missing':self.missing_weight}

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        return {'missing':self.missing_weight}

    # Concatenate into two strings without whitespaces  - - - - - - - - - - - -
    #
    str_a = ''.join(fields_a)
    str_b = ''.join(fields_b)

    # Check if strings are the same - - - - - - - - - - - - - - - - - - - - - -
    #
    if (str_a == str_b):

      if ((self.frequency_table != None) and \
          (self.frequency_table.has_key(str_a))):

        # String is found in frequency table, so compute value dependent weight
        #
        str_a_count = float(self.frequency_table.get(str_a))
        str_a_freq =  float(str_a_count / self.frequency_table.sum)

        # Agreement weight, computed according to L. Gill (2001), page 66
        #
        agree_weight = mymath.log2(1.0 / str_a_freq)

        if (self.freq_table_max_weight != None):  # Limit agreement weight
          if (agree_weight > self.freq_table_max_weight):
            agree_weight = self.freq_table_max_weight

        return {'agree':agree_weight}

      else:  # String not in frequency table, return generic agreement weight

        return {'agree':self.agree_weight}

    # Strings differ, so compute agreement and disagreement weights - - - - - -
    #
    if (self.frequency_table != None):
      if (self.frequency_table.has_key(str_a)):
        str_a_count = float(self.frequency_table.get(str_a))
        str_a_freq =  float(str_a_count / self.frequency_table.sum)

        agree_weight_a =    mymath.log2(1.0 / str_a_freq)
        disagree_weight_a = mymath.log2((1.0 - str_a_freq) / \
                                        (1.0 - str_a_freq*str_a_freq))
        if (self.freq_table_max_weight != None):  # Limit agreement weight
          if (agree_weight_a > self.freq_table_max_weight):
            agree_weight_a = self.freq_table_max_weight
        if (self.freq_table_min_weight != None):  # Limit disagreement weight
          if (disagree_weight_a < self.freq_table_min_weight):
            disagree_weight_a = self.freq_table_min_weight

      else:
        agree_weight_a =    self.agree_weight  # Take generic weights
        disagree_weight_a = self.disagree_weight

      if (self.frequency_table.has_key(str_b)):
        str_b_count = float(self.frequency_table.get(str_b))
        str_b_freq =  float(str_b_count / self.frequency_table.sum)

        agree_weight_b =    mymath.log2(1.0 / str_b_freq)
        disagree_weight_b = mymath.log2((1.0 - str_b_freq) / \
                                        (1.0 - str_b_freq*str_b_freq))
        if (self.freq_table_max_weight != None):  # Limit agreement weight
          if (agree_weight_b > self.freq_table_max_weight):
            agree_weight_b = self.freq_table_max_weight
        if (self.freq_table_min_weight != None):  # Limit disagreement weight
          if (disagree_weight_b < self.freq_table_min_weight):
            disagree_weight_b = self.freq_table_min_weight

      else:
        agree_weight_b =    self.agree_weight  # Take generic weights
        disagree_weight_b = self.disagree_weight

      # Now select minimal agreement and maximal disagreement weight - - - - -
      # ########### ??????? Is this correct ???????????? ############
      #
      agree_weight =    min(agree_weight_a, agree_weight_b)
      disagree_weight = max(disagree_weight_a, disagree_weight_b)

    # No frequency tables given so take generic weights - - - - - - - - - - - -
    #
    else:
      agree_weight =    self.agree_weight
      disagree_weight = self.disagree_weight

    return {'disagree':[str_a, str_b, agree_weight, disagree_weight]}

# =============================================================================

class FieldComparatorExactString(FieldComparator):
  """A field comparator based on exact string comparison.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No special attributes needed, just call base class constructor.
    """

    FieldComparator.__init__(self, kwargs)

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "ExactString" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    if (self.frequency_table != None):
      print '1:  Frequency table:          %s' % \
            (str(self.frequency_table.name))
    if (self.freq_table_min_weight != None):
      print '2:  Frequency minimal weight: %f' % (self.freq_table_min_weight)
    if (self.freq_table_max_weight != None):
      print '2:  Frequency maximal weight: %f' % (self.freq_table_max_weight)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two field lists using exact string comparator.
    """

    comparison_dict = self.compute_basic_weights(fields_a, fields_b)

    key =   comparison_dict.keys()[0]
    value = comparison_dict[key]

    if (key == 'agree'):  # Strings agree
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'missing'):  # Missing values
      print '3:%s    Missing values (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'disagree'):
      print '3:%s    Disagreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), value[3])
      return value[3]  # Return disagreement weight

    else:
      print 'error:%s: Illegal key returned from "compute_basic_weights"' % \
            (records_id)
      raise Exception

# =============================================================================

class FieldComparatorTruncateString(FieldComparator):
  """A field comparator based on exact string comparison with the limitation in
     number of characters that are compared (strings are truncated).

     The additional argument (besides the base class arguments) is
     'max_string_length' which gives the number of characters that are
     compared.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_string_length' argument first, then call
       the base class constructor.
    """

    self.max_string_length = None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in['max_string_len','max_string_length']):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "max_string_length" must be a positive '+ \
                'integer number'
          raise Exception
        self.max_string_length = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_string_length' attribute is set  - - - - - - - - - - - - -
    #
    if (self.max_string_length == None):
      print 'error:Argument "max_string_length" is not set'

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "TruncateString" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    if (self.frequency_table != None):
      print '1:  Frequency table:          %s' % \
            (str(self.frequency_table.name))
    if (self.freq_table_min_weight != None):
      print '2:  Frequency minimal weight: %f' % (self.freq_table_min_weight)
    if (self.freq_table_max_weight != None):
      print '2:  Frequency maximal weight: %f' % (self.freq_table_max_weight)
    print '1:  Maximum string length:    %i' % (self.max_string_length)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two field lists using exact string comparator with truncated
       strings.
    """

    comparison_dict = self.compute_basic_weights(fields_a, fields_b)

    key =   comparison_dict.keys()[0]
    value = comparison_dict[key]

    if (key == 'agree'):  # Strings agree
      print '3:%s    Agreement (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'missing'):  # Missing values
      print '3:%s    Missing values (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'disagree'):
      str_a =      value[0]
      str_b =      value[1]
      agree_w =    value[2]
      disagree_w = value[3]

      # Now check if truncated strings are the same - - - - - - - - - - - - - -
      #
      if (str_a[:self.max_string_length] == str_b[:self.max_string_length]):
        print '3:%s    Agreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), agree_w)
        return agree_w
      else:
        print '3:%s    Disagreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), disagree_w)
        return disagree_w

    else:
      print 'error:%s Illegal key returned from ' % (records_id) + \
            '"compute_basic_weights"'
      raise Exception

# =============================================================================

class FieldComparatorApproxString(FieldComparator):
  """Uses an approximate string comparator and frequency table if given.

     Approximate string comparator routines are implemented in the module
     'stringcmp.py'. Possible methods are:
       jaro
       winkler
       bigram
       editdist
       seqmatch
     They all return a value between 0.0 (strings differ totally) and 1.0
     (strings are the same).

     The additional arguments (besides the base class arguments) are the name
     of the comparison method (see above) 'compare_method' and the minimum
     value 'min_approx_value' (between 1.0 and 0.0) that is tolerated. If the
     approximate string comparator method returns a value between 1.0 and 
     'min_approx_value', then the partial agreement weight is computed using
     the following formula:

       weight = agree_weight - (1.0 - approx_str_val) /
                               (1.0 - min_approx_value) * 
                               (agree_weight+abs(disagree_weight))

     If the approximate string comparator method returns a value less than the
     'min_approx_value' then the disagreement weight is returned.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'compare_method' and 'min_approx_value'
       arguments first, then call the base class constructor.
    """

    self.compare_method =   None
    self.min_approx_value = None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'compare_method'):
        if (value not in ['jaro','winkler','bigram','editdist','seqmatch']):
          print 'error:Illegal approximate string comparator method: "%s"' % \
                (str(value))
          raise Exception
        self.compare_method = value

      elif (keyword == 'min_approx_value'):
        if (not isinstance(value,float)) or (value <= 0.0) or (value >= 1.0):
          print 'error:Argument "min_approc_value" must be a number '+ \
                'between 0.0 and 1.0'
          raise Exception
        self.min_approx_value = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'compare_method' and 'min_approx_value' attributes are set  - -
    #
    if (self.compare_method == None):
      print 'error:Argument "compare_method" is not set'
      raise Exception

    if (self.min_approx_value == None):
      print 'error:Argument "min_approx_value" is not set'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "ApproxString" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    if (self.frequency_table != None):
      print '1:  Frequency table:          %s' % \
            (str(self.frequency_table.name))
    if (self.freq_table_min_weight != None):
      print '2:  Frequency minimal weight: %f' % (self.freq_table_min_weight)
    if (self.freq_table_max_weight != None):
      print '2:  Frequency maximal weight: %f' % (self.freq_table_max_weight)
    print '1:  Compare method:           %s' % (str(self.compare_method))
    print '1:  Minimum comparison value: %f' % (self.min_approx_value)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two lists with strings using an approximate string comparator.
    """

    comparison_dict = self.compute_basic_weights(fields_a, fields_b)

    key =   comparison_dict.keys()[0]
    value = comparison_dict[key]

    if (key == 'agree'):  # Strings agree
      print '3:%s    Agreement (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'missing'):  # Missing values
      print '3:%s    Missing values (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'disagree'):
      str_a =      value[0]
      str_b =      value[1]
      agree_w =    value[2]
      disagree_w = value[3]

      if (self.compare_method == 'jaro'):
        approx_str_val = stringcmp.jaro(str_a, str_b)
      elif (self.compare_method == 'winkler'):
        approx_str_val = stringcmp.winkler(str_a, str_b)
      elif (self.compare_method == 'bigram'):
        approx_str_val = stringcmp.bigram(str_a, str_b)
      elif (self.compare_method == 'editdist'):
        approx_str_val = stringcmp.editdist(str_a, str_b)
      elif (self.compare_method == 'seqmatch'):
        approx_str_val = stringcmp.seqmatch(str_a, str_b)
      else:
        print 'error:%s Illegal approximate string comparator ' % \
              (records_id) + 'method: "%s"' % (self.compare_method)
        raise Exception

      # Check if strings are totally different  - - - - - - - - - - - - - - - -
      #
      if (approx_str_val < self.min_approx_value):
        print '3:%s    Disagreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), disagree_w)
        return disagree_w

      # Compute final adjusted weight - - - - - - - - - - - - - - - - - - - - -
      # Modified after formula in Winkler and Thibaudeau, 1991, page 12
      #
      partagree_w = agree_w - (1.0 - approx_str_val) / \
                    (1.0 - self.min_approx_value) * \
                    (agree_w+abs(disagree_w))
      print '3:%s    Partial agreement (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), partagree_w)

      return partagree_w

    else:
      print 'error:%s Illegal key returned from ' % (records_id) + \
            '"compute_basic_weights"'
      raise Exception

# =============================================================================

class FieldComparatorEncodeString(FieldComparator):
  """Uses a string encoding comparator and frequency table if given.

     String encoding routines are implemented in the module 'encode.py'.
     Possible methods are:
       soundex      Soundex
       mod_soundex  Modified Soundex
       phonex       Phonex
       nysiis       NYSIIS
       dmetaphone   Double-Metaphone

     If both strings are encoded the same way, the agreement weight is
     returned, otherwise the disagreement weight.

     The additional arguments (besides the base class arguments) are the name
     of the encoding method (see above) 'encode_method', a flag 'reverse' which
     if set to 'True' the strings are reversed first before they are encoded.
     The default value for 'reverse' is 'False' (do not reverse strings). The
     'max_code_length' argument can be used to set the maximal length (in
     characters) of the codes. Default is 4.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'encode_method' and 'reverse' arguments first,
       then call the base class constructor.
    """

    self.encode_method =   None
    self.reverse =         False
    self.max_code_length = 4

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'encode_method'):
        if (value not in ['soundex','mod_soundex','phonex','nysiis', \
                          'dmetaphone']):
          print 'error:Illegal string encoding method: "%s"' % (str(value))
          raise Exception
        self.encode_method = value

      elif (keyword == 'reverse'):
        if (value not in [True, False]):
          print 'error:Argument "reverse" must be "True" or "False"'
          raise Exception
        self.reverse = value

      elif (keyword == 'max_code_length'):
        if (not isinstance(value, int)) or (value <= 0):
          print 'error:Argument "max_code_length" must be a positive '+ \
                'integer number'
          raise Exception
        self.max_code_length = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'encode_method' attribute is set  - - - - - - - - - - - - - - -
    #
    if (self.encode_method == None):
      print 'error:Argument "encode_method" is not set'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "EncodeString" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    if (self.frequency_table != None):
      print '1:  Frequency table:          %s' % \
            (str(self.frequency_table.name))
    if (self.freq_table_min_weight != None):
      print '2:  Frequency minimal weight: %f' % (self.freq_table_min_weight)
    if (self.freq_table_max_weight != None):
      print '2:  Frequency maximal weight: %f' % (self.freq_table_max_weight)
    print '1:  Encode method:            %s' % (str(self.encode_method))
    print '1:  Maximum code length:      %i' % (self.max_code_length)
    print '1:  String reversed:          %s' % (str(self.reverse))

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two list with strings using a string encoding comparator.
    """

    comparison_dict = self.compute_basic_weights(fields_a, fields_b)

    key =   comparison_dict.keys()[0]
    value = comparison_dict[key]

    if (key == 'agree'):  # Strings agree
      print '3:%s    Agreement (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'missing'):  # Missing values
      print '3:%s    Missing values (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'disagree'):
      str_a =      value[0]
      str_b =      value[1]
      agree_w =    value[2]
      disagree_w = value[3]

      if (self.reverse == True):  # Reverse both strings
        rev = list(str_a)
        rev.reverse()
        str_a = ''.join(rev)
        rev = list(str_b)
        rev.reverse()
        str_b = ''.join(rev)

      if (self.encode_method == 'soundex'):
        code_a = encode.soundex(str_a, maxlen=self.max_code_length)
        code_b = encode.soundex(str_b, maxlen=self.max_code_length)
      elif (self.encode_method == 'mod_soundex'):
        code_a = encode.mod_soundex(str_a, maxlen=self.max_code_length)
        code_b = encode.mod_soundex(str_b, maxlen=self.max_code_length)
      elif (self.encode_method == 'phonex'):
        code_a = encode.phonex(str_a, maxlen=self.max_code_length)
        code_b = encode.phonex(str_b, maxlen=self.max_code_length)
      elif (self.encode_method == 'nysiis'):
        code_a = encode.nysiis(str_a, maxlen=self.max_code_length)
        code_b = encode.nysiis(str_b, maxlen=self.max_code_length)
      elif (self.encode_method == 'dmetaphone'):
        code_a = encode.dmetaphone(str_a, maxlen=self.max_code_length)
        code_b = encode.dmetaphone(str_b, maxlen=self.max_code_length)
      else:
        print 'error:%s Illegal string encoding ' % (records_id) + \
              'method: "%s"' % (self.encode_method)
        raise Exception

      # Check if encodings are the same or different  - - - - - - - - - - - - -
      #
      if (code_a == code_b):
        print '3:%s    Agreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), agree_w)
        return agree_w

      else:
        print '3:%s    Disagreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), disagree_w)
        return disagree_w

    else:
      print 'error:%s Illegal key returned from ' % (records_id) + \
            '"compute_basic_weights"'
      raise Exception

# =============================================================================

class FieldComparatorKeyDiff(FieldComparator):
  """A field comparator that compares the fields character-wise with a maximum
     number of errors (different characters) being tolerated. This comparator
     can be used to compare numerical fields, such as date of birth, telephone
     numbers, etc.

     The additional argument (besides the base class arguments) is
     'max_key_diff' which must be an integer that gives the maximal number of
     maximal errors tolerated (i.e. different characters).

     For example, if 'max_key_diff' is set to 'X' and 'Y' (with 0 < Y <= X)
     different characters are encountered in the two strings, then the
     resulting weight will be computed as:

       weight = agree_weight - (Y/(X+1))*(agree_weight+abs(disagree_weight))

     If the number of errors encountered is larger than X, the disagreement
     weight will be returned.     
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_key_diff' argument first, then call the
       base class constructor.
    """

    self.max_key_diff = None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_key_diff'):
        if (not isinstance(value, int)) or (value < 0):
          print 'error:Argument "max_key_diff" must be an integer number '+ \
                'equal to or larger than 0'
          raise Exception
        self.max_key_diff = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_key_diff' attribute is set - - - - - - - - - - - - - - - -
    #
    if (self.max_key_diff == None):
      print 'Argument "max_key_diff" is not set'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "KeyDiff" field comparator: "%s"' % (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    if (self.frequency_table != None):
      print '1:  Frequency table:          %s' % \
            (str(self.frequency_table.name))
    if (self.freq_table_min_weight != None):
      print '2:  Frequency minimal weight: %f' % (self.freq_table_min_weight)
    if (self.freq_table_max_weight != None):
      print '2:  Frequency maximal weight: %f' % (self.freq_table_max_weight)
    print '1:  Maximum key difference:   %i' % (self.max_key_diff)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two field lists using the key error field comparator.
    """

    comparison_dict = self.compute_basic_weights(fields_a, fields_b)

    key =   comparison_dict.keys()[0]
    value = comparison_dict[key]

    if (key == 'agree'):  # Strings agree
      print '3:%s    Agreement (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'missing'):  # Missing values
      print '3:%s    Missing values (%s, %s): %f' % \
            (records_id, str(fields_a), str(fields_b), value)
      return value

    elif (key == 'disagree'):
      str_a =      value[0]
      str_b =      value[1]
      agree_w =    value[2]
      disagree_w = value[3]

      # The initial number of errors is the difference in the strings' length
      #
      num_err = abs(len(str_a) - len(str_b))

      check_len = min(len(str_a), len(str_b))

      for i in range(check_len):  # Loop over positions in strings
        if (str_a[i] != str_b[i]):
          num_err += 1

      if (num_err > self.max_key_diff):
        print '3:%s    Disagreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), disagree_w)
        return disagree_w  # Too many errors, return disagree weight

      else:
        partagree_w = agree_w - (float(num_err)/(self.max_key_diff+1.0)) * \
                      (agree_w + abs(disagree_w))
        print '3:%s    Partial agreement (%s, %s): %f' % \
              (records_id, str(fields_a), str(fields_b), partagree_w)
        return partagree_w

    else:
      print 'error:%s Illegal key returned from ' % (records_id) + \
            '"compute_basic_weights"'
      raise Exception

# =============================================================================

class FieldComparatorNumericPerc(FieldComparator):
  """A field comparator for numeric fields, where a given percentage difference
     can be tolerated. The agreement weight is returned if the numbers are the
     same, and the disagreement weight if the percentage difference is larger
     than the 'max_perc_diff' value.

     The additional argument (besides the base class arguments) is
     'max_perc_diff' which must be a floating-point number between 0.0 and
     100.0. Default is 0.0 - in which case only an exact comparison is
     performed.

     If the percentage difference is smaller than 'max_perc_diff' the weight is
     calculated according to the following formula:

       weight = agree_weight - (perc_diff / (max_perc_diff + 1.0)) * 
                                (agree_weight+abs(disagree_weight))

     where the percentage difference is calculated as:

       perc_diff = abs(value_a - value_b) / min(value_a, value_b) * 100.0

     Currently, no frequency table can be used for this comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_perc_diff' argument first, then call the
       base class constructor.
    """

    self.max_perc_diff = 0.0  # Default only exact numeric comparison

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorNumericPerc" field comparator'
      raise Exception

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_perc_diff'):
        if ((not (isinstance(value, int) or isinstance(value, float))) or \
            (value < 0.0) or (value >= 100.0)):
          print 'error:Argument "max_perc_diff" must be a number between '+ \
                '0.0 and 100.0'
          raise Exception
        self.max_perc_diff = float(value)

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "NumericPerc" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    print '1:  Maximum percentage difference: %f' % (self.max_perc_diff)

  # ---------------------------------------------------------------------------

  def compare(self, fields_a, fields_b, records_id):
    """Compare two fields with numeric values.
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Concatenate fields into two strings without whitespaces - - - - - - - - -
    #
    str_a = ''.join(fields_a)
    str_b = ''.join(fields_b)

    # Check if field values are numbers - - - - - - - - - - - - - - - - - - - -
    #
    try:
      value_a = float(str_a)
    except:
      value_a = ''  # Set to empty string
    try:
      value_b = float(str_b)
    except:
      value_b = ''  # Set to empty string

    # If one or both are not numbers, return diagreement weight - - - - - - - -
    #
    if (value_a == '') or (value_b == ''):
      print '3:%s    Disagreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.disagree_weight)
      return self.disagree_weight

    # Check if numbers are the same - - - - - - - - - - - - - - - - - - - - - -
    #
    elif (value_a == value_b):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    # Calculate percentage difference and weight  - - - - - - - - - - - - - - -
    #
    else:
      if (self.max_perc_diff == 0.0):  # No percentage tolerance allowed
        print '3:%s    Disagreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.disagree_weight)
        return self.disagree_weight  # Because values are different

      else:
        perc_diff = abs(value_a - value_b) / min(value_a, value_b) * 100.0

        if (perc_diff > self.max_perc_diff):  # Percentage difference too large
          print '3:%s    Disagreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), \
                 self.disagree_weight)
          return self.disagree_weight

        else:
          partagree_w = self.agree_weight - (perc_diff / \
                        (self.max_perc_diff+1.0)) * \
                        (self.agree_weight+ \
                        abs(self.disagree_weight))
          print '3:%s    Partial agreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), partagree_w)
          return partagree_w

# =============================================================================

class FieldComparatorNumericAbs(FieldComparator):
  """A field comparator for numeric fields, where a given absolute difference
     can be tolerated. The agreement weight is returned if the numbers are the
     same, and the disagreement weight if the absolute difference is larger
     than the 'max_abs_diff' value. Default value for 'max_abs_diff' is 0.

     The additional argument (besides the base class arguments) is
     'max_abs_diff' which can be a positive integer or floating-point number.

     If the percentage difference is equal or smaller than 'max_abs_diff' the
     weight is calculated according to the following formula:

       weight = agree_weight - (abs_diff / (max_abs_diff + 1.0)) *
                                (agree_weight+abs(disagree_weight))

     where the absolute difference is calculated as:

       abs_diff = abs(value_a - value_b)

     Currently, no frequency table can be used for this comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_abs_diff' argument first, then call the
       base class constructor.
    """

    self.max_abs_diff = 0.0  # Default only exact numeric comparison

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorNumericAbs" field comparator'
      raise Exception

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_abs_diff'):
        if ((not (isinstance(value, int) or isinstance(value, float))) or \
            (value < 0.0)):
          print 'error:Argument "max_abs_diff" must be a positive number'
          raise Exception
        self.max_abs_diff = float(value)

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "NumericAbs" field comparator: "%s"' % \
          (str(self.name))
    print '1:  Fields in data set A:     %s' % (str(self.fields_a))
    print '1:  Fields in data set B:     %s' % (str(self.fields_b))
    print '1:  M-Probability:            %f' % (self.m_probability)
    print '1:  U-Probability:            %f' % (self.u_probability)
    print '2:  Agreement weight:         %f' % (self.agree_weight)
    print '2:  Disagreement weight:      %f' % (self.disagree_weight)
    print '1:  Missing weight:           %f' % (self.missing_weight)
    print '1:  Maximum absolute difference: %f' % (self.max_abs_diff)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two fields with numeric values.
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Concatenate fields into two strings without whitespaces - - - - - - - - -
    #
    str_a = ''.join(fields_a)
    str_b = ''.join(fields_b)

    # Check if field values are numbers - - - - - - - - - - - - - - - - - - - -
    #
    try:
      value_a = float(str_a)
    except:
      value_a = ''  # Set to empty string
    try:
      value_b = float(str_b)
    except:
      value_b = ''  # Set to empty string

    # If one or both are not numbers, return disagreement weight  - - - - - - -
    #
    if (value_a == '') or (value_b == ''):
      print '3:%s    Disagreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.disagree_weight)
      return self.disagree_weight

    # Check if numbers are the same - - - - - - - - - - - - - - - - - - - - - -
    #
    elif (value_a == value_b):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    # Calculate absolute difference and weight  - - - - - - - - - - - - - - - -
    #
    else:
      if (self.max_abs_diff == 0.0):  # No absolute difference allowed
        print '3:%s    Disagreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.disagree_weight)
        return self.disagree_weight  # Because values are different

      else:
        abs_diff = abs(value_a - value_b)

        if (abs_diff > self.max_abs_diff):  # Absolute difference too large
          print '3:%s    Disagreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), \
                 self.disagree_weight)
          return self.disagree_weight

        else:
          partagree_w = self.agree_weight - (abs_diff / \
                        (self.max_abs_diff+1.0)) * \
                        (self.agree_weight+ \
                        abs(self.disagree_weight))
          print '3:%s    Partial agreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), partagree_w)
          return partagree_w

# =============================================================================

class FieldComparatorDate(FieldComparator):
  """A field comparator that compares two dates, which must be given in three
     separate fields (day,month,year). Two additional arguments can be given:

     'max_day_a_before_b' is the maximal tolerated number that day A can be
     before day B
     'max_day_b_before_a' is the maximal tolerated number that day B can be
     before day A

     Default values for both are 0, which results in exact day comparison only
     (i.e. if the dates are not the same the disagreement value will be
     returned).

     A partial agreement will be calculated as follows:

       weight = agree_weight - (date_diff/(max_day_a_before_b+1)) *
                               (agree_weight+abs(disagree_weight))

     and similar for 'max_day_a_before_b'.

     If the day difference is larger than 'max_day_a_before_b' or
     'max_day_b_before_a' the disagreement weight will be returned.

     Different m- and u-probabilities have to be given for day, month and year.
     The general 'm_probability' and 'u_probability' argument can not be used.
     Instead, use:
       m_probability_day
       u_probability_day
       m_probability_month
       u_probability_month
       m_probability_year
       u_probability_year

     Currently, no frequency table can be used for this comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the date arguments first, then call the base
       class constructor.
    """

    self.max_day_a_before_b =    0  # Default only exact date comparison
    self.max_day_b_before_a =    0  # Default only exact date comparison
    self.m_probability_day =     None
    self.u_probability_day =     None
    self.m_probability_month =   None
    self.u_probability_month =   None
    self.m_probability_year =    None
    self.u_probability_year =    None
    self.agree_weight_day =      None
    self.disagree_weight_day =   None
    self.agree_weight_month =    None
    self.disagree_weight_month = None
    self.agree_weight_year =     None
    self.disagree_weight_year =  None

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorDate" field comparator'
      raise Exception

    # Check for arguments 'm_probability' and 'u_probability' - - - - - - - - -
    #
    if (('m_probability' in kwargs.keys()) or \
        ('u_probability' in kwargs.keys())):
      print 'error:Arguments "m_probability" and "u_probablity" can not be '+ \
            'use with "FieldComparatorDate" field comparator. Use day, '+ \
            'month and year specific m- and u-probabilities instead'
      raise Exception

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_day_a_before_b'):
        if (not isinstance(value, int) or (value < 0)):
          print 'error:Argument "max_day_a_before_b" must be a positive '+ \
                'integer number'
          raise Exception
        self.max_day_a_before_b = value
      elif (keyword == 'max_day_b_before_a'):
        if (not isinstance(value, int) or (value < 0)):
          print 'error:Argument "max_day_b_before_a" must be a positive '+ \
                'integer number'
          raise Exception
        self.max_day_b_before_a = value

      elif (keyword in ['m_day', 'm_prob_day', 'm_probability_day']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_day" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_day" value: %f' % (value)
          raise Exception
        self.m_probability_day = value

      elif (keyword in ['u_day', 'u_prob_day', 'u_probability_day']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_day" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_day" value: %f' % (value)
          raise Exception
        self.u_probability_day = value

      elif (keyword in ['m_month', 'm_prob_month', 'm_probability_month']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_month" is not a number'
          raise Exception
        if (value < 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_month" value: %f' % (value)
          raise Exception
        self.m_probability_month = value

      elif (keyword in ['u_month', 'u_prob_month', 'u_probability_month']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_month" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_month" value: %f' % (value)
          raise Exception
        self.u_probability_month = value

      elif (keyword in ['m_year', 'm_prob_year', 'm_probability_year']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_year" is not a number'
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_year" value: %f' % (value)
          raise Exception
        self.m_probability_year = value

      elif (keyword in ['u_year', 'u_prob_year', 'u_probability_year']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_year" is not a number'
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_year" value: %f' % (value)
          raise Exception
        self.u_probability_year = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Compute weights if probabilities are given - - - - - - - - - - - - - - -
    #
    if (self.m_probability_day != None) and (self.u_probability_day != None):
      self.agree_weight_day =    mymath.log2(self.m_probability_day / \
                                             self.u_probability_day)
      self.disagree_weight_day = mymath.log2((1.0-self.m_probability_day) /
                                             (1.0-self.u_probability_day))

    if (self.m_probability_month != None) and (self.u_probability_month!=None):
      self.agree_weight_month =    mymath.log2(self.m_probability_month / \
                                               self.u_probability_month)
      self.disagree_weight_month = mymath.log2((1.0-self.m_probability_month) /
                                               (1.0-self.u_probability_month))

    if (self.m_probability_year != None) and (self.u_probability_year != None):
      self.agree_weight_year =    mymath.log2(self.m_probability_year / \
                                              self.u_probability_year)
      self.disagree_weight_year = mymath.log2((1.0-self.m_probability_year) /
                                              (1.0-self.u_probability_year))

    self.agree_weight = self.agree_weight_day + self.agree_weight_month + \
                        self.agree_weight_year
    self.disagree_weight = self.disagree_weight_day + \
                           self.disagree_weight_month + \
                           self.disagree_weight_year

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "Date" field comparator: "%s"' % (str(self.name))
    print '1:  Fields in data set A:       %s' % (str(self.fields_a))
    print '1:  Fields in data set B:       %s' % (str(self.fields_b))
    print '1:  M-Probability day:          %f' % (self.m_probability_day)
    print '1:  U-Probability day:          %f' % (self.u_probability_day)
    print '1:  M-Probability month:        %f' % (self.m_probability_month)
    print '1:  U-Probability month:        %f' % (self.u_probability_month)
    print '1:  M-Probability year:         %f' % (self.m_probability_year)
    print '1:  U-Probability year:         %f' % (self.u_probability_year)
    print '2:  Agreement weight day:       %f' % (self.agree_weight_day)
    print '2:  Disagreement weight day:    %f' % (self.disagree_weight_day)
    print '2:  Agreement weight month:     %f' % (self.agree_weight_month)
    print '2:  Disagreement weight month:  %f' % (self.disagree_weight_month)
    print '2:  Agreement weight year:      %f' % (self.agree_weight_year)
    print '2:  Disagreement weight yead:   %f' % (self.disagree_weight_year)
    print '1:  Missing weight:             %f' % (self.missing_weight)
    print '1:  Maximal day A before day B: %i' % (self.max_day_a_before_b)
    print '1:  Maximal day B before day A: %i' % (self.max_day_b_before_a)

  # ---------------------------------------------------------------------------

  def compare(self, fields_a, fields_b, records_id):
    """Compare two date field lists. The format of fields_a and fields_b must
       be [day,month,year].
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Check if fields contain three elements, then check day and month values -
    #
    if (not (isinstance(fields_a,list) or isinstance(fields_a,tuple))):
      print 'error:%s Fields for data A are not a list or a tuple' % \
            (records_id)
      raise Exception
    if (not (isinstance(fields_b,list) or isinstance(fields_b,tuple))):
      print 'error:%s Fields for data B are not a list or a tuple' % \
            (records_id)
      raise Exception

    if (len(fields_a) != 3) or (len(fields_b) != 3):
      print 'error:%s Fields are not in the right format ' % (records_id) + \
            '[day,month,year]'
      raise Exception

    # Integer conversion in comparison by Dan Sabath
    #
    if ((int(fields_a[0]) < 1) or (int(fields_a[0]) > 31) or \
        (int(fields_a[1]) < 1) or (int(fields_a[1]) > 12)):
      print 'error:%s Fields A does not have valid day or month values' % \
            (records_id)
      raise Exception

    if ((int(fields_b[0]) < 1) or (int(fields_b[0]) > 31) or \
        (int(fields_b[1]) < 1) or (int(fields_b[1]) > 12)):
      print 'error:%s Fields B does not have valid day or month values' % \
            (records_id)
      raise Exception

    # Compute difference between dates in days  - - - - - - - - - - - - - - - -
    #
    date_a_epoch = date.date_to_epoch(fields_a[0], fields_a[1], fields_a[2])
    date_b_epoch = date.date_to_epoch(fields_b[0], fields_b[1], fields_b[2])

    day_diff = date_a_epoch - date_b_epoch

    # Check if dates are the same - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (day_diff == 0):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    elif (day_diff > 0):  # Day A after day B

      if (day_diff > self.max_day_b_before_a):  # To many days before
        print '3:%s    Disagreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.disagree_weight)
        return self.disagree_weight

      else:
        partagree_w = self.agree_weight - \
                      (float(day_diff) / (self.max_day_b_before_a+1.0)) * \
                      (self.agree_weight + abs(self.disagree_weight))
        print '3:%s    Partial agreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), partagree_w)
        return partagree_w

    else:  # Day B after day A (day_diff < 0)

      if (-day_diff > self.max_day_a_before_b):  # To many days before
        print '3:%s    Disagreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.disagree_weight)
        return self.disagree_weight

      else:
        partagree_w = self.agree_weight - \
                      (float(-day_diff) / (self.max_day_a_before_b+1.0)) * \
                      (self.agree_weight + abs(self.disagree_weight))
        print '3:%s    Partial agreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), partagree_w)
        return partagree_w

# =============================================================================

class FieldComparatorAge(FieldComparator):
  """A field comparator that compares ages given either as numbers (assumed to
     be ages) or dates, which are converted into ages. The comparator allows
     for a certain percentage difference between them. A 'fixe_date' can be
     given relative to which an age is computed if it is given as a date.
     Default value is 'None' in which case ages are calculated relative to the
     current day. The 'fix_date' can either be given as a date string, as a
     date tuple (in the form [day,month.year]) or the string 'today'.

     An additional argument can be 'max_perc_diff' which gives the maximum
     tolarated age difference in percentage. The value must be a floating-point
     number between 0.0 and 100.0. Default is 0.0 - in which case an exact
     comparison is performed.

     If the percentage difference is smaller than 'max_perc_diff' the weight is
     calculated according to the following formula:

       weight = agree_weight - (perc_diff / (max_perc_diff + 1.0)) * 
                                (agree_weight+abs(disagree_weight))

     where the percentage difference is calculated as:

       perc_diff = abs(age_a - age_b) / min(age_a, age_b) * 100.0

     Different m- and u-probabilities have to be given for day, month and year.
     The general 'm_probability' and 'u_probability' argument can not be used.
     Instead, use:
       m_probability_day
       u_probability_day
       m_probability_month
       u_probability_month
       m_probability_year
       u_probability_year

     Currently, no frequency table can be used for this comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the date arguments first, then call the base
       class constructor.
    """

    self.fix_date =              None
    self.max_perc_diff =         0.0  # Default only exact comparison
    self.m_probability_day =     None
    self.u_probability_day =     None
    self.m_probability_month =   None
    self.u_probability_month =   None
    self.m_probability_year =    None
    self.u_probability_year =    None
    self.agree_weight_day =      None
    self.disagree_weight_day =   None
    self.agree_weight_month =    None
    self.disagree_weight_month = None
    self.agree_weight_year =     None
    self.disagree_weight_year =  None

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorAge" field comparator'
      raise Exception

    # Check for arguments 'm_probability' and 'u_probability' - - - - - - - - -
    #
    if (('m_probability' in kwargs.keys()) or \
        ('u_probability' in kwargs.keys())):
      print 'error:Arguments "m_probability" and "u_probablity" can not be '+ \
            'use with "FieldComparatorAge" field comparator. Use day, '+ \
            'month and year specific m- and u-probabilities instead'
      raise Exception

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'fix_date'):
        self.fix_date = value

      elif (keyword == 'max_perc_diff'):
        if (not ((isinstance(value, int) or isinstance(value, float))) or \
            (value <= 0.0) or (value > 100.0)):
          print 'error:Argument "max_perc_diff" must be a number between '+ \
                '0.0 and 100.0'
          raise Exception
        self.max_perc_diff = float(value)

      elif (keyword in ['m_day', 'm_prob_day', 'm_probability_day']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_day" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_day" value: %f' % (value)
          raise Exception
        self.m_probability_day = value

      elif (keyword in ['u_day', 'u_prob_day', 'u_probability_day']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_day" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_day" value: %f' % (value)
          raise Exception
        self.u_probability_day = value

      elif (keyword in ['m_month', 'm_prob_month', 'm_probability_month']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_month" is not a number'
          raise Exception
        if (value < 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_month" value: %f' % (value)
          raise Exception
        self.m_probability_month = value

      elif (keyword in ['u_month', 'u_prob_month', 'u_probability_month']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_month" is not a number'
          raise Exception
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_month" value: %f' % (value)
          raise Exception
        self.u_probability_month = value

      elif (keyword in ['m_year', 'm_prob_year', 'm_probability_year']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "m_probability_year" is not a number'
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "m_probability_year" value: %f' % (value)
          raise Exception
        self.m_probability_year = value

      elif (keyword in ['u_year', 'u_prob_year', 'u_probability_year']):
        if (not (isinstance(value, int) or isinstance(value, float))):
          print 'error:Argument "u_probability_year" is not a number'
        if (value <= 0.0) or (value > 1.0):
          print 'error:Illegal "u_probability_year" value: %f' % (value)
          raise Exception
        self.u_probability_year = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Process fix date  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.fix_date == None) or (self.fix_date == 'today'):
      today_tuple = time.localtime(time.time())
      self.fix_date = [today_tuple[2], today_tuple[1], today_tuple[0]]

    elif (isinstance(self.fix_date, str)):
      tmp_date = date.str_to_date(self.fix_date)  # Parse string into tuple
      if (tmp_date == []):
        print 'error:Illegal "fix_date" value, no valid date: "%s"' % \
              (str(self.fix_date))
        raise Exception
      self.fix_date = tmp_date

    elif (isinstance(self.fix_date,list)) or (isinstance(self.fix_date,tuple)):
      if (len(self.fix_date) != 3):
        print 'error:Illegal list format of "fix_date", must have three ' + \
              'elements (day,month,year): "%s"' % (str(self.fix_date))
        raise Exception
    else:
      print 'error:Illegal "fix_date" format: "%s"' % (str(self.fix_date))
      raise Exception

    # Now convert the fix date tuple into corresponding epoch number
    #
    self.fix_date = date.date_to_epoch(self.fix_date[0], self.fix_date[1],
                                       self.fix_date[2])

    # Compute weights if probabilities are given  - - - - - - - - - - - - - - -
    #
    if (self.m_probability_day != None) and (self.u_probability_day != None):
      self.agree_weight_day =    mymath.log2(self.m_probability_day / \
                                             self.u_probability_day)
      self.disagree_weight_day = mymath.log2((1.0-self.m_probability_day) /
                                             (1.0-self.u_probability_day))

    if (self.m_probability_month != None) and (self.u_probability_month!=None):
      self.agree_weight_month =    mymath.log2(self.m_probability_month / \
                                               self.u_probability_month)
      self.disagree_weight_month = mymath.log2((1.0-self.m_probability_month) /
                                               (1.0-self.u_probability_month))

    if (self.m_probability_year != None) and (self.u_probability_year != None):
      self.agree_weight_year =    mymath.log2(self.m_probability_year / \
                                              self.u_probability_year)
      self.disagree_weight_year = mymath.log2((1.0-self.m_probability_year) /
                                              (1.0-self.u_probability_year))

    self.agree_weight = self.agree_weight_day + self.agree_weight_month + \
                        self.agree_weight_year
    self.disagree_weight = self.disagree_weight_day + \
                           self.disagree_weight_month + \
                           self.disagree_weight_year

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "Age" field comparator: "%s"' % (str(self.name))
    print '1:  Fields in data set A:          %s' % (str(self.fields_a))
    print '1:  Fields in data set B:          %s' % (str(self.fields_b))
    print '1:  M-Probability day:             %f' % (self.m_probability_day)
    print '1:  U-Probability day:             %f' % (self.u_probability_day)
    print '1:  M-Probability month:           %f' % (self.m_probability_month)
    print '1:  U-Probability month:           %f' % (self.u_probability_month)
    print '1:  M-Probability year:            %f' % (self.m_probability_year)
    print '1:  U-Probability year:            %f' % (self.u_probability_year)
    print '2:  Agreement weight day:          %f' % (self.agree_weight_day)
    print '2:  Disagreement weight day:       %f' % (self.disagree_weight_day)
    print '2:  Agreement weight month:        %f' % (self.agree_weight_month)
    print '2:  Disagreement weight month:     %f' % \
          (self.disagree_weight_month)
    print '2:  Agreement weight year:         %f' % (self.agree_weight_year)
    print '2:  Disagreement weight yead:      %f' % (self.disagree_weight_year)
    print '1:  Missing weight:                %f' % (self.missing_weight)
    print '1:  Fix date:                      %s' % (str(self.fix_date))
    print '1:  Maximal percentage difference: %s' % (str(self.max_perc_diff))

  # ---------------------------------------------------------------------------

  def compare(self, fields_a, fields_b, records_id):
    """Compare two age fields or date field lists. The format of fields_a and
       fields_b must either be a number (which is taken to be an age in years)
       or a [day,month,year] list.
    """

    if ('00' in fields_a) or ('00' in fields_b):
      print 'warning:%s fields_a: %s, fields_b: %s' \
            % (records_id, str(fields_a),str(fields_a))
 
    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Check if fields are either a date or an age - - - - - - - - - - - - - - -
    #
    if (len(fields_a) == 1):
      try:
        age_a = float(fields_a[0])
      except:
        print 'warning:%s Illegal age in fields A: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_a[0]))
        return self.missing_weight

    elif (len(fields_a) == 3):
      age_a = None
      try:
        day_a = int(fields_a[0])
      except:
        print 'warning:%s Illegal day in fields A: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_a[0]))
        return self.missing_weight

      try:
        month_a = int(fields_a[1])
      except:
        print 'warning:%s Illegal month in fields A: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_a[1]))
        return self.missing_weight

      try:
        year_a = int(fields_a[2])
      except:
        print 'warning:%s Illegal year in fields A: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_a[2]))
        return self.missing_weight

    else:
      print 'warning:%s Fields A is not an age or a date with ' % \
            (records_id) + 'format [day,month,year]: %s' % (str(fields_a)) + \
            ', set weight to missing value'
      return self.missing_weight

    if (len(fields_b) == 1):
      try:
        age_b = float(fields_b[0])
      except:
        print 'warning:%s Illegal day in fields B: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_b[0]))
        return self.missing_weight

    elif (len(fields_b) == 3):
      age_b = None
      try:
        day_b = int(fields_b[0])
      except:
        print 'warning:%s Illegal day in fields B: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_b[0]))
        return self.missing_weight

      try:
        month_b = int(fields_b[1])
      except:
        print 'warning:%s Illegal month in fields B: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_b[1]))
        return self.missing_weight

      try:
        year_b = int(fields_b[2])
      except:
        print 'warning:%s Illegal year in fields B: ' % (records_id) + \
              '%s, set weight to missing value' % (str(fields_b[2]))
        return self.missing_weight

    else:
      print 'warning:%s Fields B is not an age or a date with ' % \
            (records_id) + 'format [day,month,year]: %s' % (str(fields_b)) + \
            ', set weight to missing value'
      return self.missing_weight

    # Check some values for dates and calculate age - - - - - - - - - - - - - -
    #
    if (age_a == None):
      if ((day_a < 1) or (day_a > 31) or (month_a < 1) or (month_a > 12)):
        print 'warning:%s Fields A do not have valid day or' % (records_id) + \
              ' month values: %s, set weight to missing value' % \
              (str(fields_a))
        return self.missing_weight

      age_a = date.date_to_age(day_a, month_a, year_a, self.fix_date)

    if (age_b == None):
      if ((day_b < 1) or (day_b > 31) or (month_b < 1) or (month_b > 12)):
        print 'warning:%s Fields A do not have valid day or' % (records_id) + \
              ' month values: %s, set weight to missing value' % \
              (str(fields_b))
        return self.missing_weight

      age_b = date.date_to_age(day_b, month_b, year_b, self.fix_date)

    # Check age values  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (age_a < 0.0) or (age_a > 150.0):
      print 'warning:%s Illegal value for age A: %s' % \
            (records_id, str(age_a)) + ', set weight to missing value'
      return self.missing_weight

    if (age_b < 0.0) or (age_b > 150.0):
      print 'warning:%s Illegal value for age B: %s' % \
            (records_id, str(age_b)) + ', set weight to missing value'
      return self.missing_weight

    # Check if ages are the same  - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (age_a == age_b):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    # Calculate percentage difference and weight  - - - - - - - - - - - - - - -
    #
    else:
      if (self.max_perc_diff == 0.0):  # No percentage tolerance allowed
        print '3:%s    Disagreement (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.disagree_weight)
        return self.disagree_weight  # Because values are different

      else:
        perc_diff = abs(age_a - age_b) / min(age_a, age_b) * 100.0

        if (perc_diff > self.max_perc_diff):  # Percentage difference too large
          print '3:%s    Disagreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), \
                 self.disagree_weight)
          return self.disagree_weight

        else:
          partagree_w =  self.agree_weight - (perc_diff / \
                         (self.max_perc_diff+1.0)) * \
                         (self.agree_weight+ abs(self.disagree_weight))
          print '3:%s    Partial agreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), partagree_w)
          return partagree_w

# =============================================================================

class FieldComparatorTime(FieldComparator):
  """A field comparator for time fields, which must be given in 24 hours format
     (00:00 is midnight and 23:59 is 11:59 pm).

     The concatenated fields must be of the form 'HHMM' or 'HH:MM', otherwise
     an error is triggered.

     Three additional arguments can be given:

     'max_time_a_before_b' is the maximal tolerated time in minutes that time A
     can be before time B
     'max_time_b_before_a' is the maximal tolerated time in minutes that time B
     can be before time A
     'day_start' is the time (HHMM) where a 24-hours period starts. It is
     assumed that both times are within the same 24-hours period. The default
     value is midnight (00:00).

     Default values for both 'max_time_a_before_b' and 'max_time_b_before_a'
     are 0, which results in exact time comparison only (i.e. if the times are
     not the same the disagreement value will be returned).

     A partial agreement will be calculated as follows:

       weight = agree_weight - (time_diff/(max_time_a_before_b+1)) * 
                               (agree_weight+abs(disagree_weight))

     and similar for 'max_time_a_before_b'.

     If the time difference is larger than 'max_time_a_before_b' or
     'max_time_b_before_a' the disagreement weight will be returned.

     Currently, no frequency table can be used for this comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the time arguments first, then call the base class
       constructor.
    """

    self.max_time_a_before_b = 0     # Default only exact time comparison
    self.max_time_b_before_a = 0     # Default only exact time comparison
    self.day_start =           0000  # Default value midnight

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorTime" field comparator'
      raise Exception

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_time_a_before_b'):
        if (not isinstance(value, int)) or (value < 0) or (value > 1440):
          print 'error:Argument "max_time_a_before_b" must be between 0 '+ \
          'and 1440'
          raise Exception
        self.max_time_a_before_b = value

      elif (keyword == 'max_time_b_before_a'):
        if (not isinstance(value, int)) or (value < 0) or (value > 1440):
          print 'error:Argument "max_time_b_before_a" must be between 0 '+ \
          'and 1440'
          raise Exception
        self.max_time_b_before_a = value

      elif (keyword == 'day_start'):
        if (len(value) == 5):  # Remove ':' between hours and minutes
          value = value[:2]+value[3:]
        if ((len(value) != 4) or (not value.isdigit())):
          print 'error:Illegal day start time string format or value: "%s"' % \
                (str(value))
          raise Exception

        hours =   int(value[:2])
        minutes = int(value[2:])
        if ((hours < 0) or (hours > 23) or (minutes < 0) or (minutes > 59)):
          print 'error:Illegal day start value: "%s"' % (str(value))
          raise Exception

        self.day_start = (hours * 60) + minutes

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "Time" field comparator: "%s"' % (str(self.name))
    print '1:  Fields in data set A:         %s' % (str(self.fields_a))
    print '1:  Fields in data set B:         %s' % (str(self.fields_b))
    print '1:  M-Probability:                %f' % (self.m_probability)
    print '1:  U-Probability:                %f' % (self.u_probability)
    print '2:  Agreement weight:             %f' % (self.agree_weight)
    print '2:  Disagreement weight:          %f' % (self.disagree_weight)
    print '1:  Missing weight:               %f' % (self.missing_weight)
    print '1:  Maximal time A before time B: %i' % (self.max_time_a_before_b)
    print '1:  Maximal time B before time A: %i' % (self.max_time_b_before_a)
    print '1:  Day start:                    %i' % (self.day_start)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two field lists with time values.
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Concatenate fields into two strings without whitespaces - - - - - - - - -
    #
    str_a = ''.join(fields_a)
    str_b = ''.join(fields_b)

    # Check if field values are times - - - - - - - - - - - - - - - - - - - - -
    #
    if (len(str_a) == 5):  # Remove ':' between hours and minutes
      str_a = str_a[:2]+str_a[3:]
    elif (len(str_a) != 4):
      print 'error:%s Illegal time string A format: "%s"' % \
            (records_id, str_a)
      raise Exception
    if (len(str_b) == 5):  # Remove ':' between hours and minutes
      str_b = str_b[:2]+str_b[3:]
    elif (len(str_b) != 4):
      print 'error:%s Illegal time string B format: "%s"' % \
            (records_id, str_b)
      raise Exception

    if (not str_a.isdigit()):
      print 'error:%s Illegal time string A value: "%s"' % (str_a)
      raise Exception
    if (not str_b.isdigit()):
      print 'error:%s Illegal time string B value: "%s"' % (str_b)
      raise Exception

    hours_a =   int(str_a[:2])
    minutes_a = int(str_a[2:])
    hours_b =   int(str_b[:2])
    minutes_b = int(str_b[2:])

    if (hours_a < 0) or (hours_a > 23) or (minutes_a < 0) or (minutes_a > 59):
      print 'error:%s Illegal time A value: "%s"' % (records_id, str_a)
      raise Exception
    if (hours_b < 0) or (hours_b > 23) or (minutes_b < 0) or (minutes_b > 59):
      print 'error:%s Illegal time B value: "%s"' % (records_id, str_b)
      raise Exception

    # Convert into minute values  - - - - - - - - - - - - - - - - - - - - - - -
    #
    time_a = (hours_a * 60) + minutes_a
    time_b = (hours_b * 60) + minutes_b

    # Check if times are the same - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (time_a == time_b):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    elif (self.max_time_a_before_b == 0) and (self.max_time_b_before_a == 0):
      print '3:%s    Disagreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.disagree_weight)
      return self.disagree_weight  # Because values are different

    else:

      # Convert into times according to 'day_start' value - - - - - - - - - - -
      #
      time_a -= self.day_start
      if (time_a < 0):
        time_a += 1440  # Adjust into 24-hours period
      time_b -= self.day_start
      if (time_b < 0):
        time_b += 1440  # Adjust into 24-hours period

      # Times are now in a 24-hours period (values between 0 and 1439)  - - - -
      #
      time_diff = (time_a - time_b)

      if (time_diff > 0):  # Time A after time B

        if (time_diff > self.max_time_b_before_a):  # To many minutes before
          print '3:%s    Disagreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), \
                 self.disagree_weight)
          return self.disagree_weight
        else:
          partagree_w = self.agree_weight - \
                        (float(time_diff) / (self.max_time_b_before_a+1.0)) * \
                        (self.agree_weight + abs(self.disagree_weight))
          print '3:%s    Partial agreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), partagree_w)
          return partagree_w

      else:  # Time B after time A (time_diff < 0)

        if (-time_diff > self.max_time_a_before_b):  # To many minutes before
          print '3:%s    Disagreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), \
                 self.disagree_weight)
          return self.disagree_weight

        else:
          partagree_w = self.agree_weight - \
                        (float(-time_diff) / (self.max_time_a_before_b+1.0)) *\
                        (self.agree_weight + abs(self.disagree_weight))
          print '3:%s    Partial agreement (%s,%s): %f' % \
                (records_id, str(fields_a), str(fields_b), partagree_w)
          return partagree_w

# =============================================================================

class FieldComparatorDistance(FieldComparator):
  """A field comparator that computes the geographical distance between the
     given fields. Obviously, a geocode look-up table needs to be provided.

     For field values that are not found in the geocide look-up table the
     missing value weight is returned.

     The additional arguments (besides the base class arguments) are
     'geocode_table', a reference ot the geocodelook-up table, and
     'max_distance' a positive number that gives the maximum distance (in
     kilometers) tolerated. Default value is zero.

     If the computed sitance between the two field values is smaller or equal
     to 'max_distance', the partial agreement weight is calculated using the
     following formula.

       weight = agree_weight - (distance/(max_distance+1)) * 
                               (agree_weight+abs(disagree_weight))

     Currently, no frequency table can be used for this comparator.
  """

  # Class constants
  #
  earth_radius =      6372.0  # Approximate radius of earth in kilometers
  degrees_2_radians = math.pi / 180.0

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the distance arguments first, then call the base
       class constructor.
    """

    self.geocode_table = None
    self.max_distance =  0.0  # Default only same location

    # Check if frequency table arguments are given  - - - - - - - - - - - - - -
    #
    if (('frequency_table' in kwargs.keys()) or \
        ('freq_table_min_w' in kwargs.keys()) or \
        ('freq_table_min_weight' in kwargs.keys()) or \
        ('freq_table_max_w' in kwargs.keys()) or \
        ('freq_table_max_weight' in kwargs.keys())):
      print 'error:Frequency tables can currently not be used with the'+ \
            '"FieldComparatorDistance" field comparator'
      raise Exception

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'max_distance'):
        if (not (isinstance(value, int) or isinstance(value, float))) or \
            (value < 0):
          print 'error:Argument "max_distance" must be a number equal to '+ \
                'or larger than zero'
          raise Exception
        self.max_distance = value

      elif (keyword == 'geocode_table'):
        self.geocode_table = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure geocode look-up table is set  - - - - - - - - - - - - - - - - -
    #
    if (self.geocode_table == None):
      print 'error:Geocode look-up table is not set'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "Distance" field comparator: "%s"' % (str(self.name))
    print '1:  Fields in data set A:  %s' % (str(self.fields_a))
    print '1:  Fields in data set B:  %s' % (str(self.fields_b))
    print '1:  M-Probability:         %f' % (self.m_probability)
    print '1:  U-Probability:         %f' % (self.u_probability)
    print '2:  Agreement weight:      %f' % (self.agree_weight)
    print '2:  Disagreement weight:   %f' % (self.disagree_weight)
    print '1:  Missing weight:        %f' % (self.missing_weight)
    print '1:  Geocode look-up table: %s' % (str(self.geocode_table.name))
    print '1:  Maximal distance:      %f' % (self.max_distance)

  # ---------------------------------------------------------------------------
  
  def compare(self, fields_a, fields_b, records_id):
    """Compare two field lists according to their location.

       Distance computation code based on Ole Nielson's 'point.py' and
       'distances.py' modules.
    """

    # Check if fields contain missing values  - - - - - - - - - - - - - - - - -
    #
    for field in fields_a:
      if (field in self.dataset_a.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    for field in fields_b:
      if (field in self.dataset_b.missing_values):
        print '3:%s    Missing values (%s,%s): %f' % \
              (records_id, str(fields_a), str(fields_b), self.missing_weight)
        return self.missing_weight

    # Concatenate fields into two strings without whitespaces - - - - - - - - -
    #
    str_a = ''.join(fields_a)
    str_b = ''.join(fields_b)

    # Check if field values are in geocode look-up table  - - - - - - - - - - -
    #
    loc_a = self.geocode_table.get(str_a)
    loc_b = self.geocode_table.get(str_b)

    if (loc_a in [[], '', None]) or (loc_b in [[], '', None]):
      print '3:%s    Missing values (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.missing_weight)
      return self.missing_weight  # Not both fields are in look-up table

    # Check if both locations are the same  - - - - - - - - - - - - - - - - - -
    #
    if (loc_a == loc_b):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    lon_a = loc_a[0] * self.degrees_2_radians
    lat_a = loc_a[1] * self.degrees_2_radians
    lon_b = loc_b[0] * self.degrees_2_radians
    lat_b = loc_b[1] * self.degrees_2_radians

    alpha = math.cos(lon_a - lon_b)
    x     = alpha * math.cos(lat_a)*math.cos(lat_b) + \
                    math.sin(lat_a)*math.sin(lat_b)
    dist  = self.earth_radius * math.acos(x)

    # Check if distance is zero or too far  - - - - - - - - - - - - - - - - - -
    #
    if (dist == 0.0):
      print '3:%s    Agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.agree_weight)
      return self.agree_weight

    if (dist > self.max_distance):
      print '3:%s    Disagreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), self.disagree_weight)
      return self.disagree_weight

    else:  # Compute partial agrement weight

      partagree_w = self.agree_weight - (dist / (self.max_distance+1.0)) * \
                    (self.agree_weight + abs(self.disagree_weight))
      print '3:%s    Partial agreement (%s,%s): %f' % \
            (records_id, str(fields_a), str(fields_b), partagree_w)
      return partagree_w

# =============================================================================
