# =============================================================================
# standardisation.py - Classes for cleaning and standardisations
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
# The Original Software is "standardisation.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module standardisation.py - Classes for cleaning and standardisations.

   This module provides classes for record cleaning and standardisations, 
   either based on rules or a machine learning approach (Hidden Markov models)

   TODO
   - PC 29/11/2002 Check for missing values in date fields (data standardiser)
                   And how to deal with them?
"""

# =============================================================================
# Imports go here

import os
import string
import time

import date     # Module with routines for dates
import name     # Module with name standardisation routines
import address  # Module with routines for address standardisation

# =============================================================================

def clean_component(raw_str, correction_list, record_id):
  """Clean a component (name or address) input string using a correction list.

  USAGE:
    cleaned_str = clean_component(raw_str, correction_list)

  ARGUMENTS:
    raw_str          A string containing the raw component (or parts of it)
    correction_list  A correction list as defined in 'lookup.py'

  DESCRIPTION:
    This routines cleans the input string by using the given correction list.
    It also strips off all leading and trailing spaces. A cleaned string is
    returned.
  """

  # Check if the string only contains whitespaces - - - - - - - - - - - - - - -
  #
  if (raw_str.strip() == ''):
    return ''

  # First add a trailing and leading space  - - - - - - - - - - - - - - - - - -
  # (this is to make sure replacement strings do match at beginning and end)
  #
  tmp_str = ' '+raw_str+' '

  # Check for strings in the correction-list  - - - - - - - - - - - - - - - - -
  #
  for (org,repl) in correction_list:
    tmp_str = tmp_str.replace(org,repl)

  # Make sure commas are separated from words so they become list elements  - -
  #
  tmp_str = tmp_str.replace(',', ' , ')
  tmp_str = tmp_str.strip()

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:%s    Raw input string:     "%s"' % (record_id, raw_str)
  print '3:%s    Cleaned input string: "%s"' % (record_id, tmp_str)

  return tmp_str

# =============================================================================

def check_field_spilling(str1, str2, tag_lookup_table, record_id, fields_str):
  """Routine to check if a known word is spilling from one string into another.

     Return 'True' if a word spills over, otherwise 'False'
  """

  # A list of tags that should not be checked for word spilling
  #
  no_spill_tags = ['TI','II','HY','SL','CO','VB','RU','PC']

  org_str1 = str1  # Keep copies of the original input string
  org_str2 = str2

  alphanum = string.letters+string.digits  # String with all letters and digits

  if (str1 == ''):
    return False

  # Check if the characters at the 'boundary' are either letters or digits
  #
  elif (str1[-1] in alphanum) and (str2[0] in alphanum):
    test_str1 = ''
    while (str1 != '') and (str1[-1] in alphanum): # Get last char. from str1
      test_str1 =  str1[-1] + test_str1
      str1 = str1[:-1]
    test_str2 = ''
    while (str2 != '') and (str2[0] in alphanum): # Get first char. from str2
      test_str2 = test_str2 + str2[0]
      str2 = str2[1:]

    check_word = test_str1 + test_str2

    if (tag_lookup_table.has_key((check_word,))):
      list_entry = tag_lookup_table[(check_word,)]

      entry_tags = list_entry[-1]  # Get tags for this entry

      do_spilling = True
      while entry_tags != '':  # Check if the tag(s) are in the 'no_spill_tags'

        if (entry_tags[:2] in no_spill_tags):
          do_spilling = False
        entry_tags = entry_tags[3:]  # Remove first tag and '/'

      if (do_spilling == True):
        print '2:%s    Found (and corrected) word  spilling:' % (record_id) + \
              ' "%s","%s" -> "%s"%s' % (org_str1.strip(), org_str2.strip(), \
              check_word, fields_str)
        return True  # Found a word spilling
      else:
        return False

  else:
    return False

# =============================================================================

class RecordStandardiser:
  """Main class for the standardisation process. Implements routines for record
     cleaning, tagging and standardisations, both block wise and for single
     records.

     A record standardiser is made of component standardisers, which can either
     be standardisers for names, addresses or dates. These component
     standardisers need to be given when a record standarsider is initialised.

     Each component standardiser needs the following attributes:
       input_dataset   Reference to the input data set.
       output_dataset  Reference to the output data set.
       input_fields    A list of one or more field names from the input data
                       set.
       output_fields   A list of field names (depending upon what standardiser)
                       from the output data set.

     A method 'standard' is needed, which takes as inputs one or more fields
     (strings) and returns a number of fields containing the standardised
     values (for example, a date standardiser will return three fields
     [day, month, year]).

     The attributes 'input_dataset' and 'output_dataset' will be set by the
     record standardiser, so they don't need to be set by the component
     standardiser itself.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.name =                    ''
    self.description =             ''
    self.input_dataset =           None
    self.output_dataset =          None
    self.component_standardisers = None  # A list of one or more standardisers
                                         # for names, addresses or dates.

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'input_dataset'):
        self.input_dataset = value
      elif (keyword == 'output_dataset'):
        self.output_dataset = value

      elif (keyword in ['comp_std','component_standardisers']):
        if (not isinstance(value, list)):
          print 'error:Argument "component_standardisers" is not a list'
          raise Exception
        self.component_standardisers = value

      else:
        print 'error:Illegal constructor argument keyword: '+keyword
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (self.input_dataset == None):
      print 'error:Input data set is not defined'
      raise Exception

    if (self.output_dataset == None):
      print 'error:Output data set is not defined'
      raise Exception

    if (self.component_standardisers == None):
      print'error:No component standardiser defined'
      raise Exception

    # Set the data set attributes in the component standardisers and check  - -
    # if the fields are proper fields of the data sets.
    #
    for cs in self.component_standardisers:
      cs.input_dataset =  self.input_dataset
      cs.output_dataset = self.output_dataset

      for field in cs.input_fields:
        if (not self.input_dataset.fields.has_key(field)):
          print 'error:Input data set "%s" does not have a field "%s"' % \
                (str(self.input_dataset.name), str(field))
          raise Exception

      for field in cs.output_fields:
        if (field != None):  # Only check fields that are not set to None
          if (not self.output_dataset.fields.has_key(field)):
            print 'error:Output data set "%s" does not have a field "%s"' % \
                  (str(self.output_dataset.name), str(field))
            raise Exception

    # Check if there is no conflict in the output fields definition - - - - - -
    #
    output_fields_dict = {}

    for cs in self.component_standardisers:

      for field in cs.output_fields:
        if (field != None):  # Only check fields that are not set to None
          if (output_fields_dict.has_key(field)):
            print 'error:Output fields definition conflict with field "%s"' % \
                  (str(field))
            raise Exception
          output_fields_dict[field] = 1

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised record standardiser'
    print '1:  Input data set:  %s' % (str(self.input_dataset.name))
    print '1:  Output data set: %s' % (str(self.output_dataset.name))
    print '1:  Component standardisers:'
    for cs in self.component_standardisers:
      print '1:    Name: %s' % (str(cs.name))
      print '1:      Input fields:  %s' % (str(cs.input_fields))
      print '1:      Output fields: %s' % (str(cs.output_fields))

  # ---------------------------------------------------------------------------

  def standardise(self, record):
    """Standardise a record.
    """

    if (not isinstance(record, dict)):
      print 'error:Input record is not a dictionary: %s' % (str(record))
      raise Exception

    output_record = {}

    record_id = '[RecID: %s/%s]' % \
                (str(record['_rec_num_']), str(record['_dataset_name_']))

    # Standardise the input fields by passing them on to the component  - - - -
    # standardisers
    #
    for cs in self.component_standardisers:

      # Extract field values from input data set
      #
      fields = []

      for field in cs.input_fields:
        fields.append(record.get(field, cs.input_dataset.fields_default))

      # A string with the input fields, for logging
      #
      fields_str = os.linesep+'[Fields: %s]' % (str(fields))

      output_record.update(cs.standardise(fields, record_id, fields_str))

      # Insert record identification fields
      #
      output_record['_rec_num_'] =      record['_rec_num_']
      output_record['_dataset_name_'] = record['_dataset_name_']

    # A log message for medium volume log output (level 2)  - - - - - - - - - -
    #
    print '2:%s  Record standardisation:' % (record_id)
    print '2:%s    Input record:  %s' % (str(record), record_id)
    print '2:%s    Output record: %s' % (str(output_record), record_id)

    return output_record

  # ---------------------------------------------------------------------------

  def standardise_block(self, record_list):
    """Standardise a block of records.
    """

    if (not isinstance(record_list, list)):
      print 'error:Input record list is not a list'
      raise Exception

    output_list = []

    for record in record_list:

      output_record = self.standardise(record)

      output_list.append(output_record)

    return output_list

# =============================================================================

class ComponentStandardiser:
  """Base class for component standardisers.

     The following arguments must be given to the constructor of the base class
     and all derived classes:
       input_fields
       output_fields

     The following optional arguments can be used with all component
     standardisers:
       name
       description
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all component standardisers.
    #
    self.name =                  ''
    self.description =           ''
    self.input_dataset =  None  # Reference to the input data set, will be set
                                # by the record standardiser.
    self.output_dataset = None  # Reference to the output data set, will be set
                                # by the record standardiser.
    self.input_fields = None    # A string with a field name or a list of field
                                # names form the input data set, that will be
                                # standardised (if a list is given they will be
                                # concatenated into one string).
    self.output_fields = None   # A string with a field name or a list of field
                                # names form the output data set.

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value

      elif (keyword == 'input_fields'):
        if (isinstance(value, str)):
          value = [value]
        if (not isinstance(value, list)):
          print 'error:Argument "input_fields" is not a string or a ' + \
                'list: %s' % (str(value))
          raise Exception
        self.input_fields = value

      elif (keyword == 'output_fields'):
        if (isinstance(value, str)):
          value = [value]
        if (not isinstance(value, list)):
          print 'error:Argument "output_fields" is not a string or a ' + \
                'list: %s' % (str(value))
          raise Exception
        self.output_fields = value

      else:
        print 'error:Illegal constructor argument keyword: %s' % (str(keyword))
        raise Exception

    # Check if fields are defined and not empty - - - - - - - - - - - - - - - -
    #
    if (self.input_fields == None) or (self.output_fields == None):
      print 'error:Input or output fields are not defined'
      raise Exception

    if (len(self.input_fields) == 0) or (len(self.output_fields) == 0):
      print 'error:Fields must not be empty lists'
      raise Exception

    output_field_def = 0
    for f in self.output_fields:
      if (f != None):
        output_field_def = 1  # An output field is defined (not set to None)
        break
    if (output_field_def == 0):
      print 'error:At least one output field must be defined'
      raise Exception

  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Standardise the component defined in the input fields for the given
       record.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

# =============================================================================

class PassFieldStandardiser(ComponentStandardiser):
  """Dummy standardiser used to simply pass fields from input to output data
     set without doing any standardisation. Values are simply copied directly
     from the input field(s) to the output field(s).

     The number of input fields must be the same as the number of output fields
     as values are copied directly from an input field to the corresponding
     output field.

     No additional argument besides the base class arguments can be given to
     this field standardiser.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    ComponentStandardiser.__init__(self, kwargs)  # Initialise base class

    # Check if the number of input fields equals to the number of output fields
    #
    if (len(self.input_fields) != len(self.output_fields)):
      print 'error:Different number of input and output fields for ' + \
            '"PassField" standardiser'
      raise Exception

    self.num_fields = len(self.input_fields)  # Save number of fields

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "PassField" standardiser: "%s"' % \
          (str(self.name))
    print '1:  Input fields:            %s' % (str(self.input_fields))
    print '1:  Output fields:           %s' % (str(self.output_fields))
 
  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Copy the values from the input fields into corresponding output fields.
    """

    result = {}  # Output results dictionary

    if (self.num_fields == 1) and (isinstance(fields, str)):
      fields = [fields]

    if (not isinstance(fields, list)):
      print 'error:Input fields are not a list or a string: %s' % (str(fields))
      raise Exception

    if (self.num_fields != len(fields)):
      print 'error:Wrong number of input fields: %i (should be %i)' % \
            (len(fields), self.num_fields)
      raise Exception

    for i in range(self.num_fields):
      if (fields[i].strip() != ''):  # Only copy non-empty fields
        result[self.output_fields[i]] = fields[i]  # Copy into output fields

    return result

# =============================================================================

class DateStandardiser(ComponentStandardiser):
  """Routines for standardising dates into tuples (day,month,year).

     The 'output_fields' must be a list of three fields, the first field is for
     day, the second for month and the third for year. Fields can be set to
     'None' if no output is to be written, as long as at least one field is
     set.

     The additional arguments (besides the base class arguments) are
       parse_formats  A string with a date parse format or a list with date
                      parsing format strings.
       pivot_year     Value of pivot year (between 00 and 99) that controls
                      expansion of two-digit year values into four-digit year
                      values. Two-digits years smaller than the pivot year will
                      be expanded into 20XX, years larger and equal than the
                      pivot year will be expanded into 19xx
                      For example: pivot_year = 03:  68 -> 1968
                                                     03 -> 1903
                                                     02 -> 2002
                      The default value is the current year plus 1
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.parse_formats = None
    self.pivot_year =    None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'parse_formats'):
        if (isinstance(value, str)):
          value = [value]
        if (not isinstance(value, list)):
          print 'error:Argument "parse_formats" is not a string or a list'
          raise Exception
        self.parse_formats = value

      elif (keyword == 'pivot_year'):
        if (not isinstance(value, int)) or (value < 0) or (value > 99):
          print 'error:Argument "pivot_year" is not an integer number '+ \
                'between 0 and 99: %s' % (str(value))
          raise Exception
        self.pivot_year = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.output_fields) != 3):
      print 'error:Attribute "output_fields" is not a list with three '+ \
            'elements: %s' % (str(self.output_fields))
      raise Exception

    if (self.parse_formats == None):
      print 'error:Date parse format strings not defined'
      raise Exception

    # Check if parse format strings are in a valid form - - - - - - - - - - - -
    #
    for format_str in self.parse_formats:
      if (not isinstance(format_str, str)):
        print 'error:Format string "%s" is not a string' % (str(format_str))
        raise Exception
      tmp_str = format_str

      # Check if beginning of format string is a valid directive
      #
      if (tmp_str[:2] not in ['%b','%B','%d','%m','%M','%y','%Y']):
        print 'error:Illegal first directive in format string "%s"' % \
              (format_str)
        raise Exception

      tmp_str = tmp_str[2:]  # Remove first diective

      if (tmp_str[0] == ' '):
        tmp_str = tmp_str[1:]  # Remove space between directive

      # Check if middle of format string is a valid directive
      #
      if (tmp_str[:2] not in ['%b','%B','%d','%m','%M','%y','%Y']):
        print 'error:Illegal second directive in format string "%s"' % \
              (format_str)
        raise Exception

      tmp_str = tmp_str[2:]  # Remove second diective

      if (tmp_str[0] == ' '):
        tmp_str = tmp_str[1:]  # Remove space between directive

      # Check if end of format string is a valid directive
      #
      if (tmp_str[:2] not in ['%b','%B','%d','%m','%M','%y','%Y']):
        print 'error:Illegal third directive in format string "%s"' % \
              (format_str)
        raise Exception

      # Check if there is nothing after the directive
      #
      if (len(tmp_str) != 2):
        print 'error:Illegal format string "%s"' % (format_str)
        raise Exception

    # Save the number of date parse format strings  - - - - - - - - - - - - - -
    #
    self.num_date_formats = len(self.parse_formats)

    # Check if pivot year is set, if not set to current year plus one - - - - -
    #
    if (self.pivot_year == None):
      today = date.get_today()  # Get today's date

      next_year = (int(today[2]) % 100) + 1  # Get next year as integer number

      if (next_year < 0) or (next_year > 99):
        print 'error:Illegal value for "next_year" (not between 0 and 99): '+ \
              str(next_year)
        raise Exception

      self.pivot_year = next_year

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "Date" component standardiser: "%s"' % \
          (str(self.name))
    print '1:  Input fields:            %s' % (str(self.input_fields))
    print '1:  Output fields:           %s' % (str(self.output_fields))
    print '1:  Number of parse formats: %s' % (str(self.num_date_formats))
    print '1:  Pivot year:              %s' % (str(self.pivot_year))

  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Standardise the date defined in the input fields for the given record.
    """

    # Concatenate into a strings without whitespaces  - - - - - - - - - - - - -
    #
    date_str = ''.join(fields)

    if (date_str.strip() == ''):  # No date given
      return {}

    parsed = 0  # Flag for successful date parsing
    iter =   0  # Iteration counter

    # Try one date format after the other until success or none worked  - - - -
    #
    while (parsed != 1) and (iter < self.num_date_formats):

      date_try = date.str_to_date(date_str, self.parse_formats[iter],
                                  self.pivot_year)
      if (len(date_try) == 3):
        parsed = 1  # Parsed successful
      else:
        iter += 1  # Try next format string

    if (parsed == 1):  # Date successfully parsed - - - - - - - - - - - - - - -
      day =   date_try[0]
      month = date_try[1]
      year =  date_try[2]

      # A log message for high volume log output (level 3)  - - - - - - - - - -
      #
      print '3:%s    Date string: "%s" parsed with ' % (record_id, date_str)+ \
            'format string "%s" and pivot year %s: %s' % \
            (self.parse_formats[iter], str(self.pivot_year), \
            str([day,month,year]))
    else:

      print 'warning:%s Could not parse date string "%s"%s' % \
            (record_id, date_str, fields_str)

      # The input fields could not be parsed into a valid date, so set  - - - -
      # to missing values (first missing value in the output dat set's
      # missing value list)
      #
      return {}

    # Now copy result into defined (not None) output fields - - - - - - - - - -
    #
    result = {}

    if (self.output_fields[0] != None):  # Day is defined
      result[self.output_fields[0]] = day
    else:
      print '3:%s    Throw away "day" value: %s' % (record_id, str(day))

    if (self.output_fields[1] != None):  # Month is defined
      result[self.output_fields[1]] = month
    else:
      print '3:%s    Throw away "month" value: %s' % (record_id, str(month))

    if (self.output_fields[2] != None):  # Year is defined
      result[self.output_fields[2]] = year
    else:
      print '3:%s    Throw away "year" value: %s' % (record_id, str(year))

    print '2:%s  Input date string "%s" parsed into date: %s' % \
          (record_id, date_str, str(result))

    return result

# =============================================================================

class NameRulesStandardiser(ComponentStandardiser):
  """Routines for standardising names into tuples (title, gender guess, given
     names, alternative given names, surnames, alternative surnames) using
     rules.

     The 'output_fields' must be a list of six field names from the output data
     set into which the standardised names will be written. The fields are
       title
       gender guess
       given names
       alternative given names
       surnames
       alternative surnames
     Fields can be set to 'None' if no output is to be written, as long as at
     least one field is set.

     The additional arguments (besides the base class arguments) are
       name_corr_list    Reference to a correction list for names
       name_tag_table    Reference to a name tag lookup table
       male_titles       A list of male titles (like 'mr'), used to guess a
                         gender
       female_titles     A list of female titles (like 'ms'), used to guess a
                         gender
       first_name_comp   The expected first component in the namse, either
                         'gname' (for given names) or 'sname' (for surnames)
                         Default value is 'gname'
       field_separator   A character that is used to join together fields in
                         the name component in the cleaning routine. Default is
                         a whitespace.
       check_word_spill  A flag (can be 'True' or 'False'), if set to true
                         word 'spilling' between fields is be checked using the
                         name tag look-up table  (e.g. "peter pa","ul miller"
                         will become "peter paul miller").
                         This is only useful if the 'field_separator' is not an
                         empty string.
                         Default for word spill checking is 'True'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.name_corr_list =   None
    self.name_tag_table =   None
    self.male_titles =      []
    self.female_titles =    []
    self.first_name_comp =  'gname'
    self.field_separator =  ' '
    self.check_word_spill = True

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'male_titles'):
        if (not isinstance(value, list)):
          print 'error:Argument "male_titles" is not a list'
          raise Exception
        self.male_titles = value

      elif (keyword == 'female_titles'):
        if (not isinstance(value, list)):
          print 'error:Argument "female_titles" is not a list'
          raise Exception
        self.female_titles = value

      elif (keyword == 'first_name_comp'):
        if (value not in ['gname','sname']):
          print 'error:Argument "first_name_comp" must be either "gname" '+ \
                'or "sname"'
          raise Exception
        self.first_name_comp = value

      elif (keyword == 'field_separator'):
        if (not isinstance(value, str)):
          print 'error:Argument "field_separator" is not a string'
          raise Exception
        self.field_separator = value

      elif (keyword == 'check_word_spill'):
        if (value not in [True,False]):
          print 'error:Argument "check_word_spill" must be "True" or "False"'
          raise Exception
        self.check_word_spill = value
 
      elif (keyword in ['name_corr_list','name_list','corr_list']):
        self.name_corr_list = value
      elif (keyword in ['name_tag_table','tag_table','name_tag']):
        self.name_tag_table = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.output_fields) != 6):
      print 'error:Attribute "output_fields" is not a list with six '+ \
            'elements: %s' % (str(self.output_fields))
      raise Exception

    if (self.name_corr_list == None):
      print 'error:Name correction list is not defined'
      raise Exception

    if (self.name_tag_table == None):
      print 'error:Name tag table is not defined'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "NameRules" component standardiser: "%s"' % \
          (str(self.name))
    print '1:  Input fields:           %s' % (str(self.input_fields))
    print '1:  Output fields:          %s' % (str(self.output_fields))
    print '1:  Male title list:        %s' % (str(self.male_titles))
    print '1:  Female title list:      %s' % (str(self.female_titles))
    print '1:  First name component:   %s' % (str(self.first_name_comp))
    print '1:  Name correction list:   %s' % (str(self.name_corr_list.name))
    print '1:  Name tag look-up table: %s' % (str(self.name_tag_table.name))
    print '1:  Field separator:        "%s"' % (self.field_separator)
    print '1:  Check word spilling:    %s' % (str(self.check_word_spill))

  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Standardise the name defined in the input fields for the given record.
    """

    # Concate fields into one string  - - - - - - - - - - - - - - - - - - - - -
    #
    name_str = ''

    for f in fields:

      if (f != ''):
        field_sep = self.field_separator

        # Check field spilling only if field separator is not an empty string
        #
        if (field_sep != '') and (self.check_word_spill == True):

          spill_flag = check_field_spilling(name_str, f, self.name_tag_table,
                                            record_id, fields_str)
          if (spill_flag == True):
            field_sep = ''  # A word spills over, so make field separator an
                            # empty string

        name_str = name_str+field_sep+f  # Append separator and field

    # Clean the name string - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    name_str = clean_component(name_str, self.name_corr_list, record_id)

    if (name_str == ''):  # Cleaned name string is empty
      return {}

    # Tag the name string and split into elements - - - - - - - - - - - - - - -
    #
    [name_list, tag_list] = name.tag_name_component(name_str, \
                                                    self.name_tag_table,
                                                    record_id)

    gender_guess = name.get_gender_guess(name_list, tag_list, \
                                         self.male_titles, self.female_titles,
                                         record_id)

    [name_list, tag_list, title_list] = name.get_title(name_list, tag_list,
                                                       record_id)

    names_list = name.get_name_rules(name_list, tag_list, self.first_name_comp,
                                     record_id, fields_str)

    # Now copy result into defined (not None) output fields - - - - - - - - - -
    #
    result = {}

    if (self.output_fields[0] != None):  # Title is defined
      tmp_value = string.join(title_list,' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[0]] = tmp_value
    else:
      print '3:%s    Throw away "title" value: %s' % \
            (record_id, string.join(title_list,' '))

    if (self.output_fields[1] != None):  # Gender guess is defined
      tmp_value = gender_guess.strip()
      if (tmp_value != ''):
        result[self.output_fields[1]] = tmp_value
    else:
      print '3:%s    Throw away "gender guess" value: %s' % \
            (record_id, gender_guess)

    if (self.output_fields[2] != None):  # Given name is defined
      tmp_value = string.join(names_list[0],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[2]] = tmp_value
    else:
      print '3:%s    Throw away "given name" value: %s' % \
            (record_id, string.join(names_list[0],' '))

    if (self.output_fields[3] != None):  # Alternative given name is defined
      tmp_value = string.join(names_list[1],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[3]] = tmp_value
    else:
      print '3:%s    Throw away "alternative given name" value: %s' % \
            (record_id, string.join(names_list[1],' '))

    if (self.output_fields[4] != None):  # Surname is defined
      tmp_value = string.join(names_list[2],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[4]] = tmp_value
    else:
      print '3:%s    Throw away "surname" value: %s' % \
            (record_id, string.join(names_list[2],' '))

    if (self.output_fields[5] != None):  # Alternative surname is defined
      tmp_value = string.join(names_list[3],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[5]] = tmp_value
    else:
      print '3:%s    Throw away "alternative surname" value: %s'+ \
            (record_id, string.join(names_list[3],' '))

    print '2:%s  Input name string "%s" standardised into: %s' %\
          (record_id, name_str, str(result))

    return result

# =============================================================================

class NameHMMStandardiser(ComponentStandardiser):
  """Routines for standardising names into tuples (title, gender guess, given
     names, alternative given names, surnames, alternative surnames) using
     a hidden Markov model (HMM) approach.

     The 'output_fields' must be a list of six field names from the output data
     set into which the standardised names will be written. The fields are
       title
       gender guess
       given names
       alternative given names
       surnames
       alternative surnames
     Fields can be set to 'None' if no output is to be written, as long as at
     least one field is set.

     The additional arguments (besides the base class arguments) are
       name_corr_list    Reference to a correction list for names
       name_tag_table    Reference to a name tag lookup table
       male_titles       A list of male titles (like 'mr'), used to guess a
                         gender.
       female_titles     A list of female titles (like 'ms'), used to guess a
                         gender.
       first_name_comp   The expected first component in the namse, either
                         'gname' (for given names) or 'sname' (for surnames)
                         Default value is 'gname'
       name_hmm          Reference to a name hidden Markov model (HMM)
       field_separator   A character that is used to join together fields in
                         the name component in the cleaning routine. Default is
                         a whitespace.
       check_word_spill  A flag (can be 'True' or 'False'), if set to true
                         word 'spilling' between fields is be checked using the
                         name tag look-up table (e.g. "peter pa","ul miller"
                         will become "peter paul miller").
                         This is only useful if the 'field_separator' is not an
                         empty string.
                         Default for word spill checking is 'True'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.name_corr_list =   None
    self.name_tag_table =   None
    self.name_hmm =         None
    self.male_titles =      []
    self.female_titles =    []
    self.first_name_comp =  'gname'
    self.field_separator =  ' '
    self.check_word_spill =  True

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword == 'male_titles'):
        if (not isinstance(value, list)):
          print 'error:Argument "male_titles" is not a list'
          raise Exception
        self.male_titles = value

      elif (keyword == 'female_titles'):
        if (not isinstance(value, list)):
          print 'error:Argument "female_titles" is not a list'
          raise Exception
        self.female_titles = value

      elif (keyword == 'first_name_comp'):
        if (value not in ['gname','sname']):
          print 'error:Argument "first_name_comp" must be either "gname" '+ \
                'or "sname"'
          raise Exception
        self.first_name_comp = value

      elif (keyword == 'field_separator'):
        if (not isinstance(value, str)):
          print 'error:Argument "field_separator" is not a string'
          raise Exception
        self.field_separator = value

      elif (keyword == 'check_word_spill'):
        if (value not in [True,False]):
          print 'error:Argument "check_word_spill" must be "True" or "False"'
          raise Exception
        self.check_word_spill = value
 
      elif (keyword in ['name_corr_list','name_list','corr_list']):
        self.name_corr_list = value
      elif (keyword in ['name_tag_table','tag_table','name_tag']):
        self.name_tag_table = value

      elif (keyword == 'name_hmm'):
        self.name_hmm = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.output_fields) != 6):
      print 'error:Attribute "output_fields" is not a list with six '+ \
            'elements: %s' % (str(self.output_fields))
      raise Exception

    if (self.name_corr_list == None):
      print 'error:Name correction list is not defined'
      raise Exception

    if (self.name_tag_table == None):
      print 'error:Name tag table is not defined'
      raise Exception

    if (self.name_hmm == None):
      print 'error:Name hidden Markov model (HMM) is not defined'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "NameHMM" component standardiser: "%s"' % \
          (str(self.name))
    print '1:  Input fields:             %s' % (str(self.input_fields))
    print '1:  Output fields:            %s' % (str(self.output_fields))
    print '1:  Male title list:          %s' % (str(self.male_titles))
    print '1:  Female title list:        %s' % (str(self.female_titles))
    print '1:  First name component:     %s' % (str(self.first_name_comp))
    print '1:  Name correction list:     %s' % (str(self.name_corr_list.name))
    print '1:  Name tag look-up table:   %s' % (str(self.name_tag_table.name))
    print '1:  Field separator:          "%s"' % (self.field_separator)
    print '1:  Check word spilling:      %s' % (str(self.check_word_spill))
    print '1:  Name hidden Markov model: %s' % (str(self.name_hmm.name))

  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Standardise the name defined in the input fields for the given record.
    """

    # Concate fields into one string  - - - - - - - - - - - - - - - - - - - - -
    #
    name_str = ''

    for f in fields:

      if (f != ''):
        field_sep = self.field_separator

        # Check field spilling only if field separator is not an empty string
        #
        if (field_sep != '') and (self.check_word_spill == True):

          spill_flag = check_field_spilling(name_str, f, self.name_tag_table,
                                            record_id, fields_str)
          if (spill_flag == True):
            field_sep = ''  # A word spills over, so make field separator an
                            # empty string

        name_str = name_str+field_sep+f  # Append separator and field

    # Clean the name string - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    name_str = clean_component(name_str, self.name_corr_list, record_id)

    if (name_str == ''):  # Cleaned name string is empty
      return {}

    # Tag the name string and split into elements - - - - - - - - - - - - - - -
    #
    [name_list, tag_list] = name.tag_name_component(name_str, \
                                                    self.name_tag_table,
                                                    record_id)

    gender_guess = name.get_gender_guess(name_list, tag_list, \
                                         self.male_titles, self.female_titles,
                                         record_id)

    [name_list, tag_list, title_list] = name.get_title(name_list, tag_list,
                                                       record_id)

    names_list = name.get_name_hmm(name_list, tag_list, self.name_hmm,
                                   self.first_name_comp, self.name_tag_table,
                                   record_id, fields_str)

    # names_list[0] would contain title words, but they are already processed!

    # Now copy result into defined (not None) output fields - - - - - - - - - -
    #
    result = {}

    if (self.output_fields[0] != None):  # Title is defined
      tmp_value = string.join(title_list,' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[0]] = tmp_value
    else:
      print '3:%s    Throw away "title" value: %s' % \
            (record_id, string.join(title_list,' '))

    if (self.output_fields[1] != None):  # Gender guess is defined
      tmp_value = gender_guess.strip()
      if (tmp_value != ''):
        result[self.output_fields[1]] = tmp_value
    else:
      print '3:%s    Throw away "gender guess" value: %s' % \
            (record_id, gender_guess)

    if (self.output_fields[2] != None):  # Given name is defined
      tmp_value = string.join(names_list[1],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[2]] = tmp_value
    else:
      print '3:%s    Throw away "given name" value: %s' % \
            (record_id, string.join(names_list[0],' '))

    if (self.output_fields[3] != None):  # Alternative given name is defined
      tmp_value = string.join(names_list[2],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[3]] = tmp_value
    else:
      print '3:%s    Throw away "alternative given name" value: %s' % \
            (record_id, string.join(names_list[1],' '))

    if (self.output_fields[4] != None):  # Surname is defined
      tmp_value = string.join(names_list[3],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[4]] = tmp_value
    else:
      print '3:%s    Throw away "surname" value: %s' % \
            (record_id, string.join(names_list[2],' '))

    if (self.output_fields[5] != None):  # Alternative surname is defined
      tmp_value = string.join(names_list[4],' ').strip()
      if (tmp_value != ''):
        result[self.output_fields[5]] = tmp_value
    else:
      print '3:%s    Throw away "alternative surname" value: %s'+ \
            (record_id, string.join(names_list[3],' '))

    print '2:%s  Input name string "%s" standardised into: %s' %\
          (record_id, name_str, str(result))

    return result

# =============================================================================

class AddressHMMStandardiser(ComponentStandardiser):
  """Routines for standardising addresses into output field tuples  using a
     hidden Markov model (HMM) approach.

     The 'output_fields' must be a list of 17 field names from the output data
     set into which the standardised names will be written. The fields are
      wayfare_number
      wayfare_name
      wayfare_qualifier
      wayfare_type
      unit_number
      unit_type
      property_name
      institution_name
      institution_type    
      postaddress_number
      postaddress_type
      locality_name
      locality_qualifier
      postcode
      territory
      country
      address_hmm_prob
     Fields can be set to 'None' if no output is to be written, as long as at
     least one field is set.

     The additional arguments (besides the base class arguments) are
       address_corr_list   Reference to a correction list for addresses
       address_tag_table   Reference to an address tag lookup table
       address_hmm         Reference to an address hidden Markov model (HMM)
       field_separator     A character that is used to join together fields in
                           the name component in the cleaning routine. Default
                           is a whitespace.
       check_word_spill    A flag (can be 'True' or 'False'), if set to true
                           word 'spilling' between fields is be checked using
                           the address tag look-up table
                           (e.g. "miller r","d canberra" will become
                           "miller rd canberra").
                           This is only useful if the 'field_separator' is not
                           an empty string.
                           Default for word spill checking is 'True'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    # Initialise attributes
    #
    self.address_corr_list = None
    self.address_tag_table = None
    self.address_hmm =       None
    self.field_separator =  ' '
    self.check_word_spill = True

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():
      if (keyword in ['address_corr_list','address_list','corr_list']):
        self.address_corr_list = value
      elif (keyword in ['address_tag_table','tag_table','address_tag']):
        self.address_tag_table = value

      elif (keyword == 'field_separator'):
        if (not isinstance(value, str)):
          print 'error:Argument "field_separator" is not a string'
          raise Exception
        self.field_separator = value

      elif (keyword == 'check_word_spill'):
        if (value not in [True,False]):
          print 'error:Argument "check_word_spill" must be "True" or "False"'
          raise Exception
        self.check_word_spill = value
 
      elif (keyword == 'address_hmm'):
        self.address_hmm = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.output_fields) != 17):
      print 'error:Attribute "output_fields" is not a list with 17 '+ \
            'elements: %s' % (str(self.output_fields))
      raise Exception

    if (self.address_corr_list == None):
      print 'error:Address correction list is not defined'
      raise Exception

    if (self.address_tag_table == None):
      print 'error:Address tag table is not defined'
      raise Exception

    if (self.address_hmm == None):
      print 'error:Address hidden Markov model (HMM) is not defined'
      raise Exception

    # A log message for low/medium volume log output (level 1/2)  - - - - - - -
    #
    print '1:'
    print '1:Initialised "AddressHMM" component standardiser: "%s"' % \
          (str(self.name))
    print '1:  Input fields:               %s' % (str(self.input_fields))
    print '1:  Output fields:              %s' % (str(self.output_fields))
    print '1:  Address correction list:    %s' % \
          (str(self.address_corr_list.name))
    print '1:  Address tag look-up table:  %s' % \
          (str(self.address_tag_table.name))
    print '1:  Address hidden Markov mdel: %s' % \
          (str(self.address_hmm.name))
    print '1:  Field separator:            "%s"' % (self.field_separator)
    print '1:  Check word spilling:        %s' % (str(self.check_word_spill))

  # ---------------------------------------------------------------------------

  def standardise(self, fields, record_id, fields_str):
    """Standardise the address defined in input fields for the given record.
    """

    # Concate fields into one string  - - - - - - - - - - - - - - - - - - - - -
    #
    address_str = ''

    for f in fields:

      if (f != ''):
        field_sep = self.field_separator

        # Check field spilling only if field separator is not an empty string
        #
        if (field_sep != '') and (self.check_word_spill == True):

          spill_flag = check_field_spilling(address_str, f,
                                            self.address_tag_table, record_id,
                                            fields_str)
          if (spill_flag == True):
            field_sep = ''  # A word spills over, so make field separator an
                            # empty string

        address_str = address_str+field_sep+f  # Append separator and field

    # Clean the address string  - - - - - - - - - - - - - - - - - - - - - - - -
    #
    address_str = clean_component(address_str, self.address_corr_list,
                                  record_id)

    if (address_str == ''):  # Cleaned address string is empty
      return {}

    # Tag the address string and split into elements  - - - - - - - - - - - - -
    #
    [address_list, tag_list] = address.tag_address_component(address_str,
                                                        self.address_tag_table,
                                                        record_id)

    address_dict = address.get_address_hmm(address_list, tag_list,
                                           self.address_hmm,
                                           self.address_tag_table,
                                           record_id, fields_str)

    # Now copy result into defined (not None) output fields - - - - - - - - - -
    #
    result = {}

    for (key,val) in address_dict.items():
      if (key in self.output_fields):
        tmp_value = string.join(val,' ').strip()
        if (tmp_value != ''):
          result[key] = tmp_value
      else:
        print '3:%s    Throw away field "%s" with value: %s' % \
              (record_id, str(key), string.join(val,' '))

    print '2:%s  Input address string "%s" standardised into: %s' % \
          (record_id, address_str, str(result))

    return result

# =============================================================================
