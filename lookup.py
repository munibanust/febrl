# =============================================================================
# lookup.py - Classes for various types of look-up tables.
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
# The Original Software is "lookup.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module lookup.py - Classes for various types of look-up tables.

   This module contains classes for look-up table and correction lists.
"""

# =============================================================================
# Imports go here

import string
import types

# =============================================================================

class LookupTable(dict):
  """class LookupTable - Based on dictionary type.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor, set general attributes.
    """

    dict.__init__(self)  # Initialise dictionary base type
    self.name =         ''
    self.description =  ''
    self.created =      ''
    self.modified =     ''
    self.file_names =   []
    self.default =      None  # Default return value for non existing keys
    self.length =       None  # Number of entries in the look-up table

    for (keyword, value) in kwargs.items():
      if (keyword == 'default'):
        self.default = value
      elif (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value
      elif (keyword == 'created'):
        self.created = value
      elif (keyword == 'modified'):
        self.modified = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

  # ---------------------------------------------------------------------------

  def __getitem__(self, key):
    """Return an item in the look-up table with the given key. If not found,
       return the default value.
    """

    try:
      return dict.__getitem__(self, key)
    except KeyError:
      return self.default

  # ---------------------------------------------------------------------------

  def get(self, key, *args):
    """Return an item in the look-up table with the given key. If not found,
       return the default value.
    """

    if (not args):
      args = (self.default,)
    return dict.get(self, key, *args)

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files into the look-up table.
       See implementations in derived classes for details.
    """

    print 'error:Override abstract method in derived class'
    raise Exception

# =============================================================================

class TagLookupTable(LookupTable):
  """class TagLookupTable - Look-up tables for word corrections and tags.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
 
    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.max_key_length = None  # The maximum length of a key in words

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised tag look-up table "%s"' %(str(self.name))
    print '1:  With default: "%s"' % (str(self.default))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with word corrections and tags into the
       look-up table.
       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    if (not isinstance(file_names, list)):
      print 'error:Argument "file_names" is of wrong type, must be either '+ \
            'string or list'
      raise Exception

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table
    self.max_key_length = 0

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      # Open file and read all lines into a list
      #
      try:
        f = open(fn,'r')
      except:
        print 'error:Can not read from file "%s"' % (fn)
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      tag = ''  # Start with no tag

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()  # Remove line separators
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          if (l[:5] == 'tag=<'):  # It's a line with a new tag
            tag = l[5:7]

          else:  # A line with an entry

            # Make sure a tag is set
            #
            if (tag == ''):
              print 'error:Missing tag specification in file "%s"' % (fn)
              raise Exception

            ll = l.split(':')  # Separate key from values

            if (len(ll) == 2):  # Line contains a key - - - - - - - - - - - - -
              k = ll[0].strip().lower()  # Get key, make lower and strip spaces

              k_list = k.split(' ')  # Make a list of key words
              if (len(k_list) > self.max_key_length):
                self.max_key_length = len(k_list) # Update maximal key length

              val = string.join(k_list,'_')
              key = tuple(k_list)
              this_tag = tag

              if (k != ''):  # If key is non-empty insert it into dictionary
                if (self.__contains__(key)):
                  test_item = self.__getitem__(key)                   
                  test_val = test_item[0]  # Value without tag
                  test_tag = test_item[1]

                  if (val != test_val):
                    print 'warning:Key "%s" already in dictionary' % \
                          (str(val)) + ' with different value (old value ' + \
                          'will be over written)'
                    
                  if (test_tag.find(this_tag) < 0):  # This tag is new
                    this_tag = test_tag+'/'+this_tag

                this_val = (val, this_tag)
                self.__setitem__(key,this_val)  # Insert key itself

              v = ll[1].lower() # Get values in a string

            elif (len(ll) == 1):  # Line contains only values - - - - - - - - -
              v = ll[0].lower() # Get values in a string

            else:
              print 'error:Illegal file format in file: "%s" in line: %s' % \
                    (fn, l)
              raise Exception

            vv = v.split(',')  # Split values into a list

            for v in vv:  # Loop over all values  - - - - - - - - - - - - - - -
              vs = v.strip()
              if (vs != ''):  # Only append non-empty values
                k_list = vs.split(' ')  # Make a list of key words
                if (len(k_list) >  self.max_key_length):
                  self.max_key_length = len(k_list) # Update maximal key length

                key = tuple(k_list)
                this_tag = tag

                if (self.__contains__(key)):
                  test_item = self.__getitem__(key)                   
                  test_val = test_item[0]  # Value without tag
                  test_tag = test_item[1]

                  if (val != test_val):
                    print 'warning:Key "%s" already in dictionary' % \
                          (str(val)) + ' with different value (old value ' + \
                          'will be over written)'

                  if (test_tag.find(this_tag) < 0):  # This tag is new
                    this_tag = test_tag+'/'+this_tag

                this_val = (val, this_tag)
                self.__setitem__(key,this_val)

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Loaded tag look-up table "%s"' % (str(self.name))
    print '1:  From files:         %s' % (str(self.file_names))
    print '1:  Number of entries:  %i' % (self.length)
    print '1:  Maximal key length: %i' % (self.max_key_length)

# =============================================================================

class FrequencyLookupTable(LookupTable):
  """class FrequencyLookupTable - Look-up tables for words and frequencies.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
 
    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.sum = None  # The sum of all frequency counts
    self.default = 1

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised frequency look-up table "%s"' % (str(self.name))
    print '1:  With default: %s' % (str(self.default))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with words and their frequency counts into the
       look-up table.
       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    if (not isinstance(file_names, list)):
      print 'error:Argument "file_names" is of wrong type, must be either '+ \
            'string or list'
      raise Exception

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table
    self.sum = 0

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      # Open file and read all lines into a list
      #
      try:
        f = open(fn,'r')
      except:
        print 'error:Can not read from file "%s"' % (fn)
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          ll = l.split(',')  # Get fields from a line

          # Check for two columns
          #
          if (len(ll) != 2):
            print 'error:Illegal file format (not 2 columns) in file: ' + \
                  '"%s" in line: %s"' % (fn, l)
            raise Exception

          key = ll[0].strip().lower()  # Make sure it's lower case
          val = ll[1].strip().lower()

          try:
            val = int(val)  # Convert the value into an integer
          except:
            print 'error:Illegal value for frequency count: "%s"' % \
                  (str(val)) + ' in line: "%s"' % (l)
            raise Exception

          if (self.__contains__(key)):
            val += self.__getitem__(key)  # Sum up counts

          self.__setitem__(key, val)
          self.sum += val

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Loaded frequency look-up table "%s"' % (str(self.name))
    print '1:  From files:        %s' % (str(self.file_names))
    print '1:  Number of entries: %i' % (self.length)
    print '1:  Sum of all value:  %i' % (self.sum)

# =============================================================================

class GeocodeLookupTable(LookupTable):
  """class GeocodeLookupTable - Look-up tables for words and their location.

     For each entry (word) in the look-up table, its longitude and latitude
     are given. The file format is three column comma separated text file
     (CSV).

     For each entry, the key is the name of the locality, and the value is a
     list with the two entries [longitude,latitude].
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
 
    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.default = []

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised geocode look-up table "%s"' % (str(self.name))
    print '1:  With default: %s' % (str(self.default))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with words and their localities into the look-up
       table.
       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    if (not isinstance(file_names, list)):
      print 'error:Argument "file_names" is of wrong type, must be either '+ \
            'string or list'
      raise Exception

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      # Open file and read all lines into a list
      #
      try:
        f = open(fn,'r')
      except:
        print 'error:Can not read from file "%s"' % (fn)
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          ll = l.split(',')  # Get fields from a line

          # Check for three columns
          #
          if (len(ll) != 3):
            print 'error:Illegal file format (not 3 columns) in file: ' + \
                  '"%s" in line: %s' % (fn, l)
            raise Exception

          key = ll[0].strip().lower()  # Make sure it's lower case
          lon = ll[1].strip()
          lat = ll[2].strip()

          try:
            lon = float(lon)  # Convert longitute into a floating-point number
          except:
            print 'error:Longitude: "%s" is not a number in line: "%s"' % \
                  (str(lon), l)
            raise Exception
          try:
            lat = float(lat)  # Convert latitude into a floating-point number
          except:
            print 'error:Lattitude: "%s" is not a number in line: "%s"' % \
                  (str(lat), l)
            raise Exception

          if (lon < -180.0) or (lon > 180.0):
            print 'error:Illegal value for longitude: '+str(lon)
            raise Exception
          if (lat < -90.0) or (lat > 90.0):
            print 'error:Illegal value for latitude: '+str(lat)
            raise Exception

          val = [lon,lat]  # Value for dictionary

          if (self.__contains__(key)) and (self.__getitem__(key) != val):
            print 'error:Key "%s" already in look-up table with ' % \
                  (str(key)) + 'different value'
            raise Exception

          self.__setitem__(key, val)

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Loaded geocode look-up table "%s"' % (str(self.name))
    print '1:  From files:        %s' % (str(self.file_names))
    print '1:  Number of entries: %i' % (self.length)

# =============================================================================

class CorrectionList(list):
  """class CorrectionList - Based on type list
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor, set general attributes.
    """

    list.__init__(self)  # Initialise list base type
    self.name =         ''
    self.description =  ''
    self.created =      ''
    self.modified =     ''
    self.file_name =    ''
    self.length =       None  # Number of entries in the correction list

    for (keyword, value) in kwargs.items():
      if (keyword == 'name'):
        self.name = value
      elif (keyword == 'description'):
        self.description = value
      elif (keyword == 'created'):
        self.created = value
      elif (keyword == 'modified'):
        self.modified = value

      else:
        print 'error:Illegal constructor argument keyword: "%s"' % \
              (str(keyword))
        raise Exception

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Initialised correction list "%s"' % (str(self.name))

  # ---------------------------------------------------------------------------

  def load(self, file_name):
    """Load one correction list file into a sorted (decreasing length) list.
       See Febrl manual for details on the file format.
    """

    # Check input argument type and open file - - - - - - - - - - - - - - - - -
    #
    if (not isinstance(file_name, str)):
      print 'error:Argument "file_name" is of wrong type, must be a string'
      raise Exception

    self.file_name = file_name

    # Make sure the list is empty, remove all items from the correction list
    #
    while (self.__len__() > 0):
      self.pop()

    # Open file and read all lines into a list
    #
    try:
      f = open(self.file_name,'r')
    except:
      print 'error:Can not read from file "%s"' % (str(self.file_name))
      raise IOError

    file_data = f.readlines()  # Read complete file
    f.close()

    org_list  = []  # List of original strings (the ones to be replaced)
    repl_list = []  # List of replacement strings
    len_list  = []  # List of original string lengths
    repl = ''       # Set inital replacement to nothing

    # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for line in file_data:
      l = line.strip()  # Remove line separators at the end

      if (len(l) > 0) and (l[0] != '#'):  # Not an empty line and not comment
        ll = l.split(':=')  # Separate replacement from values

        if (len(ll) == 2):  # Line contains a replacement - - - - - - - - - - -

          repl = ll[0].strip().lower()  # Make replacement lower and strip

          if (not ((repl[0] == '"') and (repl[-1] == '"') or \
                   (repl[0] == "'") and (repl[-1] == "'"))):
            print 'error:Replacement string is not properly quoted: '+ \
                  '"%s" in file: "%s"' % (repl, str(self.file_name))
            raise Exception

          repl = repl[1:-1]  # Remove quotes from replacement string

          v = ll[1].lower() # Get values in a string and make lowercase

        elif (len(ll) == 1):  # Line contains only values - - - - - - - - - - -
          v = ll[0].lower() # Get values in a string and make lowercase

        else:  # More than one ':=' separator in the line - - - - - - - - - - -
          print 'error:Too many ":=" separators in line: "%s"' % (l)
          raise Exception

        # Now process the values and append them to the list  - - - - - - - - -

        vv = v.split(',')  # Split values into a list

        for v in vv:  # Loop over all values  - - - - - - - - - - - - - - - -
          org = v.strip()  # Get the original string

          if (org != ''):  # Only process non-empty values
            if (not ((org[0] == '"') and (org[-1] == '"') or \
                     (org[0] == "'") and (org[-1] == "'"))):
              print 'error:Original string is not properly quoted: '+ \
                    '"%s" in file: "%s"' % (org, str(self.file_name))
              raise Exception

            org = org[1:-1]  # Remove quotes from original string

            if (org != ''):  # Only append non-empty values
              org_list.append(org)
              repl_list.append(repl)
              len_list.append(len(org))

    tmp_list = map(None,len_list,org_list,repl_list)
    tmp_list.sort()
    tmp_list.reverse()

    for (i,org,repl) in tmp_list:
      self.append((org,repl))

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message for low volume log output (level 1) - - - - - - - - - - - -
    #
    print '1:'
    print '1:Loaded correction list "%s"' % (str(self.name))
    print '1:  From file:         %s' % (str(self.file_name))
    print '1:  Number of entries: %i' % (self.length)

# =============================================================================
