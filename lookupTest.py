# =============================================================================
# lookupTest.py - Test module for lookup.py
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
# The Original Software is "lookupTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module lookupTest.py - Test module for lookup.py.

   Define file names for various look-up table files in 'setUp'.
"""

# -----------------------------------------------------------------------------

import unittest
import lookup

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.tag_lookup_files = ['./data/post_address.tbl',
                             './data/postcode_nsw.tbl',
                             './data/address_misc.tbl',
                             './data/saints.tbl',
                             './data/givenname_f.tbl',
                             './data/surname.tbl']

    self.frequency_lookup_files = ['./data/givenname_f_freq.csv',
                                   './data/surname_nsw_freq.csv',
                                   './data/postcode_act_freq.csv']

    self.geocode_lookup_files = ['./data/postcode_nsw_geocode.csv']

    self.correction_list_files = ['./data/address_corr.lst',
                                  './data/name_corr.lst']

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testTagLookupTables(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test tag look-up tables"""

    # First load all files separately
    #
    for f in self.tag_lookup_files:
      lookup_table = lookup.TagLookupTable(name=f, default = '')

      assert (lookup_table.name == f), \
             'Look-up table "'+f+'" does not have correct name: "'+ \
             lookup_table.name+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.max_key_length > 0), \
             'Look-up table "'+f+'" has illegal key length: '+ \
             str(lookup_table.max_key_length)

      for (key,value) in lookup_table.items():
        assert (isinstance(key,tuple)), \
               'Key in look-up table "'+f+'" is not a tuple: '+str(key)
        assert (len(value) == 2), \
               'Value in look-up table "'+f+'" does not contain two '+ \
               'elements: '+str(values)

    # Now load all files into one look-up table
    #
    lookup_table = lookup.TagLookupTable(name=self.tag_lookup_files[0], \
                   default = '')

    assert (lookup_table.name == self.tag_lookup_files[0]), \
           'Combined look-up table does not have correct name: "'+ \
           str(self.tag_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    lookup_table.load(self.tag_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.max_key_length > 0), \
           'Combined look-up table has illegal key length: '+ \
           str(lookup_table.max_key_length)

    for (key,value) in lookup_table.items():
      assert (isinstance(key,tuple)), \
             'Key in combined look-up table is not a tuple: '+str(key)
      assert (len(value) == 2), \
             'Value in combined look-up table does not contain two '+ \
             'elements: '+str(values)

  def testFrequencyLookupTables(self):  # - - - - - - - - - - - - - - - - - - -
    """Test frequency look-up tables"""

    # First load all files separately
    #
    for f in self.frequency_lookup_files:
      lookup_table = lookup.FrequencyLookupTable(name=f, default = '')

      assert (lookup_table.name == f), \
             'Look-up table "'+f+'" does not have correct name: "'+ \
             lookup_table.name+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      assert (lookup_table.sum == None), \
             'Look-up table "'+f+'" has non-None sum'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table.sum > 0), \
             'Look-up table "'+f+'" has zero sum'

      assert (lookup_table.sum > lookup_table.length), \
             'Look-up table "'+f+'" has sum smaller than number of entries: '+\
             str(self.length)

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      for (key,value) in lookup_table.items():
        assert (isinstance(key,str)), \
               'Key in look-up table "'+f+'" is not a string: '+str(key)
        assert (isinstance(value,int) and (value > 0)), \
               'Value in look-up table "'+f+'" is not an integer: '+str(value)

    # Now load all files into one look-up table
    #
    lookup_table = lookup.FrequencyLookupTable( \
                   name=self.frequency_lookup_files[0], \
                   default = '')

    assert (lookup_table.name == self.frequency_lookup_files[0]), \
           'Combined look-up table does not have correct name: "'+ \
           str(self.frequency_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    assert (lookup_table.sum == None), \
           'Combined look-up table has non-None sum'

    lookup_table.load(self.frequency_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table.sum > 0), \
           'Combined look-up table has zero sum'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    for (key,value) in lookup_table.items():
      assert (isinstance(key,str)), \
             'Key in combined look-up table is not a string: '+str(key)
      assert (isinstance(value,int) and (value > 0)), \
             'Value in combined look-up table is not an integer: '+str(value)


  def testGeocodeLookupTables(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test geocode look-up tables"""

    # First load all files separately
    #
    for f in self.geocode_lookup_files:
      lookup_table = lookup.GeocodeLookupTable(name=f, default = '')

      assert (lookup_table.name == f), \
             'Look-up table "'+f+'" does not have correct name: "'+ \
             lookup_table.name+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      for (key,value) in lookup_table.items():
        assert (isinstance(key,str)), \
               'Key in look-up table "'+f+'" is not a string: '+str(key)
        assert (len(value) == 2) and \
                isinstance(value[0],float) and isinstance(value[1],float), \
               'Value in look-up table "'+f+'" is not an location: '+str(value)

        assert (value[0] >= -180.0) and (value[0] <= 180.0), \
               'Location in look-up table "'+f+'" has illegal longitude: '+ \
               str(value[0])
        assert (value[1] >= -90.0) and (value[1] <= 90.0), \
               'Location in look-up table "'+f+'" has illegal latitude: '+ \
               str(value[1])

    # Now load all files into one look-up table
    #
    lookup_table = lookup.GeocodeLookupTable( \
                   name=self.geocode_lookup_files[0], \
                   default = '')

    assert (lookup_table.name == self.geocode_lookup_files[0]), \
           'Combined look-up table does not have correct name: "'+ \
           str(self.geocode_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    lookup_table.load(self.geocode_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    for (key,value) in lookup_table.items():
      assert (isinstance(key,str)), \
             'Key in combined look-up table is not a string: '+str(key)
      assert (len(value) == 2) and \
              isinstance(value[0],float) and isinstance(value[1],float), \
             'Value in combined look-up table is not an location: '+str(value)

      assert (value[0] >= -180.0) and (value[0] <= 180.0), \
             'Location in combined look-up table has illegal longitude: '+ \
             str(value[0])
      assert (value[1] >= -90.0) and (value[1] <= 90.0), \
             'Location in combined look-up table has illegal latitude: '+ \
             str(value[1])

  def testCorrectionLists(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test correction lists"""

    for f in self.correction_list_files:
      corr_list = lookup.CorrectionList(name=f)

      assert (corr_list.name == f), \
             'Correction list "'+f+'" does not have correct name: "'+ \
             corr_list.name+'"'

      assert (corr_list.file_name == ''), \
             'Correction list "'+f+'" has non-empty file name: '+ \
             str(lookup_table.file_name)

      assert (isinstance(corr_list,list)), \
             'Correction list "'+f+'" is not a list'

      corr_list.load(f)  # Load the table

      assert (len(corr_list) > 0), \
             'Correction list "'+f+'" is empty after load'

      elem_len = 9999  # Long elements first
      old_key = ''

      for (key,value) in corr_list:
        assert (isinstance(key,str)), \
               'Key in correction list "'+f+'" is not a string: '+str(key)
        assert (isinstance(value,str)), \
               'Value in correction list "'+f+'" is not a string: '+str(value)

        assert (len(key) <= elem_len), \
               'Correction list "'+f+'" element key "'+value+ \
               '" is longer than previous key: "'+old_key+'"'

        old_key = key
        elem_len = len(key)

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
