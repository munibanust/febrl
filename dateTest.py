# =============================================================================
# dateTest.py - Test module for date.py
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
# The Original Software is "dateTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module dateTest.py - Test module for date.py.
"""

# -----------------------------------------------------------------------------

import unittest
import date

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # A list with dates and corresponding date tuples and epoch numbers
    #
    self.dates = [['Sep 1, 68',        ['01','09','1968'], 25080, 14],
                  ['18 Jan 2002',      ['18','01','2002'], 37272,  1],
                  ['17:2:2002',        ['17','02','2002'], 37302,  0],
                  ['2002-02-25',       ['25','02','2002'], 37310,  4],
                  ['18,03,2001',       ['18','03','2001'], 36966,  0],
                  ['21.12.1999',       ['21','12','1999'], 36513,  0],
                  ['February 18,19',   ['18','02','1919'],  6987, 14],
                  ['23\\July\\1968',   ['23','07','1968'], 25040,  1],
                  ['18-02-2002',       ['18','02','2002'], 37303,  0],
                  ['5/03/01',          ['05','03','2001'], 36953,  9],
                  ['19680429',         ['29','04','1968'], 24955,  6],
                  ['600810',           ['10','08','1960'], 22136, 15],
                  ['3:05:2000',        ['03','05','2000'], 36647,  0],
                  ['30.11.1989',       ['30','11','1989'], 32840,  0],
                  ["1. January '70",   ['01','01','1970'], 25567, 10],
                  ['01011970',         ['01','01','1970'], 25567,  8],
                  ['10011970',         ['10','01','1970'], 25576,  7],
                  ['31 dec  1969',     ['31','12','1969'], 25566,  1],
                  ['30 december  69',  ['30','12','1969'], 25565, 10],
                  ['01011970',         ['01','01','1970'], 25567,  8],
                  ['13 Feb 1945',      ['13','02','1945'], 16479,  1],
                  ['Feb 13, \'45',     ['13','02','1945'], 16479, 14],
                  ['April 29 1968',    ['29','04','1968'], 24955,  3],
                  ['29-4=68',          ['29','04','1968'], 24955,  9],
                  ['11-01-1972',       ['11','01','1972'], 26307,  0],
                  ['January 10. 1972', ['10','01','1972'], 26306,  3],
                  ['29 Feb 1932',      ['29','02','1932'], 11746,  1],
                  ['29 Feb 32',        ['29','02','1932'], 11746, 10],
                  ['11 Jun 1902',      ['11','06','1902'],   891,  1],
                  ['11 Jul 1989',      ['11','07','1989'], 32698,  1],
                  ['12111968',         ['12','11','1968'], 25152,  7],
                  ['      21111969  ', ['21','11','1969'], 25526,  7]]

    # Define a list of date parsing format strings
    #
    self.date_parse_formats = ['%d %m %Y',   # 24 04 2002  or  24 4 2002
                               '%d %B %Y',   # 24 Apr 2002 or  24 April 2002
                               '%m %d %Y',   # 04 24 2002  or  4 24 2002
                               '%B %d %Y',   # Apr 24 2002 or  April 24 2002
                               '%Y %m %d',   # 2002 04 24  or  2002 4 24
                               '%Y %B %d',   # 2002 Apr 24 or  2002 April 24
                               '%Y%m%d',     # 20020424        ISO standard
                               '%d%m%Y',     # 24042002
                               '%m%d%Y',     # 04242002
                               '%d %m %y',   # 24 04 02    or  24 4 02
                               '%d %B %y',   # 24 Apr 02   or  24 April 02
                               '%y %m %d',   # 02 04 24    or  02 4 24
                               '%y %B %d',   # 02 Apr 24   or  02 April 24
                               '%m %d %y',   # 04 24 02    or  4 24 02
                               '%B %d %y',   # Apr 24 02   or  April 24 02
                               '%y%m%d',     # 020424
                               '%d%m%y',     # 240402
                               '%m%d%y',     # 042402

                              #'%Y %B %d',   # 2002 Apr 24  ??? Is this needed?
                              #'%y %B %d',   # 02 Apr 24   ??? Is this needed?
                              ]
    self.pivot_year = 3

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testDateToEpoch(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'date_to_epoch' string encoding"""

    for d in self.dates:
      day =   d[1][0]
      month = d[1][1]
      year =  d[1][2]

      epoch = date.date_to_epoch(day,month,year)

      assert (epoch == d[2]), \
             'Wrong epoch number for date '+str(d[1])+' (is '+str(epoch)+ \
             ' but should be '+str(d[2])+')'

  def testEpochToDate(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'epoch_to_date' string encoding"""

    for d in self.dates:
      epoch = d[2]

      d1 = date.epoch_to_date(epoch)

      assert (d1 == d[1]), \
             'Wrong date for epoch number '+str(epoch)+' (is '+str(d1)+ \
             ' but should be '+str(d[1])+')'


  def testDateToAge(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'date_to_age' string encoding"""

    for d in self.dates:
      day =   d[1][0]
      month = d[1][1]
      year =  d[1][2]

      age1 = date.date_to_age(day,month,year)
      assert (age1 > 0), 'Illegal age for date "'+d[0]+'": '+str(age1)

      age2 = date.date_to_age(day,month,year, fix_date='today')
      assert (age2 > 0), 'Illegal age for date "'+d[0]+'": '+str(age2)

      assert (age1 == age2), 'Age 1 and age 2 differ for date "'+d[0]+ \
             '": '+str(age1)+' / '+str(age2)

      assert (age1 < 105.0), \
             'Age 1 for date "'+str(d[0])+'" is not smaller than 100.0: '+ \
             str(age1)

      age3 = date.date_to_age(day,month,year, fix_date=[1,9,2002])
      assert (age3 > 0), 'Illegal age for date "'+d[0]+'" with fix date "'+ \
             str([1,9,2002])+'": '+str(age3)

      assert (age3 < age1), \
             'Age 1 for date "'+str(d[0])+'" is not larger than age 3: '+ \
             str(age1)+' / '+str(age3)

      age4 = date.date_to_age(day,month,year, fix_date= [1, 1,1902])
      assert (age4 < 0), 'Illegal age for date "'+d[0]+'" with fix date "'+ \
             str([1,1,1902])+'": '+str(age4)

    age5 = date.date_to_age(1, 1, 1903, fix_date= [1, 1,2003])
    assert (age5 == 100.0), 'Illegal age for date "[1,1,1903]" with fix '+ \
           'date "[1,1,2003]": '+str(age5)

  def testStrToDate(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'str_to_date' string encoding"""

    for d in self.dates:

      format_str = self.date_parse_formats[d[3]]

      d1 = date.str_to_date(d[0], format_str, self.pivot_year)

      assert d1 == d[1], 'Wrongly parsed date "'+d[0]+'" into: '+ str(d1)

  def testGetToday(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'get_today' string encoding"""

    d = date.get_today()

    assert (isinstance(d,list) and len(d) == 3), \
           'No list of length 3 is returned'

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
