# =============================================================================
# encodeTest.py - Test module for encode.py
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
# The Original Software is "encodeTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module encodeTest.py - Test module for encode.py.
"""

# -----------------------------------------------------------------------------

import unittest
import encode

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.strings = ['peter','christen','ole','nielsen','markus','hegland',
                    'stephen','steve','roberts','tim','churches','xiong',
                    'ng','miller','millar','foccachio','van de hooch',
                    'xiao ching','asawakun','prapasri','von der felde','vest',
                    'west','oioi','ohio','oihcca', 'nielsen', 'kim', 'lim',
                    'computer','record','linkage','probabilistic',
                    'aa','aaaa aaa','  x   ','x']

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testSoundex(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Soundex' string encoding"""

    c = encode.soundex('')  # Test with empty string
    assert (c == '0000'), \
           '"Soundex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.soundex(s)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return a string: '+ \
             str(type(code))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+s+ \
             '" differs from original string: '+str(c)

      assert (len(c) == 4), \
             'Length of "Soundex" code for string "'+s+'" is not four '+ \
             'characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "Soundex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

      c = encode.soundex(s,maxlen=1)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return a string: '+ \
             str(type(code))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+ \
             s+'" differs from original string: '+str(c)

      assert (len(c) == 1), \
             'Length of "Soundex" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.soundex(s,maxlen=6)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return'+ \
             ' a string: '+str(type(code))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+ \
             s+'" differs from original string: '+str(c)

      assert (len(c) <= 6), \
             '"Soundex" code for string "'+s+'" is longer than'+ \
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "Soundex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)


  def testModSoundex(self):   # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'ModSoundex' string encoding"""

    c = encode.mod_soundex('')  # Test with empty string
    assert (c == '0000'), '"ModSoundex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.mod_soundex(s)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(code))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) == 4), 'Length of "ModSoundex" code for string "'+s+'" '+\
             'is not four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "ModSoundex'+ \
               '" code for string "'+s+'" are not digits: '+str(c)

      c = encode.mod_soundex(s,maxlen=1)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(code))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) == 1), 'Length of "ModSoundex" code for string "'+s+'" '+\
             'is not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.mod_soundex(s,maxlen=6)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(code))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) <= 6), '"ModSoundex" code for string "'+s+'" is longer '+\
             'than six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "ModSoundex'+ \
               '" code for string "'+s+'" are not digits: '+str(c)


  def testPhonex(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Phonex' string encoding"""

    c = encode.phonex('')  # Test with empty string
    assert (c == '0000'), '"Phonex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.phonex(s)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) == 4), 'Length of "Phonex" code for string "'+s+'" is '+ \
             'not four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

      c = encode.phonex(s,maxlen=1)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) == 1), 'Length of "Phonex" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.phonex(s,maxlen=6)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) <= 6), '"Phonex" code for string "'+s+'" is longer than'+\
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)


  def testNYSIIS(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'NYSIIS' string encoding"""

    c = encode.nysiis('')  # Test with empty string
    assert (c == ''), '"NYSIIS" of empty string is not ""'

    for s in self.strings:

      c = encode.nysiis(s)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) <= 4), 'Length of "NYSIIS" code for string "'+s+'" is '+ \
             'more than four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in "NYSIIS" '+ \
               'code for string "'+s+'" are not letters: '+str(c)

      c = encode.nysiis(s,maxlen=1)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) == 1), 'Length of "NYSIIS" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.nysiis(s,maxlen=6)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(code))

      assert (len(c) <= 6), '"NYSIIS" code for string "'+s+'" is longer than'+\
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in "NYSIIS" '+ \
               'code for string "'+s+'" are not letters: '+str(c)


  def testDoubleMetaphone(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test 'DoubleMetaphone' string encoding"""

    c = encode.dmetaphone('')  # Test with empty string
    assert (c == ''), '"DoubleMetaphone" of empty string is not ""'

    for s in self.strings:

      c = encode.dmetaphone(s)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(code))

      assert (len(c) <= 4), 'Length of "DoubleMetaphone" code for string "'+s+\
             '" is more than four characters: '+str(c)+' with length: '+ \
             str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in '+ \
               '"DoubleMetaphone" code for string "'+s+'" are not letters: '+\
               str(c)

      c = encode.dmetaphone(s,maxlen=1)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(code))

      assert (len(c) == 1), 'Length of "DoubleMetaphone" code for string "'+s+\
             '" is '+'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.dmetaphone(s,maxlen=6)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(code))

      assert (len(c) <= 6), '"DoubleMetaphone" code for string "'+s+'" is '+ \
             'longer than six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in '+ \
               '"DoubleMetaphone" code for string "'+s+'" are not letters: '+ \
               str(c)

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
