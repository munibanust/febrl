# =============================================================================
# stringcmpTest.py - Test module for stringcmp.py
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
# The Original Software is "stringcmpTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module stringcmpTest.py - Test module for stringcmp.py.
"""

# -----------------------------------------------------------------------------

import unittest
import stringcmp

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.string_pairs = [['peter',            'peter'            ],
                         ['peter',            ' peter'           ],
                         ['peter',            'xanthalope"'      ],
                         ['christensen',      'christensen'      ],
                         ['prapasi asawakun', 'prapasi asawakun',],
                         ['prapasi asawakun', 'paprasi asawakun',],
                         ['shackleford',      'shackelford'      ],
                         ['dunningham',       'cunnigham'        ],
                         ['nichleson',        'nichulson'        ],
                         ['jones',            'johnson'          ],
                         ['massey',           'massie'           ],
                         ['abroms',           'abrams'           ],
                         ['hardin',           'martinez'         ],
                         ['aa',               'aaa'              ],
                         ['aab',              'aaa'              ],
                         ['  aab  ',          ' a  a  a  '       ],
                         ['aaaaaaaaaa',       'aaaaaaaaab'       ],
                         ['aaaaaaaaaa',       'aaaaaaaa'         ],
                         ['    ',             '   '              ],
                         [' ',                '      '           ],
                         [' ',                ' #    '           ],
                         ['1234567890',       '2345678901'       ],
                         ['  re$ mkM )"- ',   ' re$ mkM )"- '    ],
                         ['  re$ mkM )"- ',   '  re$ mkM )"- '   ],
                         ['  re$% mkM %")- ', ' " or4%~~!][{ . ' ]]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testJaro(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Jaro' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.jaro(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Jaro" does not return a floating point number for: '+str(pair)

      assert (approx_str_value >= 0.0), \
             '"Jaro" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Jaro" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.jaro(pair[0],pair[1])
      approx_str_value_2 = stringcmp.jaro(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Jaro" returns different values for pair and swapped pair: '+ \
             str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Jaro" does not return 1.0 if strings are equal: '+str(pair)


  def testWinkler(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Winkler' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.winkler(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Winkler" does not return a floating point number for:'+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"Winkler" returns a negative number for:'+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Winkler" returns a number larger than 1.0 for:'+str(pair)

      approx_str_value_1 = stringcmp.winkler(pair[0],pair[1])
      approx_str_value_2 = stringcmp.winkler(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Winkler" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Winkler" does not return 1.0 if strings are equal: '+str(pair)

      # Winkler should always return a value equal to or larger than Jaro
      #
      approx_str_value_winkler = stringcmp.winkler(pair[0],pair[1])
      approx_str_value_jaro =    stringcmp.jaro(pair[0],pair[1])

      assert (approx_str_value_winkler >= approx_str_value_jaro), \
             '"Winkler" value smaller than "Jaro" value for:'+str(pair)


  def testBigram(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Bigram' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.bigram(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Bigram" does not return a floating point number for: '+str(pair)

      assert (approx_str_value >= 0.0), \
             '"Bigram" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Bigram" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.bigram(pair[0],pair[1])
      approx_str_value_2 = stringcmp.bigram(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Bigram" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Bigram" does not return 1.0 if strings are equal: '+str(pair)


  def testEditDist(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'EditDist' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.editdist(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"EditDist" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"EditDist" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"EditDist" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.editdist(pair[0],pair[1])
      approx_str_value_2 = stringcmp.editdist(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"EditDist" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"EditDist" does not return 1.0 if strings are equal: '+ \
               str(pair)


  def testSeqMatch(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'SeqMatch' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.seqmatch(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"SeqMatch" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"SeqMatch" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"SeqMatch" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.seqmatch(pair[0],pair[1])
      approx_str_value_2 = stringcmp.seqmatch(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"SeqMatch" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"SeqMatch" does not return 1.0 if strings are equal: '+ \
               str(pair)


# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
