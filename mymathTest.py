# =============================================================================
# mymathTest.py - Test module for mymath.py
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
# The Original Software is "mymathTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module mymathTest.py - Test module for mymath.py.
"""

# -----------------------------------------------------------------------------

import unittest
import mymath

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # Vectors plus their mean and standard deviations
    #
    self.vectors = [([1],             1.0,  0.0),
                    ([1,2.0],         1.5,  0.5),
                    ([10,100],       55.0, 45.0),
                    ([10.0,100.0],   55.0, 45.0),
                    ([1,2,3,4,5,6,7], 4.0, 2.0)]

    # Numbers and their log2 values
    #
    self.log2numbers = [(1,0),(2,1),(4,2),(8,3),(16,4),(32,5),(64,6),(1024,10)]

    # A list of tag sequences and their permuations (number and permutatations)
    #
    self.tag_lists = [(['a','b'],1,[['a','b']]),
                      ([['a','c'],'b'],2,[['a','b'],['c','b']]),
                      (['a',['b','c']],2,[['a','b'],['a','c']]),
                      (['a','b','c','d'],1,[['a','b','c','d']]),
                      ([['a','b'],['c','d']],4,[['a','c'],['b','c'], \
                                                ['a','d'],['b','d']]),
                      ([['a','b'],['c','d','e'],'f'],6,[['a','c','f'], \
                                                        ['b','c','f'], \
                                                        ['a','d','f'], \
                                                        ['b','d','f'], \
                                                        ['a','e','f'], \
                                                        ['b','e','f']]),
                     ]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testMean(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'mean' routine"""

    for i in self.vectors:

      m = mymath.mean(i[0])

      assert (isinstance(m,float)), \
             'Value returned from "mean" is not a float: '+str(m)

      assert m == i[1], \
             'Wrong "mean" with data: '+str(i[0])+' (should be: '+str(i[1])+ \
             '): '+str(m)

  def testStdDev(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'stddev' routine"""

    for i in self.vectors:

      s = mymath.stddev(i[0])

      assert (isinstance(s, float)), \
             'Value returned from "stddev" is not a float: '+str(s)

      assert s == i[2], \
             'Wrong "stddev" with data: '+str(i[0])+' (should be: '+ \
             str(i[2])+'): '+str(s)

  def testLog2(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'log2' routine"""

    for n in self.log2numbers:

      l = mymath.log2(n[0])

      assert (isinstance(l, float)), \
             'Value returned from "log2" is not a float: '+str(l)

      assert l == n[1], \
             'Wrong "log2" with value: '+str(n[0])+' (should be: '+ \
             str(n[1])+'): '+str(l)

  def testPermTagSeq(self):   # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'perm_tag_sequence' routine"""

    for l in self.tag_lists:

      t = mymath.perm_tag_sequence(l[0])

      assert len(t) == l[1], \
             '"perm_tag_sequence" returns wrong number of permutations with '+\
             'list: '+str(l[0])+' (should be: '+str(l[1])+'): '+str(len(t))

      for i in range(len(t)):
        assert t[i] == l[2][i], \
               '"perm_tag_sequence" returns wrong permutation: '+str(t[i])+ \
               ', should be: '+str(l[2][i])

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
