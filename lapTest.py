# =============================================================================
# lapTest.py - Test module for lap.py
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
# The Original Software is "lapTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module lapTest.py - Test module for lap.py.
"""

TEST_ALGO = 'auction'
# TEST_ALGO = 'lapmod'

# -----------------------------------------------------------------------------

import random
import unittest

import lap

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):  # Define results dictionaries and expected assignments

    self.res_dict = []  # Results dictionaries
    self.ass_thre = []  # Thresholds for assignment filter
    self.pro_type = []  # The process type
    self.ass_resu = []  # The expected results


    # Simple results dictionary with one sub-set (tests 0 and 1)
    #
    self.res_dict.append({1:{2:99.99}, 2:{3:55.55}, 3:{4:111.111}})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(1,2):99.99, (3,4):111.111})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(1,2):99.99, (2,3):55.55, (3,4):111.111})

 
    # Same, but different threshold (filter out one pair) (tests 2 and 3)
    #
    self.res_dict.append({1:{2:99.99}, 2:{3:55.55}, 3:{4:111.111}})
    self.ass_thre.append(60.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(1,2):99.99, (3,4):111.111})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(1,2):99.99, (3,4):111.111})


    # Same but different threshold (filter out two pairs) (tests 4 and 5)
    #
    self.res_dict.append({1:{2:99.99}, 2:{3:55.55}, 3:{4:111.111}})
    self.ass_thre.append(100.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(3,4):111.111})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(3,4):111.111})


    # Another small results disctionary with one sub-set (tests 6 and 7)
    #
    self.res_dict.append({  4:{ 66:175.42, 102:176.88, 120:177.93}, \
                           66:{102:176.63, 120:176.95}, \
                          102:{120:175.49}})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(4,120):177.93, (66,102):176.63})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(4,66):175.42, (66,102):176.63, (102,120):175.49})


    # Now a results dictionary with three indepdendent sub-sets (tests 8 and 9)
    #
    # Sub-set 1: rows: 4, 66, 102
    # Sub-set 2: rows 235, 305
    # Sub-set 3: rows: 269, 295 (one column only)
    #
    self.res_dict.append({  4:{ 66:175.42, 102:176.88, 120:177.93}, \
                           66:{102:176.63, 120:176.95}, \
                          102:{120:175.49}, \
                          235:{305:17.44, 364:173.01, 400:170.88}, \
                          269:{429:169.09}, \
                          295:{429:175.99}, \
                          305:{400:176.78}})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(4,120):177.93, (66,102):176.63, (295,429):175.99, \
                          (305,400):176.78, (235,364):173.01})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(4,66):175.42, (66,102):176.63, (102,120):175.49, \
                          (295,429):175.99, (305,400):176.78, \
                          (235,364):173.01})


    # Now a results dictionary with one row only (tests 10 and 11)
    #
    self.res_dict.append({99:{ 33:175.42, 66:176.88, 333:177.93, 11:175.42}})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(99,333):177.93})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(99,333):177.93})

    # Now a results dictionary with one element per row only (tests 12 and 13)
    #
    self.res_dict.append({33:{333:184.42},
                          66:{333:191.88},
                          99:{333:199.99},
                          11:{333:171.42}})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(99,333):199.99})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(99,333):199.99})


    # Now a bit more complex results dictionary from a real world data set
    # (tests 14 and 15)
    #
    self.res_dict.append({  241: {6545:32.9562},
                           1076: {6545:31.3719, 1671:33.5029, 4943:43.7987}, \
                           1267: {6545:34.5932, 6914:45.6304, 6148:33.9058}, \
                           1569: {6545:33.6152, 1671:48.1276}, \
                           2217: {6545:32.1244}, \
                           2388: {6545:34.1836}, \
                           4943: {14405:30.3479}, \
                           6148: {14617:31.1298}, \

})
    self.ass_thre.append(20.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(1569,1671):48.1276, (1267,6914):45.6304, \
                          (1076,4943):43.7987, (2388,6545):34.1836, \
                          (6148,14617): 31.1298})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(1569,1671):48.1276, (1267,6914):45.6304, \
                          (1076,4943):43.7987, (2388,6545):34.1836, \
                          (4943,14405):30.3479, (6148,14617): 31.1298})


    # And another real world results dictionary that first resulted in an
    # endless loop (tests 16 and 17)
    #
    self.res_dict.append({259: {999:27.1033, 828:25.8376, 783:20.4161}, \
                          828: {878:24.7124, 999:43.3455}, \
                          878: {999:26.0539}, \
                          783: {828:25.0573, 999:23.7814}})
    self.ass_thre.append(10.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(828,999):43.3455, (259,783):20.4161})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(878,999):26.0539, (828,878):24.7124, \
                          (783,828):25.0573, (259,783):20.4161})


    # And another real world results dictionary that first resulted in an
    # endless loop (tests 18 and 19)
    #
    self.res_dict.append({  8: {712:98},
                          670: {712:154, 938:161, 711:70},
                          711: {712:231},
                          536: {712:91, 670:140}, \
                           94: {938:54},
                          110: {938:105}})
    self.ass_thre.append(10.0)
    self.pro_type.append('deduplication')
    self.ass_resu.append({(536,670):140, (110,938):105, (711,712):231})

    self.res_dict.append(self.res_dict[-1])  # Same as above
    self.ass_thre.append(self.ass_thre[-1])  # Same as above
    self.pro_type.append('linkage')
    self.ass_resu.append({(536,670):140, (670,711):70, (711,712):231, \
                          (110,938):105})

    # And another real world results dictionary that first resulted in an
    # endless loop (tests 20 and 21)
    #
    self.res_dict.append({100549: {101172: 54.591968564734351,
                                   100756: 54.591968564734351,
                                   100549: -257.51421074870666},
                          100756: {101172: 257.51421074870666,
                                   100756: -257.51421074870666,
                                   100549: 54.591968564734351},
                          101172: {100756: 257.51421074870666,
                                   101172: -257.51421074870666,
                                   100549: 54.591968564734351}})

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def test_simple(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test lap.do_lap routine with simple dictionaries"""

    # Loop over all defined tests
    #
    do_tests = [0,2,4,6,8,10,12,14,16,18]  #  Deduplication process
    do_tests += [1,3,5,7,9,11,13,15,17,19]  # Linkage process

    for i in do_tests:
      test_lap_results = lap.do_lap(TEST_ALGO, self.res_dict[i], \
                                    self.pro_type[i], self.ass_thre[i])

      # Check results for correctness
      #
      assert (len(test_lap_results) == len(self.ass_resu[i])), \
             '"%s" returned LAP result of wrong length: %i, should be: %i' % \
             (TEST_ALGO, len(test_lap_results), len(self.ass_resu[i])) + \
             ', results returned: %s' % (str(test_lap_results))

      lap_keys = self.ass_resu[i].keys()
      lap_keys.sort()

      for k in lap_keys:
        assert (k in test_lap_results), \
               'Record pair "%s" not in LAP results: %s' % \
               (str(k), str(test_lap_results))

        assert (True == test_lap_results[k]), \
               'Record pair "%s" with wrong weight: %f, should be: %f' % \
               (str(k), test_lap_results[k], self.ass_resu[i][k])

  def test_random_simple(self):   # - - - - - - - - - - - - - - - - - - - - - -
    """Test lap.do_lap routine with random generated results dictionaries"""

    random.seed()

    dim = 1000
    thresh = 20.0
    spread = 100.0

    # First create a results dictionary with one element per row
    #
    col_nums = range(dim)
    random.shuffle(col_nums)  # Shuffle randomly

    lower = thresh+1.0
    upper = lower + spread
    res_dict = {}
    for i in range(dim):
      res_dict[i] = {col_nums[i]:random.uniform(lower, upper)}

    print res_dict

    test_lap_results = lap.do_lap(TEST_ALGO,res_dict,'deduplication',thresh)


  def test_larger1(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test lap.do_lap routine with larger dictionaries"""

    self.results_dict = { 1:{1:96.04, 6:87.13, 10:49.51}, \
                          2:{5:83.17, 6:66.34,  9:54.46}, \
                          3:{4:60.40, 5:96.04,  8:37.63}, \
                          4:{3:46.54, 5:59.41,  9:28.72}, \
                          5:{2:53.47, 5: 1.99,  6:63.37}, \
                          6:{3:91.09, 5:46.54,  6:98.02}, \
                          7:{2:46.54, 5:60.40,  8:54.46}, \
                          8:{4:75.25, 8:93.07, 10:23.77}, \
                          9:{2: 9.91, 3: 9.91,  4:80.20}, \
                         10:{5:94.06, 7:47.53,  9:97.03}}

    self.process_type = 'linkage'  # Define a process type

    self.results_pairs = {(1,1):96.04, (2,9):54.46, (3,5):96.04, (4,3):46.54, \
                          (5,2):53.47, (6,6):98.02, (7,8):54.46, \
                          (8,10):23.77, (9,4):80.20, (10,7):47.53}

    self.threshold =   10.0  # Should result in the same LAP results
    self.lap_results = {}

    test_lap_results = lap.do_lap(TEST_ALGO, self.results_dict, \
                                  self.process_type, self.threshold)

    # Check results for correctness
    #
    assert (len(test_lap_results) == 10), \
           '"%s" returned LAP result of wrong length: %i, should be: 10' % \
           (TEST_ALGO, len(test_lap_results)) + ', results returned: %s' % \
           (str(test_lap_results))

    lap_keys = self.results_pairs.keys()
    lap_keys.sort()

    for k in lap_keys:
      assert (k in test_lap_results), \
             'Record pair "%s" not in LAP results: %s' % \
             (str(k), str(test_lap_results))

      assert (True == test_lap_results[k]), \
             'Record pair "%s" with wrong weight: %f, should be: %f' % \
             (str(k), test_lap_results[k], self.results_pairs[k])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    self.process_type = 'deduplication'  # Define a process type

    self.results_pairs = {(1,1):96.04, (3,5):96.04, (6,6):98.02, \
                          (7,2):46.54, (8,4):75.25, (10,9):97.03}

    self.threshold =   10.0  # Should result in the same LAP results
    self.lap_results = {}

    test_lap_results = lap.do_lap(TEST_ALGO, self.results_dict, \
                                  self.process_type, self.threshold)

    # Check results for correctness
    #
    assert (len(test_lap_results) == 6), \
           '"%s" returned LAP result of wrong length: %i, should be: 6' % \
           (TEST_ALGO, len(test_lap_results)) + ', results returned: %s' % \
           (str(test_lap_results))

    lap_keys = self.results_pairs.keys()
    lap_keys.sort()

    for k in lap_keys:
      assert (k in test_lap_results), \
             'Record pair "%s" not in LAP results: %s' % \
             (str(k), str(test_lap_results))

      assert (True == test_lap_results[k]), \
             'Record pair "%s" with wrong weight: %f, should be: %f' % \
             (str(k), test_lap_results[k], self.results_pairs[k])

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
