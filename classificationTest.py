# =============================================================================
# classificationTest.py - Test module for classification.py
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
# The Original Software is "classificationTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module classificationTest.py - Test module for classification.py.
"""

# -----------------------------------------------------------------------------

import unittest
import classification

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    class Dataset:
      def __init__(self, name):
        self.name = name

    # For data sets only define names (enough to test classifiers)
    #
    self.dataset_a = Dataset('my_data_a')
    self.dataset_b = Dataset('my_data_b')

    name_a = self.dataset_a.name
    name_b = self.dataset_b.name

    # Define various weight vectors
    #
    self.weight_vectors = [
      [name_a, 0, name_b, 0, 0.0, 0.0, 0.0, 0.0],
      [name_a, 0, name_b, 1, 1.0, 2.0, 3.0, 4.0],
      [name_a, 2, name_b, 2, 1.0, 2.0, 3.0, 4.0],
      [name_a, 0, name_b, 3, 1.0, 2.0, 3.0, 4.0, 5.0],
      [name_a, 4, name_b, 3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      [name_a, 0, name_b, 5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
      [name_a, 4, name_b, 3, -99.99, 99.99, 99.99, -99.99],
      [name_a, 0, name_b, 0, 1.0, 2.0, 3.0, 4.0]
    ]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testFellegiSunterClassifier(self):  # - - - - - - - - - - - - - - - - - -
    """Test Fellegi and Sunter classifier"""

    my_class = classification.FellegiSunterClassifier(name = 'fell/sunter',
                                                 dataset_a = self.dataset_a,
                                                 dataset_b = self.dataset_b,
                                           lower_threshold = 0.0,
                                           upper_threshold = 5.0)

    num_vectors = len(self.weight_vectors)

    weight_results = [0.0, 10.0, 10.0, 15.0, 0.0, 9.0, 0.0, 10.0]
 
    for i in range(num_vectors):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      vec = self.weight_vectors[i]
      res = weight_results[i]

      my_class.classify(vec)

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from Fellegi and Sunter classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)

    # Now classify a block of vectors
    #
    my_class.classify_block(self.weight_vectors)

    weight_results = [10.0, 10.0, 10.0, 15.0, 0.0, 9.0]

    for i in range(len(weight_results)):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      res = weight_results[i]

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from Fellegi and Sunter classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)


  def testFlexibleClassifier(self):   # - - - - - - - - - - - - - - - - - - - -
    """Test Flexible classifier"""

    my_class = classification.FlexibleClassifier(name = 'flexible',
                                            dataset_a = self.dataset_a,
                                            dataset_b = self.dataset_b,
                                      lower_threshold = 0.0,
                                      upper_threshold = 5.0,
                                            calculate = [('max', [0,1]),
                                                         ('min', [1,2]),
                                                         ('add', [2,3])],
                                          final_funct = 'add')

    num_vectors = len(self.weight_vectors)

    weight_results = [0.0, 11.0, 11.0, 11.0, 0.0, 4.0, 199.98, 11.0]
 
    for i in range(num_vectors):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      vec = self.weight_vectors[i]
      res = weight_results[i]

      my_class.classify(vec)

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from flexible classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)

    # Now classify a block of vectors
    #
    my_class.classify_block(self.weight_vectors)

    weight_results = [11.0, 11.0, 11.0, 11.0, 199.98, 4.0]

    for i in range(len(weight_results)):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      res = weight_results[i]

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from flexible classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)

    # Now define other calculation functions
    #
    del my_class

    my_class = classification.FlexibleClassifier(name = 'flexible',
                                            dataset_a = self.dataset_a,
                                            dataset_b = self.dataset_b,
                                      lower_threshold = 0.0,
                                      upper_threshold = 5.0,
                                            calculate = [('add', [0,1,2]),
                                                         ('min', [1,2,3]),
                                                         ('mult', [0,3])],
                                          final_funct = 'avrg')

    num_vectors = len(self.weight_vectors)

    weight_results = [0.0, 4.0, 4.0, 4.0, 0.0, 1.6666666666666667, \
                      3332.6666999999998, 4.0]
 
    for i in range(num_vectors):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      vec = self.weight_vectors[i]
      res = weight_results[i]

      my_class.classify(vec)

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from flexible classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)

    # Now classify a block of vectors
    #
    my_class.classify_block(self.weight_vectors)

    weight_results = [4.0, 4.0, 4.0, 4.0, 3332.6666999999998, \
                      1.6666666666666667]

    for i in range(len(weight_results)):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      res = weight_results[i]

      class_res_dict = my_class.results[rec_id_a]
      class_res = class_res_dict[rec_id_b]

      assert (class_res == res), \
             'Wrong result returned from flexible classifier for' + \
             ' vector "%s": %f (should be: %f)' % (str(vec), class_res, res)

  def testMergingClassifiers(self):   # - - - - - - - - - - - - - - - - - - - -
    """Test the merging of classifiers"""

    my_class1 = classification.FellegiSunterClassifier(name = 'fell/sunter',
                                                  dataset_a = self.dataset_a,
                                                  dataset_b = self.dataset_b,
                                            lower_threshold = 0.0,
                                            upper_threshold = 5.0)

    my_class2 = classification.FlexibleClassifier(name = 'flexible',
                                             dataset_a = self.dataset_a,
                                             dataset_b = self.dataset_b,
                                       lower_threshold = 0.0,
                                       upper_threshold = 5.0,
                                             calculate = [('max', [0,1]),
                                                          ('min', [1,2]),
                                                          ('add', [2,3])],
                                           final_funct = 'add')

    my_class3 = classification.FlexibleClassifier(name = 'flexible',
                                             dataset_a = self.dataset_a,
                                             dataset_b = self.dataset_b,
                                       lower_threshold = 0.0,
                                       upper_threshold = 5.0,
                                             calculate = [('add', [0,1,2]),
                                                          ('min', [1,2,3]),
                                                          ('mult', [0,3])],
                                           final_funct = 'avrg')

    num_vectors = len(self.weight_vectors)

    for i in range(num_vectors):
      rec_id_a = self.weight_vectors[i][1]
      rec_id_b = self.weight_vectors[i][3]

      vec = self.weight_vectors[i]

      my_class1.classify(vec)
      my_class2.classify(vec)
      my_class3.classify(vec)

    my_class1.merge(my_class2.results)
    my_class1.merge(my_class3.results)
    my_class2.merge(my_class3.results)
    my_class3.merge(my_class2.results)

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
