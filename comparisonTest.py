# =============================================================================
# comparisonTest.py - Test module for comparison.py
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
# The Original Software is "comparisonTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module comparisonTest.py - Test module for comparison.py.

   TODO:
   - PC, 27/11/2002: Tests with frequency tables
"""

# -----------------------------------------------------------------------------

import unittest
import comparison
import mymath
import dataset  # Needed for setting up comparison vectors
import lookup

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.fields={'givenname':0, 'surname':1, 'postcode':2,'dob':3,'time':4,
                 'mob':5,'yob':6}

    self.dataset = dataset.DataSet({'name':'comparisonTest',
                                    'fields':self.fields,
                                    'access_mode':'read',
                                    'missing_values':['','missing']})

    # Define M- and U-probabilities and weights
    #
    self.m_prob = 0.95
    self.u_prob = 0.01

    self.agree_weight = mymath.log2(self.m_prob / self.u_prob)
    self.disagree_weight = mymath.log2((1.0-self.m_prob) / (1.0-self.u_prob))
    self.miss_weight = 0.0

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testExactString(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorExactString comparison routine"""

    exact_compare = \
      comparison.FieldComparatorExactString(fields_a='surname',
                                            fields_b='surname',
                                            m_prob=self.m_prob,
                                            u_prob=self.u_prob,
                                            missing_weight=self.miss_weight)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [exact_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                     {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar = float(similar[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (similar == self.disagree_weight), \
           'Wrong similar weight (should be: '+str(self.disagree_weight)+ \
           '): '+str(similar)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testTruncateString(self):   # - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorTruncateString comparison routine"""

    truncate_compare = \
      comparison.FieldComparatorTruncateString(fields_a='surname',
                                               fields_b='surname',
                                               m_prob=self.m_prob,
                                               u_prob=self.u_prob,
                                               missing_weight=self.miss_weight,
                                               max_string_length=3)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [truncate_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (similar1 == self.agree_weight), \
           'Wrong similar1 weight (should be: '+str(self.agree_weight)+ \
           '): '+str(similar1)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be: '+str(self.disagree_weight)+ \
           '): '+str(similar2)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)

    truncate_compare = \
      comparison.FieldComparatorTruncateString(fields_a='surname',
                                               fields_b='surname',
                                               m_prob=self.m_prob,
                                               u_prob=self.u_prob,
                                               missing_weight=self.miss_weight,
                                               max_string_length=4)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [truncate_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (similar1 == self.disagree_weight), \
           'Wrong similar1 weight (should be: '+str(self.disagree_weight)+ \
           '): '+str(similar1)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be: '+str(self.disagree_weight)+ \
           '): '+str(similar2)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testApproxString(self):   # - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorApproxString comparison routine"""

    approx_compare = \
      comparison.FieldComparatorApproxString(fields_a='surname',
                                             fields_b='surname',
                                             m_prob=self.m_prob,
                                             u_prob=self.u_prob,
                                             missing_weight=self.miss_weight,
                                             compare_method='jaro',
                                             min_approx_value=0.5)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [approx_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 > self.disagree_weight), \
           'Wrong similar2 weight (should be larger than disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


    approx_compare = \
      comparison.FieldComparatorApproxString(fields_a='surname',
                                             fields_b='surname',
                                             m_prob=self.m_prob,
                                             u_prob=self.u_prob,
                                             missing_weight=self.miss_weight,
                                             compare_method='jaro',
                                             min_approx_value=0.8)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [approx_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


    approx_compare = \
      comparison.FieldComparatorApproxString(fields_a='surname',
                                             fields_b='surname',
                                             m_prob=self.m_prob,
                                             u_prob=self.u_prob,
                                             missing_weight=self.miss_weight,
                                             compare_method='bigram',
                                             min_approx_value=0.5)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [approx_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petea','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testEncodeString(self):   # - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorEncodeString comparison routine"""

    encode_compare = \
      comparison.FieldComparatorEncodeString(fields_a='surname',
                                             fields_b='surname',
                                             m_prob=self.m_prob,
                                             u_prob=self.u_prob,
                                             missing_weight=self.miss_weight,
                                             encode_method='soundex',
                                             max_code_length=2)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [encode_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight == similar1), \
           'Wrong similar1 weight (should be equal to agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 == similar2), \
           'Wrong similar2 weight (should be equal to similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 == self.agree_weight), \
           'Wrong similar2 weight (should be equal to agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)

    encode_compare = \
      comparison.FieldComparatorEncodeString(fields_a='surname',
                                             fields_b='surname',
                                             m_prob=self.m_prob,
                                             u_prob=self.u_prob,
                                             missing_weight=self.miss_weight,
                                             encode_method='nysiis',
                                             max_code_length=5,
                                             reverse=True)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [encode_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'pedra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.disagree_weight == similar1), \
           'Wrong similar1 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar1)

    assert (similar1 == similar2), \
           'Wrong similar2 weight (should be equal to similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testKeyDiff(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorKeyDiff comparison routine"""

    keydiff_compare = \
      comparison.FieldComparatorKeyDiff(fields_a='surname',
                                        fields_b='surname',
                                        m_prob=self.m_prob,
                                        u_prob=self.u_prob,
                                        missing_weight=self.miss_weight,
                                        max_key_diff=1)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [keydiff_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petea','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)

    keydiff_compare = \
      comparison.FieldComparatorKeyDiff(fields_a='surname',
                                        fields_b='surname',
                                        m_prob=self.m_prob,
                                        u_prob=self.u_prob,
                                        missing_weight=self.miss_weight,
                                        max_key_diff=3)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [keydiff_compare])

    agree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petea','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'petra','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'surname':'xanthilla','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'surname':'peter','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'surname':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'surname':'tim','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 > self.disagree_weight), \
           'Wrong similar2 weight (should be larger than disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testNumericPerc(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorNumericPerc comparison routine"""

    numperc_compare = \
      comparison.FieldComparatorNumericPerc(fields_a='postcode',
                                            fields_b='postcode',
                                            m_prob=self.m_prob,
                                            u_prob=self.u_prob,
                                            missing_weight=self.miss_weight,
                                            max_perc_diff=10.0)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [numperc_compare])

    agree = rec_comparator.compare({'postcode':'12345','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'postcode':'12345','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'995','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'910.1','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    similar3 = rec_comparator.compare({'postcode':'1000.0','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'880','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar3 = float(similar3[4])

    disagree = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'99','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'postcode':'999','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'postcode':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'postcode':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'postcode':'678','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 > similar3), \
           'Wrong similar3 weight (should be smaller than similar2 weight '+ \
           str(similar2)+'): '+ str(similar3)

    assert (similar2 > self.disagree_weight), \
           'Wrong similar2 weight (should be larger than disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (similar3 == self.disagree_weight), \
           'Wrong similar3 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar3)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testNumericAbs(self):   # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorNumericAbs comparison routine"""

    numabs_compare = \
      comparison.FieldComparatorNumericAbs(fields_a='postcode',
                                            fields_b='postcode',
                                            m_prob=self.m_prob,
                                            u_prob=self.u_prob,
                                            missing_weight=self.miss_weight,
                                            max_abs_diff=10)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [numabs_compare])

    agree = rec_comparator.compare({'postcode':'12345.34','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'postcode':'12345.34','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'postcode':'1000.001','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'995','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'990','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    similar3 = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'980.9876','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    similar3 = float(similar3[4])

    disagree = rec_comparator.compare({'postcode':'1000','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                      {'postcode':'99','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'postcode':'999.098871','_rec_num_':0,\
                                    '_dataset_name_':self.dataset.name},
                                       {'postcode':'','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'postcode':'missing','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                       {'postcode':'678','_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (self.agree_weight > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(self.agree_weight)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 > similar3), \
           'Wrong similar3 weight (should be smaller than similar2 weight '+ \
           str(similar2)+'): '+ str(similar3)

    assert (similar2 > self.disagree_weight), \
           'Wrong similar2 weight (should be larger than disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar2)

    assert (similar3 == self.disagree_weight), \
           'Wrong similar3 weight (should be equal to disagree weight '+ \
           str(self.disagree_weight)+'): '+ str(similar3)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testDate(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorDate comparison routine"""

    date_compare = \
      comparison.FieldComparatorDate(fields_a=['dob','mob','yob'],
                                     fields_b=['dob','mob','yob'],
                                     m_probability_day =   self.m_prob,
                                     u_probability_day =   self.u_prob,
                                     m_probability_month = (self.m_prob-0.01),
                                     u_probability_month = (self.u_prob+0.01),
                                     m_probability_year =  (self.m_prob-0.04),
                                     u_probability_year =  (self.u_prob+0.04),
                                     missing_weight=self.miss_weight,
                                     max_day_a_before_b=10,
                                     max_day_b_before_a=2)

    agree_w = mymath.log2(self.m_prob / self.u_prob) + \
              mymath.log2((self.m_prob-0.01) / (self.u_prob+0.01)) + \
              mymath.log2((self.m_prob-0.04) / (self.u_prob+0.04))

    disagree_w = mymath.log2((1-self.m_prob) / (1-self.u_prob)) + \
                 mymath.log2((1-(self.m_prob-0.01)) / (1-(self.u_prob+0.01)))+\
                 mymath.log2((1-(self.m_prob-0.04)) / (1-(self.u_prob+0.04)))


    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [date_compare])

    agree = rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'dob':29,'mob':05,'yob':1968, \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'dob':23,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'dob':11,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    similar3 = rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':28,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar3 = float(similar3[4])

    similar4 = rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':26,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar4 = float(similar4[4])

    similar5 = rec_comparator.compare({'dob':21,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar5 = float(similar5[4])

    disagree = rec_comparator.compare({'dob':21,'mob':3,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1986, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'dob':'','mob':'','yob':'', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'dob':'missing','mob':5,'yob':1968, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'dob':29,'mob':5,'yob':1978, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == agree_w), \
           'Wrong agreement weight (should be: '+str(agree_w)+ \
           '): '+str(agree)

    assert (disagree == disagree_w), \
           'Wrong disagreement weight (should be: '+ \
           str(disagree_w)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (agree_w > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(agree_w)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar5 > similar2), \
           'Wrong similar2 weight (should be smaller than similar5 weight '+ \
           str(similar5)+'): '+ str(similar2)

    assert (similar1 > similar5), \
           'Wrong similar5 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar5)

    assert (similar2 == disagree_w), \
           'Wrong similar2 weight (should be equal to disagree weight '+ \
           str(disagree_w)+'): '+ str(similar2)

    assert (agree_w > similar3), \
           'Wrong similar3 weight (should be smaller than agree weight '+ \
           str(agree_w)+'): '+ str(similar3)

    assert (similar3 > similar4), \
           'Wrong similar4 weight (should be smaller than similar3 weight '+ \
           str(similar3)+'): '+ str(similar4)

    assert (similar4 == disagree_w), \
           'Wrong similar4 weight (should be equal to disagree weight '+ \
           str(disagree_w)+'): '+ str(similar4)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testAge(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorAge comparison routine"""

    age_compare = \
      comparison.FieldComparatorAge(fields_a=['dob','mob','yob'],
                                    fields_b=['dob','mob','yob'],
                                    m_probability_day =   self.m_prob,
                                    u_probability_day =   self.u_prob,
                                    m_probability_month = (self.m_prob-0.01),
                                    u_probability_month = (self.u_prob+0.01),
                                    m_probability_year =  (self.m_prob-0.04),
                                    u_probability_year =  (self.u_prob+0.04),
                                    missing_weight=self.miss_weight,
                                    max_perc_diff=5.0)

    agree_w = mymath.log2(self.m_prob / self.u_prob) + \
              mymath.log2((self.m_prob-0.01) / (self.u_prob+0.01)) + \
              mymath.log2((self.m_prob-0.04) / (self.u_prob+0.04))

    disagree_w = mymath.log2((1-self.m_prob) / (1-self.u_prob)) + \
                 mymath.log2((1-(self.m_prob-0.01)) / (1-(self.u_prob+0.01)))+\
                 mymath.log2((1-(self.m_prob-0.04)) / (1-(self.u_prob+0.04)))

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [age_compare])

    agree = rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'dob':29,'mob':05,'yob':1968, \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'dob':23,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'dob':11,'mob':5,'yob':1986, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1986, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    similar3 = rec_comparator.compare({'dob':29,'mob':1,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':28,'mob':9,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar3 = float(similar3[4])

    similar4 = rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':26,'mob':5,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar4 = float(similar4[4])

    disagree = rec_comparator.compare({'dob':21,'mob':3,'yob':1968, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'dob':29,'mob':5,'yob':1986, \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'dob':29,'mob':5,'yob':1968, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'dob':'','mob':'','yob':'', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'dob':'missing','mob':5,'yob':1968, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'dob':29,'mob':5,'yob':1978, \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == agree_w), \
           'Wrong agreement weight (should be: '+str(agree_w)+ \
           '): '+str(agree)

    assert (disagree == disagree_w), \
           'Wrong disagreement weight (should be: '+ \
           str(disagree_w)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (agree_w > similar1), \
           'Wrong similar1 weight (should be smaller than agree weight '+ \
           str(agree_w)+'): '+ str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight '+ \
           str(similar1)+'): '+ str(similar2)

    assert (similar2 > disagree_w), \
           'Wrong similar2 weight (should be largere than disagree weight '+ \
           str(disagree_w)+'): '+ str(similar2)

    assert (agree_w > similar3), \
           'Wrong similar3 weight (should be smaller than agree weight '+ \
           str(agree_w)+'): '+ str(similar3)

    assert (similar4 > similar3), \
           'Wrong similar3 weight (should be smaller than similar4 weight '+ \
           str(similar4)+'): '+ str(similar3)

    assert (similar4 > disagree_w), \
           'Wrong similar4 weight (should be larger than disagree weight '+ \
           str(disagree_w)+'): '+ str(similar4)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testTime(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorTime comparison routine"""

    time_compare = \
      comparison.FieldComparatorTime(fields_a='time',
                                     fields_b='time',
                                     m_prob=self.m_prob,
                                     u_prob=self.u_prob,
                                     missing_weight=self.miss_weight,
                                     max_time_a_before_b=60,
                                     max_time_b_before_a=10)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [time_compare])

    agree = rec_comparator.compare({'time':'11:59', \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'time':'1159', \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'time':'0838', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'time':'09:38', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'time':'0838', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'time':'10:38', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    similar3 = rec_comparator.compare({'time':'0942', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'time':'0938', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar3 = float(similar3[4])

    similar4 = rec_comparator.compare({'time':'17:52', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'time':'17:38', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar4 = float(similar4[4])

    disagree = rec_comparator.compare({'time':'23:58', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'time':'01:10', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'time':'1111', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'time':'', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'time':'missing', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'time':'0000', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (similar1 < self.agree_weight), \
           'Wrong similar1 weight (should be smaller than agree weight: '+ \
           str(self.agree_weight)+'): '+str(similar1)
    assert (similar1 > self.disagree_weight), \
           'Wrong similar1 weight (should be larger than disagree weight: '+ \
           str(self.disagree_weight)+'): '+str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight: '+ \
           str(similar1)+'): '+str(similar2)
    assert (similar2 == self.disagree_weight), \
           'Wrong similar2 weight (should be equal to disagree weight: '+ \
           str(self.disagree_weight)+'): '+str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


  def testDistance(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test the FieldComparatorDistance comparison routine"""

    pc_dict = lookup.GeocodeLookupTable(name='geocode test')

    pc_dict['2000'] = [151.20710, -33.87060]
    pc_dict['2007'] = [151.19761, -33.87849]
    pc_dict['2008'] = [151.19655, -33.88829]
    pc_dict['2009'] = [151.19395, -33.87014]
    pc_dict['2010'] = [151.21524, -33.88299]
    pc_dict['2011'] = [151.22422, -33.87105]

    # 2000 ->
    #   2007: 1.240km
    #   2008: 2.195km
    #   2009: 1.215km
    #   2010: 1.570km
    #   2011: 1.582km

    dist_compare = \
      comparison.FieldComparatorDistance(fields_a='postcode',
                                         fields_b='postcode',
                                         m_prob=self.m_prob,
                                         u_prob=self.u_prob,
                                         missing_weight=self.miss_weight,
                                         max_distance=2.0,
                                         geocode_table=pc_dict)

    rec_comparator = comparison.RecordComparator(self.dataset, self.dataset,
                                                 [dist_compare])

    agree = rec_comparator.compare({'postcode':'2000', \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name},
                                   {'postcode':'2000', \
                                    '_rec_num_':0, \
                                    '_dataset_name_':self.dataset.name})
    agree = float(agree[4])

    similar1 = rec_comparator.compare({'postcode':'2000', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'postcode':'2007', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar1 = float(similar1[4])

    similar2 = rec_comparator.compare({'postcode':'2011', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'postcode':'2000', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    similar2 = float(similar2[4])

    disagree = rec_comparator.compare({'postcode':'2000', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name},
                                      {'postcode':'2008', \
                                       '_rec_num_':0, \
                                       '_dataset_name_':self.dataset.name})
    disagree = float(disagree[4])

    missing1 =  rec_comparator.compare({'postcode':'2042', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'postcode':'2000', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing1 = float(missing1[4])

    missing2 =  rec_comparator.compare({'postcode':'2042', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'postcode':'', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing2 = float(missing2[4])

    missing3 =  rec_comparator.compare({'postcode':'missing', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name},
                                       {'postcode':'0000', \
                                        '_rec_num_':0, \
                                        '_dataset_name_':self.dataset.name})
    missing3 = float(missing3[4])

    assert (agree == self.agree_weight), \
           'Wrong agreement weight (should be: '+str(self.agree_weight)+ \
           '): '+str(agree)

    assert (disagree == self.disagree_weight), \
           'Wrong disagreement weight (should be: '+ \
           str(self.disagree_weight)+'): '+str(disagree)

    assert (agree > disagree), \
           'Agreement weight ('+str(agree)+') is not larger than disagree'+ \
           ' weight: '+str(disagree)

    assert (similar1 < self.agree_weight), \
           'Wrong similar1 weight (should be smaller than agree weight: '+ \
           str(self.agree_weight)+'): '+str(similar1)
    assert (similar1 > self.disagree_weight), \
           'Wrong similar1 weight (should be larger than disagree weight: '+ \
           str(self.disagree_weight)+'): '+str(similar1)

    assert (similar1 > similar2), \
           'Wrong similar2 weight (should be smaller than similar1 weight: '+ \
           str(similar1)+'): '+str(similar2)
    assert (similar2 > self.disagree_weight), \
           'Wrong similar2 weight (should be larger than disagree weight: '+ \
           str(self.disagree_weight)+'): '+str(similar2)

    assert (missing1 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing1)

    assert (missing2 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing2)

    assert (missing3 == self.miss_weight), \
           'Wrong missing weight (should be: '+ \
           str(self.miss_weight)+'): '+str(missing3)

    assert (agree > missing1), \
           'Agreement weight ('+str(agree)+') is not larger than missing1'+ \
           ' weight: '+str(missing1)


# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
