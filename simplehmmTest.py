# =============================================================================
# simpleTest.py - Test module for simplehmm.py
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
# The Original Software is "simplehmmTest.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module simpleTest.py - Test module for simplehmm.py.
"""

# -----------------------------------------------------------------------------

import unittest
import simplehmm

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    self.delta = 0.000001  # Account for floating-point rounding errors

    self.states = ['title','givenname','surname']
    self.observ = ['TI','GM','GF','SN','UN']

    self.train_data = [[('title','TI'),('givenname','GF'),('surname','SN')], \
                      [('givenname','GM'),('surname','UN')], \
                      [('title','UN'),('givenname','GM'),('surname','UN')], \
                      [('title','TI'),('givenname','SN'),('surname','SN')], \
                      [('givenname','GM'),('surname','SN')], \
                      [('title','TI'),('givenname','GF'),('surname','SN')], \
                      [('title','TI'),('surname','SN'),('givenname','GM')], \
                      [('title','TI'),('givenname','GM'),('givenname','UN'), \
                       ('surname','SN')], \
                      [('surname','UN'),('givenname','UN')], \
                      [('givenname','GF'),('givenname','GF'),('surname','SN')]]

    self.test_data = [['TI','GM','SN'], \
                      ['UN','SN'], \
                      ['TI','UN','UN'], \
                      ['TI','GF','UN'], \
                      ['UN','UN','UN','UN'], \
                      ['TI','GM','UN','SN'], \
                      ['SN','SN'], \
                      ['UN','UN'], \
                      ['TI','GM'], \
                      ['TI','GF'], \
                      ['TI','SN'], \
                      ['TI','UN']]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testHMM(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test basic HMM functionality"""

    hmm1 = simplehmm.hmm('Test HMM', self.states, self.observ)

    assert (hmm1.N == len(self.states)), \
           'Illegal number of states in HMM ('+str(hmm1.N)+'), should be: '+ \
           str(len(self.states))
    assert (len(hmm1.S_ind) == len(self.states)), \
           'Illegal number of states in HMM state dictionary ('+ \
           str(len(hmm1.S_ind))+'), should be: '+str(len(self.states))

    assert (hmm1.M == len(self.observ)), \
           'Illegal number of observations in HMM ('+str(hmm1.M)+ \
           '), should be: '+str(len(self.observ))
    assert (len(hmm1.O_ind) == len(self.observ)), \
           'Illegal number of observations in HMM observation dictionary ('+ \
           str(len(hmm1.O_ind))+'), should be: '+ str(len(self.observ))

    for i in range(hmm1.N):
      assert (hmm1.pi[i] == 0.0), \
             'Initial probability in HMM 1 is not 0.0 at location ['+ \
             str(i)+']: '+str(hmm1.pi[i])

      for j in range(hmm1.N):
         assert (hmm1.A[i][j] == 0.0), \
                'Transition probability in HMM 1 is not 0.0 at location ['+ \
                str(i)+','+str(j)+']: '+str(hmm1.A[i][j])
      for j in range(hmm1.M):
         assert (hmm1.B[i][j] == 0.0), \
                'Observation probability in HMM 1 is not 0.0 at location ['+ \
                str(i)+','+str(j)+']: '+str(hmm1.B[i][j])


    hmm1.train(self.train_data)
    hmm1.check_prob()
    hmm1.print_hmm()

    for i in range(hmm1.N):
      assert ((hmm1.pi[i] >= 0.0) and (hmm1.pi[i] <= 1.0)), \
             'Initial probability in HMM 1 is not between 0.0 and 1.0 at '+ \
             'location ['+str(i)+']: '+str(hmm1.pi[i])

      for j in range(hmm1.N):
         assert ((hmm1.A[i][j] >= 0.0) and (hmm1.A[i][j] <= 1.0)), \
                'Transition probability in HMM 1 is not between 0.0 and 1.0'+ \
                ' at location ['+str(i)+','+str(j)+']: '+str(hmm1.A[i][j])
      for j in range(hmm1.M):
         assert ((hmm1.B[i][j] >= 0.0) and (hmm1.B[i][j] <= 1.0)), \
                'Observation probability in HMM 1 is not between 0.0 and '+ \
                '1.0 at location ['+str(i)+','+str(j)+']: '+str(hmm1.B[i][j])

    for test_rec in self.test_data:
      [state_seq, seq_prob] = hmm1.viterbi(test_rec)

      for state in state_seq:
        assert (state in self.states), \
               'Returned state "'+state+'" not in tate list'
      assert ((seq_prob >= 0.0) and (seq_prob <= 1.0)), \
            'Sequence probability is not between 0.0 and 1.0:'+ str(seq_prob)


    hmm1.train(self.train_data,smoothing='laplace')
    hmm1.check_prob()
    hmm1.print_hmm()

    for i in range(hmm1.N):
      assert ((hmm1.pi[i] >= 0.0) and (hmm1.pi[i] <= 1.0)), \
             'Initial probability in HMM 1 is not between 0.0 and 1.0 at '+ \
             'location ['+str(i)+']: '+str(hmm1.pi[i])

      for j in range(hmm1.N):
         assert ((hmm1.A[i][j] >= 0.0) and (hmm1.A[i][j] <= 1.0)), \
                'Transition probability in HMM 1 is not between 0.0 and 1.0'+ \
                ' at location ['+str(i)+','+str(j)+']: '+str(hmm1.A[i][j])
      for j in range(hmm1.M):
         assert ((hmm1.B[i][j] >= 0.0) and (hmm1.B[i][j] <= 1.0)), \
                'Observation probability in HMM 1 is not between 0.0 and '+ \
                '1.0 at location ['+str(i)+','+str(j)+']: '+str(hmm1.B[i][j])

    for test_rec in self.test_data:
      [state_seq, seq_prob] = hmm1.viterbi(test_rec)

      for state in state_seq:
        assert (state in self.states), \
               'Returned state "'+state+'" not in tate list'
      assert ((seq_prob >= 0.0) and (seq_prob <= 1.0)), \
            'Sequence probability is not between 0.0 and 1.0:'+ str(seq_prob)


    hmm1.train(self.train_data,smoothing='absdiscount')
    hmm1.check_prob()
    hmm1.print_hmm()

    for i in range(hmm1.N):
      assert ((hmm1.pi[i] >= 0.0) and (hmm1.pi[i] <= 1.0)), \
             'Initial probability in HMM 1 is not between 0.0 and 1.0 at '+ \
             'location ['+str(i)+']: '+str(hmm1.pi[i])

      for j in range(hmm1.N):
         assert ((hmm1.A[i][j] >= 0.0) and (hmm1.A[i][j] <= 1.0)), \
                'Transition probability in HMM 1 is not between 0.0 and 1.0'+ \
                ' at location ['+str(i)+','+str(j)+']: '+str(hmm1.A[i][j])
      for j in range(hmm1.M):
         assert ((hmm1.B[i][j] >= 0.0) and (hmm1.B[i][j] <= 1.0)), \
                'Observation probability in HMM 1 is not between 0.0 and '+ \
                '1.0 at location ['+str(i)+','+str(j)+']: '+str(hmm1.B[i][j])

    for test_rec in self.test_data:
      [state_seq, seq_prob] = hmm1.viterbi(test_rec)

      for state in state_seq:
        assert (state in self.states), \
               'Returned state "'+state+'" not in tate list'
      assert ((seq_prob >= 0.0) and (seq_prob <= 1.0)), \
            'Sequence probability is not between 0.0 and 1.0:'+ str(seq_prob)


    hmm1.save_hmm('testhmm.hmm')
    hmm2 = hmm1

    hmm2 = simplehmm.hmm('Test2 HMM', ['dummy'], ['dummy'])

    hmm2.load_hmm('testhmm.hmm')

    assert (hmm1.N == hmm2.N), \
           'Loaded HMM has differnt number of states'
    assert (hmm1.M == hmm2.M), \
           'Loaded HMM has differnt number of observations'

    for i in range(hmm1.N):
      assert (abs(hmm1.pi[i]- hmm2.pi[i]) < self.delta), \
             'Initial probability in HMM 1 is different from HMM 2: '+ \
             str(hmm1.pi[i])+' / '+str(hmm2.pi[i])

      for j in range(hmm1.N):
         assert (abs(hmm1.A[i][j] - hmm2.A[i][j]) < self.delta), \
                'Transition probability in HMM 1 is different from HMM 2 '+ \
                'at location ['+str(i)+','+str(j)+']: '+str(hmm1.A[i][j])+ \
                ' / '+str(hmm2.A[i][j])

      for j in range(hmm1.M):
         assert (abs(hmm1.B[i][j] - hmm1.B[i][j]) < self.delta), \
                'Observation probability in HMM 1 is different from HMM 2 '+ \
                'at location ['+str(i)+','+str(j)+']: '+str(hmm1.B[i][j])+ \
                ' / '+str(hmm2.B[i][j])

# -----------------------------------------------------------------------------
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

  # The following code does the same as 'unittest.main()'
  #
  # mysuite = unittest.makeSuite(TestCase,'test')
  # testrunner = unittest.TextTestRunner(verbosity=1)
  # testrunner.run(mysuite)
