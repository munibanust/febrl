# =============================================================================
# project-linkage.py - Configuration for a record linkage project.
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
# The Original Software is "project-linkage.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module project-linkage.py - Configuration for a record linkage project

   Briefly, what needs to be defined for a linkage project is:
   - A Febrl object, a project, plus a project logger
   - Two input data set
   - Two corresponding temporary data set (with 'readwrite' access)
   - Lookup tables to be used
   - Standardisers for names, addresses and dates for both input data sets
     (these standardisers can differ for the two input data sets)
   - Field comparator functions and a record comparator
   - Two blocking indexes (they should be the same, i.e. lead to similar
     blocking variables)
   - A classifier

   and then the 'link' method can be called.

   For more information see chapter

   "Configuration and Running Febrl using a Module derived from 'project.py'"

   in the Febrl manual.

   This project module will standardised and then link the example data sets
   'dataset4a.csv' and 'dataset4b.csv' given in the 'dbgen' directory.
"""

# =============================================================================
# Imports go here

import sys
import time

from febrl import *            # Main Febrl classes
from dataset import *          # Data set routines
from standardisation import *  # Standardisation routines
from comparison import *       # Comparison functions
from lookup import *           # Look-up table routines
from indexing import *         # Indexing and blocking routines
from simplehmm import *        # Hidden Markov model (HMM) routines
from classification import *   # Classifiers for weight vectors

# =============================================================================
# Set up Febrl and create a new project (or load a saved project)

myfebrl = Febrl(description = 'Example linkage Febrl instance',
                 febrl_path = '.')

myproject = myfebrl.new_project(name = 'example-link',
                         description = 'Link example data sets 4',
                           file_name = 'example-linkage.fbr',
                          block_size = 1000,
                      parallel_write = 'host')

# =============================================================================
# Define a project logger

mylog = ProjectLog(file_name = 'example-linkage.log',
                     project = myproject,
                   log_level = 1,
               verbose_level = 1,
                   clear_log = True,
                     no_warn = False,
              parallel_print = 'host')

# =============================================================================
# Define original input data set(s)
# Two data sets are needed for linkage

indata_a = DataSetCSV(name = 'example4a-in',
               description = 'Example data set 4a',
               access_mode = 'read',
              header_lines = 1,
               file_name = './dbgen/dataset4a.csv',
                  fields = {'rec_id':0,
                            'given_name':1,
                            'surname':2,
                            'street_num':3,
                            'address_part_1':4,
                            'address_part_2':5,
                            'suburb':6,
                            'postcode':7,
                            'state':8,
                            'date_of_birth':9,
                            'soc_sec_id':10},
          fields_default = '',
            strip_fields = True,
          missing_values = ['','missing'])

indata_b = DataSetCSV(name = 'example4b-in',
               description = 'Example data set 4b',
               access_mode = 'read',
              header_lines = 1,
              write_header = True,
               file_name = './dbgen/dataset4b.csv',
                  fields = {'rec_id':0,
                            'given_name':1,
                            'surname':2,
                            'street_num':3,
                            'address_part_1':4,
                            'address_part_2':5,
                            'suburb':6,
                            'postcode':7,
                            'state':8,
                            'date_of_birth':9,
                            'soc_sec_id':10},
          fields_default = '',
            strip_fields = True,
          missing_values = ['','missing'])


# =============================================================================
# Define temporary data set(s) (one per input data set)
# Commented lines are only needed for the disk based Shelve data set

tmpdata_a = DataSetMemory(name = 'example4a-tmp',
#tmpdata_a = DataSetShelve(name = 'example4a-tmp',
#                     file_name = './example4a-shelve',
#                         clear = True,
                   description = 'Temporary example 4a data set',
                   access_mode = 'readwrite',
                        fields = {'title':1,
                                  'gender_guess':2,
                                  'given_name':3,
                                  'alt_given_name':4,
                                  'surname':5,
                                  'alt_surname':6,
                                  'wayfare_number':7,
                                  'wayfare_name':8,
                                  'wayfare_qualifier':9,
                                  'wayfare_type':10,
                                  'unit_number':11,
                                  'unit_type':12,
                                  'property_name':13,
                                  'institution_name':14,
                                  'institution_type':15,
                                  'postaddress_number':16,
                                  'postaddress_type':17,
                                  'locality_name':18,
                                  'locality_qualifier':19,
                                  'postcode':20,
                                  'territory':21,
                                  'country':22,
                                  'dob_day':23,
                                  'dob_month':24,
                                  'dob_year':25,
# The following are output fields that are passed without standardisation
                                  'rec_id':0,
                                  'soc_sec_id':26,
# The last output field contains the probability of the address HMM
                                  'address_hmm_prob':27,
                                 },
                missing_values = ['','missing'])

tmpdata_b = DataSetMemory(name = 'example4b-tmp',
#tmpdata_b = DataSetShelve(name = 'example4b-tmp',
#                     file_name = './example4b-shelve',
#                         clear = True,
                   description = 'Temporary example 4b data set',
                   access_mode = 'readwrite',
                        fields = {'title':1,
                                  'gender_guess':2,
                                  'given_name':3,
                                  'alt_given_name':4,
                                  'surname':5,
                                  'alt_surname':6,
                                  'wayfare_number':7,
                                  'wayfare_name':8,
                                  'wayfare_qualifier':9,
                                  'wayfare_type':10,
                                  'unit_number':11,
                                  'unit_type':12,
                                  'property_name':13,
                                  'institution_name':14,
                                  'institution_type':15,
                                  'postaddress_number':16,
                                  'postaddress_type':17,
                                  'locality_name':18,
                                  'locality_qualifier':19,
                                  'postcode':20,
                                  'territory':21,
                                  'country':22,
                                  'dob_day':23,
                                  'dob_month':24,
                                  'dob_year':25,
# The following are output fields that are passed without standardisation
                                  'rec_id':0,
                                  'soc_sec_id':26,
# The last output field contains the probability of the address HMM
                                  'address_hmm_prob':27,
                                 },
                missing_values = ['','missing'])

# =============================================================================
# Define and load lookup tables

name_lookup_table = TagLookupTable(name = 'Name lookup table',
                                default = '')

name_lookup_table.load(['./data/givenname_f.tbl',
                        './data/givenname_m.tbl',
                        './data/name_prefix.tbl',
                        './data/name_misc.tbl',
                        './data/saints.tbl',
                        './data/surname.tbl',
                        './data/title.tbl'])

name_correction_list = CorrectionList(name = 'Name correction list')

name_correction_list.load('./data/name_corr.lst')

surname_freq_table = FrequencyLookupTable(name = 'Surname frequency table',
                                       default = 1)

surname_freq_table.load('./data/surname_nsw_freq.csv')

address_lookup_table = TagLookupTable(name = 'Address lookup table',
                                   default = '')

address_lookup_table.load(['./data/country.tbl',
                           './data/address_misc.tbl',
                           './data/address_qual.tbl',
                           './data/institution_type.tbl',
                           './data/locality_name_act.tbl',
                           './data/locality_name_nsw.tbl',
                           './data/post_address.tbl',
                           './data/postcode_act.tbl',
                           './data/postcode_nsw.tbl',
                           './data/saints.tbl',
                           './data/territory.tbl',
                           './data/unit_type.tbl',
                           './data/wayfare_type.tbl'])

address_correction_list = CorrectionList(name = 'Address correction list')

address_correction_list.load('./data/address_corr.lst')

pc_geocode_table = GeocodeLookupTable(name = 'NSW postcode locations',
                                   default = [])

pc_geocode_table.load('./data/postcode_nsw_geocode.csv')

# =============================================================================
# Define and load hidden Markov models (HMMs)

name_states = ['titl','baby','knwn','andor','gname1','gname2','ghyph',
               'gopbr','gclbr','agname1','agname2','coma','sname1','sname2',
               'shyph','sopbr','sclbr','asname1','asname2','pref1','pref2',
               'rubb']
name_tags = ['NU','AN','TI','PR','GF','GM','SN','ST','SP','HY','CO','NE','II',
             'BO','VB','UN','RU']

myname_hmm = hmm('Name HMM', name_states, name_tags)
myname_hmm.load_hmm('./hmm/name-absdiscount.hmm')
# myname_hmm.load_hmm('./hmm/name.hmm')
# myname_hmm.load_hmm('./hmm/name-laplace.hmm')

address_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty','prna1',
                  'prna2','inna1','inna2','inty','panu','paty','hyph','sla',
                  'coma','opbr','clbr','loc1','loc2','locql','pc','ter1',
                  'ter2','cntr1','cntr2','rubb']
address_tags = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT','LQ','WT',
                'WN','UT','HY','SL','CO','VB','PA','UN','RU']

myaddress_hmm = hmm('Address HMM', address_states, address_tags)
myaddress_hmm.load_hmm('./hmm/address-absdiscount.hmm')
# myaddress_hmm.load_hmm('./hmm/address.hmm')
# myaddress_hmm.load_hmm('./hmm/address-laplace.hmm')

# =============================================================================
# Define a list of date parsing format strings

date_parse_formats = ['%d %m %Y',   # 24 04 2002  or  24 4 2002
                      '%d %B %Y',   # 24 Apr 2002 or  24 April 2002
                      '%m %d %Y',   # 04 24 2002  or  4 24 2002
                      '%B %d %Y',   # Apr 24 2002 or  April 24 2002
                      '%Y %m %d',   # 2002 04 24  or  2002 4 24
                      '%Y %B %d',   # 2002 Apr 24 or  2002 April 24
                      '%Y%m%d',     # 20020424                   ISO standard
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
                     ]

# =============================================================================
# Define standardisers for dates (one for each data set needed)

dob_std_a = DateStandardiser(name = 'DOB-std-a',
                      description = 'Date of birth a standardiser',
                     input_fields = 'date_of_birth',
                    output_fields = ['dob_day','dob_month', 'dob_year'],
                    parse_formats = date_parse_formats)

dob_std_b = DateStandardiser(name = 'DOB-std-b',
                      description = 'Date of birth b standardiser',
                     input_fields = 'date_of_birth',
                    output_fields = ['dob_day','dob_month', 'dob_year'],
                    parse_formats = date_parse_formats)

# =============================================================================
# Define standardisers for names based on rules (one for each data set needed)

name_rules_std_a = NameRulesStandardiser(name = 'Name-Rules-a',
                                 input_fields = ['given_name','surname'],
                                output_fields = ['title',
                                                 'gender_guess',
                                                 'given_name',
                                                 'alt_given_name',
                                                 'surname',
                                                 'alt_surname'],
                               name_corr_list = name_correction_list,
                               name_tag_table = name_lookup_table,
                                  male_titles = ['mr'],
                                female_titles = ['ms'],
                              field_separator = ' ',
                             check_word_spill = True)

name_rules_std_b = NameRulesStandardiser(name = 'Name-Rules-b',
                                 input_fields = ['given_name','surname'],
                                output_fields = ['title',
                                                 'gender_guess',
                                                 'given_name',
                                                 'alt_given_name',
                                                 'surname',
                                                 'alt_surname'],
                               name_corr_list = name_correction_list,
                               name_tag_table = name_lookup_table,
                                  male_titles = ['mr'],
                                female_titles = ['ms'],
                              field_separator = ' ',
                             check_word_spill = True)

# =============================================================================
# Define standardisers for names based on HMM (one for each data set needed)

name_hmm_std_a = NameHMMStandardiser(name = 'Name-HMM-a',
                             input_fields = ['given_name','surname'],
                            output_fields = ['title',
                                             'gender_guess',
                                             'given_name',
                                             'alt_given_name',
                                             'surname',
                                             'alt_surname'],
                           name_corr_list = name_correction_list,
                           name_tag_table = name_lookup_table,
                              male_titles = ['mr'],
                            female_titles = ['ms'],
                                 name_hmm = myname_hmm,
                          field_separator = ' ',
                         check_word_spill = True)

name_hmm_std_b = NameHMMStandardiser(name = 'Name-HMM-b',
                             input_fields = ['given_name','surname'],
                            output_fields = ['title',
                                             'gender_guess',
                                             'given_name',
                                             'alt_given_name',
                                             'surname',
                                             'alt_surname'],
                           name_corr_list = name_correction_list,
                           name_tag_table = name_lookup_table,
                              male_titles = ['mr'],
                            female_titles = ['ms'],
                                 name_hmm = myname_hmm,
                          field_separator = ' ',
                         check_word_spill = True)

# =============================================================================
# Define standardisers for addresses based on HMM (one for each data set)

address_hmm_std_a = AddressHMMStandardiser(name = 'Address-HMM-a',
                                   input_fields = ['street_num',
                                                   'address_part_1',
                                                   'address_part_2','suburb',
                                                   'postcode', 'state'],
                                  output_fields = ['wayfare_number',
                                                   'wayfare_name',
                                                   'wayfare_qualifier',
                                                   'wayfare_type',
                                                   'unit_number',
                                                   'unit_type',
                                                   'property_name',
                                                   'institution_name',
                                                   'institution_type',
                                                   'postaddress_number',
                                                   'postaddress_type',
                                                   'locality_name',
                                                   'locality_qualifier',
                                                   'postcode',
                                                   'territory',
                                                   'country',
                                                   'address_hmm_prob'],
                              address_corr_list = address_correction_list,
                              address_tag_table = address_lookup_table,
                                    address_hmm = myaddress_hmm)

address_hmm_std_b = AddressHMMStandardiser(name = 'Address-HMM-b',
                                   input_fields = ['street_num',
                                                   'address_part_1',
                                                   'address_part_2','suburb',
                                                   'postcode', 'state'],
                                  output_fields = ['wayfare_number',
                                                   'wayfare_name',
                                                   'wayfare_qualifier',
                                                   'wayfare_type',
                                                   'unit_number',
                                                   'unit_type',
                                                   'property_name',
                                                   'institution_name',
                                                   'institution_type',
                                                   'postaddress_number',
                                                   'postaddress_type',
                                                   'locality_name',
                                                   'locality_qualifier',
                                                   'postcode',
                                                   'territory',
                                                   'country',
                                                   'address_hmm_prob'],
                              address_corr_list = address_correction_list,
                              address_tag_table = address_lookup_table,
                                    address_hmm = myaddress_hmm)

# =============================================================================
# Define pass field standardisers for all fields that should be passed from
# the input to the output data set without any cleaning or standardisdation.

pass_fields_a = PassFieldStandardiser(name = 'Pass fields-a',
                              input_fields = ['rec_id', 'soc_sec_id'],
                             output_fields = ['rec_id', 'soc_sec_id'])

pass_fields_b = PassFieldStandardiser(name = 'Pass fields-b',
                              input_fields = ['rec_id', 'soc_sec_id'],
                             output_fields = ['rec_id', 'soc_sec_id'])

# =============================================================================
# Define record standardisers (one for each data set)

comp_stand_a = [dob_std_a, name_rules_std_a, address_hmm_std_a, pass_fields_a]

comp_stand_b = [dob_std_b, name_rules_std_b, address_hmm_std_b, pass_fields_b]

# The HMM based name standardisation is not used in the above standardisers,
# uncomment the lines below (and comment the ones above) to use HMM
# standardisation for names.
#
#comp_stand_a = [dob_std_a, name_hmm_std_a, address_hmm_std_a, pass_fields_a]

#comp_stand_b = [dob_std_b, name_hmm_std_b, address_hmm_std_b, pass_fields_b]

example_standardiser_a = RecordStandardiser(name = 'Example-std-a',
                                     description = 'Example standardiser a',
                                   input_dataset = indata_a,
                                  output_dataset = tmpdata_a,
                                        comp_std = comp_stand_a)

example_standardiser_b = RecordStandardiser(name = 'Example-std-b',
                                     description = 'Example standardiser b',
                                   input_dataset = indata_b,
                                  output_dataset = tmpdata_b,
                                        comp_std = comp_stand_b)

# =============================================================================
# Define blocking indexes (one per temporary data set)

myblock_def_a = [[('surname','dmetaphone', 4),('dob_year','direct')],
                 [('given_name','truncate', 3), ('postcode','direct')],
                 [('locality_name','nysiis'),('dob_month','direct')],
                ]

myblock_def_b = [[('surname','dmetaphone', 4),('dob_year','direct')],
                 [('given_name','truncate', 3), ('postcode','direct')],
                 [('locality_name','nysiis'),('dob_month','direct')],
                ]

# Define one or more indexes (to be used in the classifier further below)

example_index_a = BlockingIndex(name = 'Index-blocking-a',
                             dataset = tmpdata_a,
                           index_def = myblock_def_a)

example_sorting_index_a = SortingIndex(name = 'Index-sorting-a',
                                    dataset = tmpdata_a,
                                  index_def = myblock_def_a,
                                window_size = 3)

example_bigram_index_a = BigramIndex(name = 'Index-bigram-a',
                                  dataset = tmpdata_a,
                                index_def = myblock_def_a,
                                threshold = 0.75)

example_index_b = BlockingIndex(name = 'Index-blocking-b',
                             dataset = tmpdata_b,
                           index_def = myblock_def_b)

example_sorting_index_b = SortingIndex(name = 'Index-sorting-b',
                                    dataset = tmpdata_b,
                                  index_def = myblock_def_b,
                                window_size = 3)

example_bigram_index_b = BigramIndex(name = 'Index-bigram-b',
                                  dataset = tmpdata_b,
                                index_def = myblock_def_b,
                                threshold = 0.75)

# =============================================================================
# Define comparison functions for linkage

given_name_nysiis = FieldComparatorEncodeString(name = 'Given name NYSIIS',
                                            fields_a = 'given_name',
                                            fields_b = 'given_name',
                                              m_prob = 0.95,
                                              u_prob = 0.001,
                                      missing_weight = 0.0,
                                       encode_method = 'nysiis',
                                             reverse = False)

surname_dmetaphone = FieldComparatorEncodeString(name = 'Surname D-Metaphone',
                                             fields_a = 'surname',
                                             fields_b = 'surname',
                                               m_prob = 0.95,
                                               u_prob = 0.001,
                                       missing_weight = 0.0,
                                        encode_method = 'dmetaphone',
                                              reverse = False)

wayfare_name_winkler = FieldComparatorApproxString(name = 'Wayfare name ' + \
                                                          'Winkler',
                                               fields_a = 'wayfare_name',
                                               fields_b = 'wayfare_name',
                                                 m_prob = 0.95,
                                                 u_prob = 0.001,
                                         missing_weight = 0.0,
                                         compare_method = 'winkler',
                                       min_approx_value = 0.7)

locality_name_key = FieldComparatorKeyDiff(name = 'Locality name key diff',
                                       fields_a = 'locality_name',
                                       fields_b = 'locality_name',
                                         m_prob = 0.95,
                                         u_prob = 0.001,
                                 missing_weight = 0.0,
                                   max_key_diff = 2)

postcode_distance = FieldComparatorDistance(name = 'Postcode distance',
                                        fields_a = 'postcode',
                                        fields_b = 'postcode',
                                          m_prob = 0.95,
                                          u_prob = 0.001,
                                  missing_weight = 0.0,
                                   geocode_table = pc_geocode_table,
                                    max_distance = 50.0)

age = FieldComparatorAge(name = 'Age',
                     fields_a = ['dob_day','dob_month', 'dob_year'],
                     fields_b = ['dob_day','dob_month', 'dob_year'],
            m_probability_day = 0.95,
            u_probability_day = 0.03333,
          m_probability_month = 0.95,
          u_probability_month = 0.083,
           m_probability_year = 0.95,
           u_probability_year = 0.01,
                max_perc_diff = 10.0,
                     fix_date = 'today')

field_comparisons = [given_name_nysiis, surname_dmetaphone, \
                     wayfare_name_winkler, locality_name_key, \
                     postcode_distance, age]

example_comparator = RecordComparator(tmpdata_a, tmpdata_b, field_comparisons)

# =============================================================================
# Define a classifier for classifying the matching vectors

example_fs_classifier = FellegiSunterClassifier(name = 'Fellegi and Sunter',
                                           dataset_a = tmpdata_a,
                                           dataset_b = tmpdata_b,
                                     lower_threshold = 0.0,
                                     upper_threshold = 10.0)

example_flex_classifier = FlexibleClassifier(name = 'Example flex classifier',
                                        dataset_a = tmpdata_a,
                                        dataset_b = tmpdata_b,
                                  lower_threshold = 0.0,
                                  upper_threshold = 10.0,
                                        calculate = [('avrg', [0,1]),
                                                     ('max',  [2,3,4]),
                                                     ('min',  [5])],
                                      final_funct = 'avrg')

# =============================================================================
# Start the linkage task

myproject.link(input_dataset_a = indata_a,
               input_dataset_b = indata_b,
                 tmp_dataset_a = tmpdata_a,
                 tmp_dataset_b = tmpdata_b,
            rec_standardiser_a = example_standardiser_a,
            rec_standardiser_b = example_standardiser_b,
              blocking_index_a = example_index_a,
              blocking_index_b = example_index_b,
                rec_comparator = example_comparator,
                    classifier = example_fs_classifier,
                first_record_a = 0,
              number_records_a = 5000,
                first_record_b = 0,
              number_records_b = 5000,
              output_histogram = 'link-example-histogram.res',
       output_rec_pair_details = 'link-example-details.res',
       output_rec_pair_weights = 'link-example-weights.res',
              output_threshold = 10.0,
             output_assignment = 'one2one')

# =============================================================================

myfebrl.finalise()

# =============================================================================
