# =============================================================================
# config.py - Various system configuration settings and information
#
# Freely extensible biomedical record linkage (Febrl) Version 0.1
# See http://datamining.anu.edu.au/projects/linkage.html
#
# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.0
#
# The contents of this file are subject to the ANUOS License Version 1.0 (the
# "License"); you may not use this file except in compliance with the License.
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
# the specific language governing rights and limitations under the License.
# The Original Software is "config.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module config.py - Various system configuration settings and information

   This module contains system wide settings and information that is not user
   adjustable, but needed in various other modules.

   This module also loads the 'project.py' module and tests and processes all
   settings in there.
"""

# -----------------------------------------------------------------------------

import sys
import types

import tcsv  # Tim Churches slow but flexible CSV parser

import inout
import mymath
import simplehmm

# =============================================================================
# Load the project module as given as first command line argument.

if (len(sys.argv) < 2):
  print '***** Error: %s needs at least one argument:'% (sys.argv[0])
  print '*****        (project module missing)'
  raise Exception()

project_mod = sys.argv[1].strip()
if (project_mod[-3:] == '.py'):  # Strip off .py extension
  project_mod = project_mod[:-3]

exec('import '+project_mod)
exec('project = '+project_mod)

options = sys.argv[2:]  # Save command line options (after project argument)

# =============================================================================
# Process setting from project and make wrappers around them

verbose  = project.verbose
logging  = project.logging
log_file = project.log_file
nowarn   = project.nowarn
proc_ind = project.proc_ind

in_file_name = project.in_file_name
in_file_type = project.in_file_type
out_file_name = project.out_file_name
out_file_type = project.out_file_type
input_component = project.input_component
output_field = project.output_field
input_space_sep = project.input_space_sep
input_check_spilling = project.input_check_spilling
output_quote_character = project.output_quote_character
name_female_title = project.name_female_title
name_male_title = project.name_male_title
name_standard_method = project.name_standard_method
geoloc_standard_method = project.geoloc_standard_method

date_pivot_year = project.date_pivot_year
date_parse_formats = project.date_parse_formats
date_perc_fix_date = project.date_perc_fix_date
date_age_fix_date = project.date_age_fix_date
date_day_m_prob = project.date_day_m_prob
date_day_u_prob = project.date_day_u_prob
date_month_m_prob = project.date_month_m_prob
date_month_u_prob = project.date_month_u_prob
date_year_m_prob = project.date_year_m_prob
date_year_u_prob = project.date_year_u_prob
date_comp_max_subst = project.date_comp_max_subst
date_comp_max_trans = project.date_comp_max_trans
date_comp_max_day_before = project.date_comp_max_day_before
date_comp_max_day_after =  project.date_comp_max_day_after
date_comp_max_perc_before = project.date_comp_max_perc_before
date_comp_max_perc_after = project.date_comp_max_perc_after
date_linkage_weight_comp = project.date_linkage_weight_comp

# =============================================================================
# Content and number of current line in the input file that is processed

curr_line       = ''  # The original input string
curr_line_list  = []  # The input split into fields
curr_line_no = 0

# =============================================================================
# Number of warnings and number of corrected word spillings

num_warning = 0
num_word_spills = 0

# =============================================================================
# Definition of Hidden Markov Model states and observations

name_hmm_states = ['titl','baby','knwn','andor','gname1','gname2','ghyph', \
                   'gopbr','gclbr','agname1','agname2','coma','sname1', \
                   'sname2','shyph','sopbr','sclbr','asname1','asname2', \
                   'pref1','pref2','rubb']
name_hmm_obser  = ['NU','AN','TI','PR','GF','GM','SN','ST','SP','HY','CO', \
                   'NE','II','BO','VB','UN','RU']

geoloc_hmm_states = ['wfnu','wfna1','wfna2','wfql','wfty','unnu','unty', \
                     'prna1','prna2','inna1','inna2','inty','panu','paty', \
                     'hyph','sla','coma','opbr','clbr','loc1','loc2', \
                     'locql','pc','ter1','ter2','cntr1','cntr2','rubb']
geoloc_hmm_obser  = ['PC','N4','NU','AN','TR','CR','LN','ST','IN','IT', \
                     'LQ','WT','WN','UT','HY','SL','CO','VB','PA','UN', \
                     'RU']

# =============================================================================
# Dictionary of month name abbreviations (used in date.str2date() routine)

month_abbrev_dict = {'jan':1, 'feb':2, 'mar':3, 'apr':4, 'may':5, 'jun':6, \
                     'jul':7, 'aug':8, 'sep':9, 'oct':10, 'nov':11, 'dec': 12}

# =============================================================================
# If Hidden Markov Model standardisation methods are activate load HMM(s)

if (project.name_standard_method == 'hmm'):
  name_hmm = simplehmm.hmm([],[])  # Create new empty HMM object
  name_hmm.load_hmm(project.name_hmm_file_name)

if (project.geoloc_standard_method == 'hmm'):
  geoloc_hmm = simplehmm.hmm([],[])  # Create new empty HMM object
  geoloc_hmm.load_hmm(project.geoloc_hmm_file_name)

# =============================================================================
# List of all supported data file types
#
# File type names must have a length of 3 characters, or 4 characters if the
# file type is quoted (in which case the last character must be a 'Q')
#
# Currently supported file types are:
#   CSV  - Comma separated values, fields separated by commas
#   CSVQ - Comma separated values, where each field starts and ends with
#          a quote character
#   TAB  - Tabulator separated values, fields separated by commas
#   TABQ - Tabulator separated values, where each field starts and ends with
#          a quote character
#   COL  - Column wise, fields within specific column ranges
#
#   A database access file type (SQL) is planned to be included in a future
#   release of this software.
#
file_types = ['CSV','CSVQ','TAB','TABQ','COL']

    ################################
    # Add later: SQL database access
    ################################

# -----------------------------------------------------------------------------
# Test input and output file types as defined in project module

if (project.in_file_type not in file_types):
  print '***** Error: Illegal input file type:', project.in_file_type
  print '*****        File type must be in:', file_types
  raise Exception()

if (project.out_file_type not in file_types):
  print '***** Error: Illegal output file type:', project.out_file_type
  print '*****        File type must be in:', file_types
  raise Exception()

# -----------------------------------------------------------------------------
# Check if input and output file are not the same
#
if (project.in_file_name == project.out_file_name):
  print '***** Error: Input and output files must differ'
  print '*****        Input file name: ', project.in_file
  print '*****        Output file name:', project.out_file
  raise Exception()

# -----------------------------------------------------------------------------
# Check if definition of input components is correct with file types
#
input_values = input_component.values()
input_len = -1  # Length of the input (either in number of fields (CSV and TAB
                # files) or in characters (COL files)

output_keys = output_field.keys() # Check if 'original_input' is in output
                                  # fields, and if so check for correctness
for k in output_keys:
  if (k[:14] == 'original_input'):
    v = k[14:].strip()
    if (v != ''):  # There is a field or column range given
      if (v[0] == '[') and (v[-1] == ']'):
        v = v[1:-1]  # Remove brackets
      else:
        inout.log_message('Wrong input component definition: '+str(k) + \
                          ' for "original_input" output field','err')
        raise Exception()
      if (v[0] == '(') and (v[-1] == ')'):  # It's a tuple
        v = v[1:-1]  # Remove tuple brackets
      v = v.split(',')  # Make a list
      for i in range(len(v)):
        v[i] = int(v[i])  # Make integers
      if (len(v) == 1):  # One integer only, must be a field number
        input_values.append(v)  # Append 'original_input' field number
      elif (len(v) == 2):  # Two integers, must be a column range
        input_values.append([(v[0],v[1])])  # Append as a tuple
      else:
        inout.log_message('Wrong input component value: '+str(k) + \
                          ' for "original_input" output field','err')
        raise Exception()

for v in input_values:
  if (v != []):
    for e in v:  # For each element in this list
      if (in_file_type == 'COL') and (type(e) == types.TupleType):
        if (len(e) == 2):
          if ((int(e[0])+int(e[1])) > input_len):
            input_len = int(e[0])+int(e[1])
        else:
          inout.log_message(['Wrong input component definition: '+str(v), \
                 'for COL input file type (wrong tuple size)'],'err')
          raise Exception()
      elif (in_file_type == 'COL') and (type(e) != types.TupleType):
        inout.log_message(['Illegal input component definition: '+str(v), \
                 'for COL input file type (elements must be tuples)'],'err')
        raise Exception()
      elif (in_file_type in ['CSV','CSVQ','TAB','TABQ']) and \
         (type(e) == types.IntType):
        if (e > input_len):
          input_len = e
      elif (in_file_type in ['CSV','CSVQ','TAB','TABQ']) and \
         (type(e) != types.IntType):
        inout.log_message(['Illegal input component definition: '+str(v), \
         'for CSV or TAB input file type (elements must be integers)'],'err')
        raise Exception()

if (in_file_type in ['CSV','CSVQ','TAB','TABQ']):
  input_len += 1  # Input field numbering starts with zero (for fields)

# -----------------------------------------------------------------------------
# Check if definition of output components is correct with file types
#
output_values = output_field.values()

for v in output_values:
  if (v != []):
    for e in v:  # For each element in this list
      if (out_file_type == 'COL') and (type(e) != types.TupleType):
        inout.log_message(['Illegal output component definition: '+str(v), \
                 'for COL output file type (elements must be tuples)'],'err')
        raise Exception()
      elif (out_file_type in ['CSV','CSVQ','TAB','TABQ']) and \
         (type(e) != types.IntType):
        inout.log_message(['Illegal output component definition: '+str(v), \
        'for CSV or TAB output file type (elements must be integers)'],'err')
        raise Exception()

# -----------------------------------------------------------------------------
# Set up a CSV parser (with field separator according to input file type)
#
if (project.in_file_type in ['TAB','TABQ']):
  input_line_sep = '\t'
elif (project.in_file_type in ['CSV','CSVQ']):
  input_line_sep = ','

# Use Tim Churches flexible CSV parser
#
if (project.in_file_type in ['CSV','CSVQ','TAB','TABQ']):
  line_parser = tcsv.delimited_parser(delimiter_chars=input_line_sep, \
                                      as_strings=1)

# -----------------------------------------------------------------------------
# Create output separators and field quotes
#
if (out_file_type in ['CSVQ','TABQ']):
  out_file_quoted = 1
else:
  out_file_quoted = 0

if (out_file_type in ['CSV','CSVQ']):
  out_field_sep = ','  # Fields separator is a comma
elif (out_file_type in ['TAB','TABQ']):
  out_field_sep = '\t'  # Field separator is a tabulator

# =============================================================================
# Load correction-list files into lists & lookup-table files into dictionaries

name_corr_list   = inout.load_corr_list(project.name_corr_list_file)
geoloc_corr_list = inout.load_corr_list(project.geoloc_corr_list_file)

[name_lookup_dict, name_dict_seq_len] = \
  inout.load_lookup_tables(project.name_lookup_table_files)
[geoloc_lookup_dict, geoloc_dict_seq_len] = \
  inout.load_lookup_tables(project.geoloc_lookup_table_files)

inout.log_message(['Loaded correction-lists and lookup-tables:', \
                   '  name_corr_list has     '+str(len(name_corr_list))+ \
                   ' entries', \
                   '  geoloc_corr_list has   '+str(len(geoloc_corr_list))+ \
                   ' entries', \
                   '  name_lookup_dict has   '+str(len(name_lookup_dict))+ \
                   ' entries',\
                   '  name_lookup_dict maximal sequence length is '+ \
                   str(name_dict_seq_len), \
                   '  geoloc_lookup_dict has '+str(len(geoloc_lookup_dict))+ \
                   ' entries', \
                   '  geoloc_lookup_dict maximal sequence length is '+ \
                   str(geoloc_dict_seq_len)],'v1')

#=============================================================================
# Create a sorted list of the output fields (according to column numbers)

output_field_names     = []
output_field_positions = []

for k in output_field.keys():  # Only extract fields that are not empty
  if (output_field[k] != []):
    output_field_names.append(k)
    output_field_positions.append(output_field[k])

output_field_list = map(None, output_field_positions, output_field_names)
output_field_list.sort()

# Make sure fields are numbered in sequence for CVS and TAB output file types -
#
if (project.out_file_type in ['CSV','CSVQ','TAB','TABQ']):
  i = 0
  for f in output_field_list:  # Loop over output fields
    if (int(f[0][0]) != i):
      inout.log_message('Illegale output field numbering (not in sequence)' + \
                        ' at position '+ str(i)+': '+ str(output_field_list), \
                        'err')
      raise Exception()
    i += 1

elif (project.out_file_type in ['COL']):  # And check columns for column files
  i = 0
  for f in output_field_list:  # Loop over output fields
    s_col = int(f[0][0][0])   # Start column
    length = int(f[0][0][1])  # Length (number of characters)
    if (s_col != i):
      inout.log_message('Illegale output column sequence (not continuous)' + \
                        ' at columns '+ str(s_col)+' with length '+ \
                        str(lenght)+': ' + str(output_field_list), 'err')
      raise Exception()
    i += length  # Start column of next field

# =============================================================================
