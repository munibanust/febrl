# =============================================================================
# project.py - Configuration and settings related to a data standardisation
#              and record linkage project
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
# The Original Software is "project.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module project.py - Configuration and settings related to a data
                       standardisation and record linkage project

   This module contains all user adjustable parameters andsettings that are
   related to a particular project:
   - File names and file types
   - Definitions of input and output fields
   - Hidden Markov Model files used
   - File names of lookup- and frequency table files
   - Parameters for the data cleaning, standardisation and linkage processes

   This is the ONLY file a user has to edit before she/he can start a data
   standardisation and linkage process.
"""

# =============================================================================
# Verbose output and file logging options.
# - Possible values for 'verbose' are:
#     0  - No verbose output
#     1  - Moderate verbose output
#     2  - Extensive verbose output
# - Possible values for 'logging' are:
#     0  - No logging to file
#     1  - Logging to file activated
# - 'log_file' is the name of the file logging information will be append to.
#   If this file doesn't exist it will be created.
# - Possible values for 'nowarn' are:
#     0  - Do print warning messages
#     1  - Do not print warning messages
#
# NOTE:
# - The 'nowarn' value determines only if warning messages are printed or not.
# - Logging is done according to the value of the verbose output.
# - Warning and error messages are always logged if logging is activated.
# - Both verbose and logging settings can be overwritten with command line
#   arguments (-l, -v1, -v2)

verbose  = 0
logging  = 0
log_file = './logging.txt'
nowarn   = 0

# Number of records between a process indication message is printed
# Set to -1 to disable process indication messages
#
proc_ind = 100

# =============================================================================
# File names and file types.
# in_file_name   - Name of the input file containing the original data
# in_file_type   - The type of the input file (see below for a list of possible
#                  file types)
# out_file_name  - Name of the output file containing the cleaned and
#                  standardised data
# out_file_type  - The type of the output file (see below for a list of
#                  possible file types)
#
# Possible file types are:
#   CSV  - Comma separated values, fields separated by commas
#   CSVQ - Comma separated values, where each field starts and ends with
#          a quote character
#   TAB  - Tabulator separated values, fields separated by commas
#   TABQ - Tabulator separated values, where each field starts and ends with
#          a quote character
#   COL  - Column wise, fields within specific column ranges
#
#   Note: - CSVQ and TABQ are only useful for the output file type, the input
#           parsing routine automatically handles both quoted and non-quoted
#           files.
#         - A database access file type (SQL) is planned to be included in a
#           future release of this software.

in_file_name  = '../../data/nswhealth_mdc/mdc.csv'
#in_file_name = './nametest10000.csv'
in_file_type  = 'CSV'
out_file_name = './test.csv'
out_file_type = 'CSV'

# =============================================================================
# Definition and settings for the input components (name, geocode, locality and
# dates).
# - The 'input_component' dictionary contains information on which field(s) or
#   which columns from an input record should be assigned to which component.
#   With field oriented input files (comma or tabulator separated fields) each
#   component is specified as a list of field numbers (starting with 0).
#   For column oriented input files, a list of column ranges (start,length)
#   needs to be specified for each component.
# - The 'input_space_sep' dictionary contains a flag (set to 0 or 1) for each
#   component. If set, one space character is inserted between the fields of
#   a component before they are concatenated.
# - Only one name, geocode and locality component is supported, but two date
#   components (date1 and date2) are possible.
# - If a component is not available or not used set the corresponding entry
#   to an empty list.
# - For the name component, if given names and surnames are already available
#   two separated input fuelds, it is possible to define the 'givenname' and
#   'surname' components (and then use name.parse_name_component' routines).
#   If you want seperate givenname and surname components, set the 'name'
#   in 'input_component' to an empty list []
#
# Examples:
# - For file what contains comma or tabulator seperated fields:
#   input_component = {'name':[4,5],    # Name component in fields 4 and 5
#                      'geocode':[8,9]} # Geocode component in fields 8 and 9
# - For file what contains column oriented fields:
#   input_component = {'name':[(0,20),(22,9)]}
#              # Name component is in columns 0 to 19 and 22 to 31

# dth8502.csv: name: 5,6,7; locality: 2,3,4; geocode: 1
#              gname: 5, sname: 6,7
# mdc: name: 4,5; locality: 10,11,12; geocode: 8,9, dates: 6,14
#      gname: 4, sname: 5

input_component = {'name':      [4,5], \
                   'givenname': [], \
                   'surname':   [], \
                   'geocode':   [8,9], \
                   'locality':  [10,11,12], \
                   'date1':     [6], \
                   'date2':     [14]}

# -----------------------------------------------------------------------------
# If a space separator should be added between the fields, set the according
# value for a component to 1, otherwise set it to 0.

input_space_sep = {'name':      1, \
                   'givenname': 1, \
                   'surname':   1, \
                   'geocode':   1, \
                   'locality':  1, \
                   'date1':     0, \
                   'date2':     0}

# -----------------------------------------------------------------------------
# If word 'spilling' between fields should be checked set the value for a
# component to 1 (e.g. "miller r","d canberra"), set to 0 if no checking
# should be performed

input_check_spilling = {'name':      1, \
                        'givenname': 1, \
                        'surname':   1, \
                        'geocode':   1, \
                        'locality':  1, \
                        'date1':     0, \
                        'date2':     0}

# =============================================================================
# Definition and settings for the output fields to be written into the output
# file.
# - Similar to the definition of the input components, for each output field a
#   field number (for comma and tabulator separated files) or a start column
#   and length of the field (for column wise files) can be given.
# - An empty list [] for an output fields means it will not be written into
#   the output file. Column and field numbers start with 0.
# - Two special output fields are 'name_hmm_proba' and 'geoloc_hmm_proba' which
#   are the probabilities returned by the Viterbi algorithm for the most likely
#   HMM state sequence (which is the one chosen for the standardisation)
# - The record number for each processed record can be added to the output
#   by using the special field 'record_id'. Record numbers are counted from
#   the beginning of the input file, starting with zero. Each line in the
#   input file is counted as one record. 
# - It is possible to write parts or all of the original input record
#   unmodified into the output file by using the special output field
#   'original_input'. Three forms of 'original_input' are possible:
#   - 'original_input[in_field_num]' for comma/tabulator separated input files
#   - 'original_input[in_start_col,in_end_col]' for column wise input files
#   - 'original_input' for the whole input record (the original input line)
#   Examples:
#   'original_input[1]':[2]   writes the second (field numbering is starting
#                             input field into the third output field.
#   'original_input':[0]      writes the complete original input line into the
#                             first output field.

output_field = {'record_id':          [0], \
                'original_input[1]':  [], \
                'original_input[4]':  [], \
                'original_input':     [], \
                'name_hmm_proba':     [], \
                'geoloc_hmm_proba':   [], \
                'title':              [], \
                'gender_guess':       [1], \
                'givenname':          [2], \
                'alt_givenname':      [3], \
                'surname':            [4], \
                'alt_surname':        [5], \
                'wayfare_number':     [6], \
                'wayfare_name':       [7], \
                'wayfare_type':       [8], \
                'wayfare_qualifier':  [9], \
                'unit_type':          [10], \
                'unit_number':        [11], \
                'property_name':      [12], \
                'institution_name':   [13], \
                'institution_type':   [14], \
                'postaddress_number': [15], \
                'postaddress_type':   [16], \
                'postcode':           [17], \
                'locality_name':      [18], \
                'locality_qualifier': [19], \
                'territory':          [20], \
                'country':            [], \
                'day1':               [21], \
                'month1':             [22], \
                'year1':              [23], \
                'day2':               [24], \
                'month2':             [25], \
                'year2':              [26]}

# -----------------------------------------------------------------------------
# For quoted data files, the quote character has to be given below.
# For unquoted output files, if a comma is embedded in an output field, such
# a field will be quoted with the output quote character given below.

output_quote_character = '"'

# =============================================================================
# File names of correction-list and lookup-table files

name_corr_list_file     = './data/name_corr.lst'
geoloc_corr_list_file   = './data/geoloc_corr.lst'

name_lookup_table_files   = ['./data/givenname_f.tbl', \
                             './data/givenname_m.tbl', \
                             './data/name_prefix.tbl', \
                             './data/name_misc.tbl', \
                             './data/saints.tbl', \
                             './data/surname.tbl', \
                             './data/title.tbl']

geoloc_lookup_table_files = ['./data/country.tbl', \
                             './data/geoloc_misc.tbl', \
                             './data/geoloc_qual.tbl', \
                             './data/institution_type.tbl', \
                             './data/locality_name_act.tbl', \
                             './data/locality_name_nsw.tbl', \
                             './data/post_address.tbl', \
                             './data/postcode_act.tbl', \
                             './data/postcode_nsw.tbl', \
                             './data/saints.tbl', \
                             './data/territory.tbl', \
                             './data/unit_type.tbl', \
                             './data/wayfare_type.tbl']

# =============================================================================
# Options for name component
#
# -----------------------------------------------------------------------------
# Definition of which method to use: 'hmm' or 'rules'
# Currently, only 'hmm' is supported for the geocode/locality components

name_standard_method = 'hmm'
name_hmm_file_name = './hmm/name-absdiscount.hmm'

name_female_title = ['ms']
name_male_title   = ['mr']

# =============================================================================
#
# Options for locality component
#
# -----------------------------------------------------------------------------
# Definition of which method to use: 'hmm' or 'rules'
# (currently, only 'hmm' is supported, rule based routines not yet implemented)

geoloc_standard_method = 'hmm'
geoloc_hmm_file_name = './hmm/geoloc-absdiscount.hmm'

# =============================================================================
# Options for date component
#
# -----------------------------------------------------------------------------
# Definition of the date parsing format strings
#
# Date parsing and standardisation is done using a format string (similar to
# the Python time.strptime() function).
# - A format string, must be made of three directives, which can either be
#   written directly one after the other or being separated by a space between
#   them (e.g. "%d%b%Y" or "%m %d %y")
# - The list of given formats is ordered, and the first one that can parse
#   a given date string correctly will be applied.
# - The following format directives are possible:
#     %b  Abbreviated month name (Jan, Feb, Mar, etc.)
#     %B  Full month name (January, February, etc.)
#     %d  Day of the month as a decimal number [01,31]
#     %m  Month as a decimal number [01,12]
#     %y  Year without century as a decimal number [00,99]
#     %Y  Year with century as a decimal number

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

                      #'%Y %B %d',   # 2002 Apr 24    ??? Is this needed?
                      #'%y %B %d',   # 02 Apr 24      ??? Is this needed?
                     ]

# -----------------------------------------------------------------------------
# Pivot year: Two-digits years smaller than the pivot year will be expanded
# into 20XX, years larger and equal than the pivot year will be expanded into
# 19xx
# Example: pivot_year = 03:  68 -> 1968, 03 -> 1903 and 02 -> 2002

date_pivot_year = 03

# -----------------------------------------------------------------------------
# Date percentage comparison relative fix date (how do we call this??)
#
date_perc_fix_date = 'today'  # Set to a date or to string 'today'

# Date age computation relative fix date
#
date_age_fix_date = 'today'  # Set to a date or to string 'today'

# -----------------------------------------------------------------------------
# Date matching and unmatching probabilities
#
# As an example, we set the matching (m) and unmatching (u) probabilities
# similar to the example given in L. Gill, 2002, pages 64/65.
# Days:  m=0.95  (assume 5% errors)
#        u=1/30  (assume average month length 30 days)
# Month: m=0.95  (assume 5% errors)
#        u=1/12  (random agreement for month in 1/12 of all pairs)
# Years: m=0.95  (assume 5% errors)
#        u=1/100 (assume simplified uniform age distribution over 100 years)

# Probability that days are the same in a linked pair
date_day_m_prob = 0.95

# Probability that days are the same in an unlinked pair
date_day_u_prob = 0.0328

# Probability that months are the same in a linked pair
date_month_m_prob = 0.95

# Probability that months are the same in an unlinked pair
date_month_u_prob = 0.0833

# Probability that years are the same in a linked pair
date_year_m_prob = 0.95

# Probability that years are the same in an unlinked pair
date_year_u_prob = 0.01

# -----------------------------------------------------------------------------
# Date comparison tolerations for weight computation
# Set any of these variables to an empty string, i.e. '', if they shouldn't be
# considered for the final date weight computation.
# 
date_comp_max_subst = [1,1,2]  # Maximum number of substitutions
                               # tolerated (for day,month, year)
date_comp_max_trans = [1,0,1] # Maximum number of transpositions
                              # tolerated (for day,month, year)
date_comp_max_day_before = ''
date_comp_max_day_after =  ''

#date_comp_max_day_before = 10  # Number of days tolerated for the first date
                               # being before the second date
#date_comp_max_day_after =  10  # Number of days tolerated for the first date
                               # being after the second date

date_comp_max_perc_before = 10.0  # Numerical value in percent for the first
                                  # date being before the second date
date_comp_max_perc_after = 20.0   # Numerical value in percent for the first
                                  # date being after the second date

# -----------------------------------------------------------------------------
# Final date linkage weight computation, which can be the minimum, maximum or
# a combination of the four weights: substitution_weight, transposition_weight,
# day_diff_weight and perc_diff_weight.
#
#date_linkage_weight_comp = 'min'
#date_linkage_weight_comp = 'max'
date_linkage_weight_comp = [0.1,0.2,0.3,0.4]  # Fractional weights, sum = 1.0 
  # [substitution_weight,transposition_weight,day_diff_weight,perc_diff_weight]

# =============================================================================
# END
# =============================================================================
