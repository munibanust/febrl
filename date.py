# =============================================================================
# date.py - Routines for date parsing, comparing and linkage.
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
# The Original Software is "date.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module date.py - Routines for date parsing, comparing and linkage.

   PUBLIC FUNCTIONS:
     parse_datestr        Parse a date string and split into [day,month,year]
     date_diff            Compare two dates and return their difference in
                          days and percentage
     date_comp            Compare two dates, and return number of
                          transpositions and substitutions
     date_linkage_weight  Compute the linkage weight for a pair of dates
     epoch_to_date        Convert a Unix epoch day number into a date
     date_to_epoch        Convert a date into a Unix epoch day integer
     date_to_age          Convert a date into an age (relative to a fix date)
     str2date             A routine that converts a string into a date using
                          a format string
     test                 Test routine which does basic tests for all functions

   See doc strings of individual functions for detailed documentation.

   See also the relevant section in the config.py module.
"""

# -----------------------------------------------------------------------------

import types
import time
import string

import config
import inout
import mymath

# -----------------------------------------------------------------------------

def parse_datestr(datestr):
  """Normalise a date string and split into year, month and day.

  USAGE:
    [year,month,day,status] = parse_datestr(datestr)

  ARGUMENTS:
    datestr  Input date as a string. Various formats possible:
                'DD/MM/YYYY', 'DD/MM/YY', 'DD,MM,YYYY', 'DD,MM,YY',
                'D Mon YY', etc.
             See config.py for a list of supported date formats.

  DESCRIPTION:
    The function tries to parse date information from the input string.
    If not successful, an error status is returned.

    The formats are given in the config.py module, where a list of formats
    (date_formats) can be defined by the user.

    The return value of status is 'OK' (if the date string has been parsed
    successfully) or 'ER' (if it wasn't possible to parse the date string).
    In the latter case, the values for day, mont hand year are all set to -1.

  EXAMPLE:
    datestr = ' 1 January 2001  '
    [day,month,year,status] = parse_datestr(datestr)
  """

  if (type(datestr) != types.StringType):
    inout.log_message('Input argument is not of string type: '+str(datestr),\
                      'err')
    raise TypeError(parse_datestr)

  # Remove leading and trailing whitespaces - - - - - - - - - - - - - - - - - -
  #
  datestr = datestr.strip()

  # Apply replace table - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  string_replace = ["'.,:-=_/\\", \
                    "         "]
  # Characters in the first list are replaced by
  # the corresponding character in the second list

  replace_table = string.maketrans(string_replace[0], string_replace[1])
  tmp_str = datestr.translate(replace_table)

  # Number of defined date formats
  num_date_formats = len(config.date_parse_formats)

  parsed = 0  # Flag for successful date parsing
  i = 0  # Interation counter
  
  # Try one date format after the other until success or none worked  - - - - -
  #
  while (parsed != 1) and (i < num_date_formats):

    date_try = str2date(tmp_str,config.date_parse_formats[i])
    if (len(date_try) == 3):
      parsed = 1  # Parsed successful
    else:
      inout.log_message('Date format string "'+(config.date_parse_formats[i])+\
                        '" did not match date: '+tmp_str,'v2')
    i=i+1

  if (parsed == 1):  # Date successfully parsed - - - - - - - - - - - - - - - -
    day =   date_try[0]
    month = date_try[1]
    year =  date_try[2]

    if (year < 100): # Year was most probably only two digits, expand into four
      if (year < config.date_pivot_year):  # Map into four digit year
        year = year+2000
      else:  # For years between 04 and 99
        year = year+1900

    inout.log_message('  Date parsed: Day='+str(day)+', month='+str(month)+ \
                      ', year='+str(year),'v1')
    result = [day,month,year,'OK']

  else:
    inout.log_message(' Input argument is not in a valid date format: '+ \
                      str(datestr),'warn')
    result = [-1,-1,-1,'ER']

  return result

# -----------------------------------------------------------------------------

def date_diff(date1, date2):
  """Compate two dates and return their difference in days and percentage.

  USAGE:
    [day_diff, perc_diff] = date_diff(date1, date2)

  ARGUMENTS:
    date1  Date tuple/list in format [day,month,year]
    date2  Date tuple/list in format [day,month,year]
 
  DESCRIPTION:
    Returns the absolute day difference between the two input dates, as well as
    a percentage difference relative to a date that is defined in config.py
    (if this date is not set, the current system date ('today') is taken).

    If the absolute day difference is negative, the first date is before the
    second date, if the day difference is positive then the first date is
    after the second date.
    The percentage difference returned is always positive.
  """

  # Get fix date from config.py or make fix date 'today'  - - - - - - - - - - -
  #
  if (not config.date_perc_fix_date) or (config.date_perc_fix_date == 'today'):
    sys_time = time.localtime(time.time())  # Get current system date
    fix_date = [sys_time[2],sys_time[1],sys_time[0]]
  else:
    parsed_date = parse_datestr(config.date_perc_fix_date) # config.py fix date
    fix_date = parsed_date[0:3]

  inout.log_message('Date percentage comparison fix date: '+str(fix_date),'v2')

  # Get epoch number for both input dates - - - - - - - - - - - - - - - - - - -
  #
  date1_epoch = date_to_epoch(date1[0],date1[1],date1[2])
  date2_epoch = date_to_epoch(date2[0],date2[1],date2[2])

  day_diff = date1_epoch - date2_epoch  # Get absolute day difference
  
  if (day_diff == 0):  # Check if the dates are equal
    perc_diff = 0  # No percentage difference

  else:  # The two dates are not equal - - - - - - - - - - - - - - - - - - - -
    fix_epoch = date_to_epoch(fix_date[0],fix_date[1],fix_date[2])

    date1_abs = date1_epoch - fix_epoch
    date2_abs = date2_epoch - fix_epoch

    if (date1_abs == 0) or (date2_abs == 0):
      diff = 100.0  # Set percentage to a very high level
    elif ((date1_abs < 0) and (date2_abs >= 0)) or \
       ((date2_abs < 0) and (date1_abs >= 0)):
      date1_abs = abs(date1_abs)
      date2_abs = abs(date2_abs)
      diff = float(date1_abs + date2_abs) / float(min(date1_abs, date2_abs))
    else:
      diff = float(abs(date1_abs - date2_abs)) /  \
             float(min(abs(date1_abs), abs(date2_abs)))

    perc_diff = diff*100.0

  return [day_diff, perc_diff]

# -----------------------------------------------------------------------------

def date_comp(date1,date2): 
  """Compare two dates, and return number of transpositions and substitutions.

  USAGE:
    [trans, subst] = date_comp(date1, date2)

  ARGUMENTS:
    date1  Date tuple/list in format [day,month,year]
    date2  Date tuple/list in format [day,month,year]
 
  DESCRIPTION:
    Computes and returns the number of characters that are transposed and
    substituted between the two dates.

    Transpositions are only considered between neighbouring digits.
  """

  # Make strings out of input dates - - - - - - - - - - - - - - - - - - - - - -
  #
  day1_str = string.zfill(str(date1[0]),2)
  day2_str = string.zfill(str(date2[0]),2)
  month1_str = string.zfill(str(date1[1]),2)
  month2_str = string.zfill(str(date2[1]),2)
  year1_str = str(date1[2])
  year2_str = str(date2[2])

  trans = [0,0,0]  # Number of transpositions for days, months and years
  subst = [0,0,0]  # Number of substitutions for days, months and years

  # Check in days - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (day1_str[0] != day1_str[1]) and \
     (day1_str[0] == day2_str[1]) and (day1_str[1] == day2_str[0]):
    trans[0] = trans[0]+1
  else:  # No transposition in days, check for substitutions
    if (day1_str[0] != day2_str[0]):
      subst[0] = subst[0]+1
    if (day1_str[1] != day2_str[1]):
      subst[0] = subst[0]+1

  # Check in months - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (month1_str[0] != month1_str[1]) and \
     (month1_str[0] == month2_str[1]) and (month1_str[1] == month2_str[0]):
    trans[1] = trans[1]+1
  else:  # No transposition in months, check for substitutions
    if (month1_str[0] != month2_str[0]):
      subst[1] = subst[1]+1
    if (month1_str[1] != month2_str[1]):
      subst[1] = subst[1]+1

  # Check in years - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (year1_str[0] != year1_str[1]) and \
     (year1_str[0] == year2_str[1]) and (year1_str[1] == year2_str[0]):
    trans[2] = trans[2]+1  # Transposition between first two digits

    if (year1_str[2] != year1_str[3]) and \
       (year1_str[2] == year2_str[3]) and (year1_str[3] == year2_str[2]):
      trans[2] = trans[2]+1  # Transposition between last two digits
    else:  # No transposition in last two year digits, check for substitutions
      if (year1_str[2] != year2_str[2]):
        subst[2] = subst[2]+1
      if (year1_str[3] != year2_str[3]):
        subst[2] = subst[2]+1

  else:  # No transposition in first two year digits
    if (year1_str[0] != year2_str[0]):  # Check for substitution in first digit
      subst[2] = subst[2]+1

    if (year1_str[1] != year1_str[2]) and \
       (year1_str[1] == year2_str[2]) and (year1_str[2] == year2_str[1]):
      trans[2] = trans[2]+1  # Transposition between middle two digits

      if (year1_str[3] != year2_str[3]):  # Check substitution in last digit 
        subst[2] = subst[2]+1

    else:  # No transposition in middle two digits
      if (year1_str[1] != year2_str[1]):  # Check for substitution in 2. digit
        subst[2] = subst[2]+1

      if (year1_str[2] != year1_str[3]) and \
         (year1_str[2] == year2_str[3]) and (year1_str[3] == year2_str[2]):
        trans[2] = trans[2]+1  # Transposition between last two digits
      else: # No transposition in last two year digits, check for substitutions
        if (year1_str[2] != year2_str[2]):
          subst[2] = subst[2]+1
        if (year1_str[3] != year2_str[3]):
          subst[2] = subst[2]+1

  return [trans, subst]

# -----------------------------------------------------------------------------

def date_linkage_weight(date1,date2):
  """Compute the linkage weight for a pair of dates.

  USAGE:
    [link_weight] = date_linkage_weights(date1, date2)

  ARGUMENTS:
    date1  Date tuple/list in format [day,month,year]
    date2  Date tuple/list in format [day,month,year]
 
  DESCRIPTION:
    Compute the linkage weight according to probablistic linkage theory.
    
    Match (M) and non-match (U) probabilities are given in config.py and can
    be modified by the user. Seperate probabilities are given for day, month
    and year.

    Four different linkage wieghts are computed:
    1) According to the number of substitutions
    2) According to the number of transpositions
    3) According to the absolute day difference
    4) According to the percentage day difference

    The final linkage weight that is returned can ba any combination of these
    four weights according to setting given in config.py

    If one of the dates is not valid (status 'ER' or one of the elements is -1,
    the maximum diagreement weight is returned.
  """

  # Compute initial linking weights for date components
  #
  # TODO: Move this into config or somewhere, as this is 'static' and needs to
  #       be done once only
  #
  day_agree_weight = \
    mymath.log2(config.date_day_m_prob / config.date_day_u_prob)
  month_agree_weight = \
    mymath.log2(config.date_month_m_prob / config.date_month_u_prob)
  year_agree_weight = \
    mymath.log2(config.date_year_m_prob / config.date_year_u_prob)

  day_disagree_weight = \
    mymath.log2((1.0-config.date_day_m_prob) / (1.0-config.date_day_u_prob))
  month_disagree_weight = \
    mymath.log2((1.0-config.date_month_m_prob) /(1.0-config.date_month_u_prob))
  year_disagree_weight  = \
    mymath.log2((1.0-config.date_year_m_prob) / (1.0-config.date_year_u_prob))

  total_agree_weight = day_agree_weight +  month_agree_weight + \
                       year_agree_weight
  total_disagree_weight = day_disagree_weight + month_disagree_weight + \
                          year_disagree_weight
  
  day_weight_range = day_agree_weight + abs(day_disagree_weight)
  month_weight_range = month_agree_weight + abs(month_disagree_weight)
  year_weight_range = year_agree_weight + abs(year_disagree_weight)
  total_weight_range = total_agree_weight + abs(total_disagree_weight)

  if (config.verbose > 0):
    print '  Dates day agreement weight:', day_agree_weight
    print '  Dates day disagreement weight:', day_disagree_weight
    print '  Dates day weight range:', day_weight_range

    print '  Dates month agreement weight:', month_agree_weight
    print '  Dates month disagreement weight:', month_disagree_weight
    print '  Dates month weight range:', month_weight_range

    print '  Dates year agreement weight:', year_agree_weight
    print '  Dates year disagreement weight:', year_disagree_weight
    print '  Dates year weight range:', year_weight_range

    print '  Dates total agreement weight:', total_agree_weight
    print '  Dates total disagreement weight:', total_disagree_weight
    print '  Dates total weight range:', total_weight_range


  # Special case: The dates are equal, no further checking required - - - - - -
  #
  if (date1 == date2):
    return total_agree_weight

  # Special case: A date is not valid (unsuccessful parsing) - - - - - - - - -
  #
  if (-1 in date1) or (-1 in date2) or ('ER' in date1) or ('ER' in date2):
    return total_disagree_weight

  # Get day and percentage differences, substitutions and transpositions
  #
  [day_diff, perc_diff] = date_diff(date1,date2)
  [trans, subst] = date_comp(date1,date2)

  # Check substitutions - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (config.date_comp_max_subst != ''):  # Consider substitutions
    if (subst[0] == 0):  # No substitutions in days
      subst_day_weight = day_agree_weight  # Set to agreement weight

    elif (subst[0] <= config.date_comp_max_subst[0]):
      penalty_fac = float(subst[0]) / \
                    (1.0 + config.date_comp_max_subst[0])
      subst_day_weight = day_agree_weight - penalty_fac * day_weight_range

    else:  # Too many substitutions in day field
      subst_day_weight = day_disagree_weight  # Set to disagreement weight

    if (subst[1] == 0):  # No substitutions in month
      subst_month_weight = month_agree_weight  # Set to agreement weight

    elif (subst[1] <= config.date_comp_max_subst[1]):
      penalty_fac = float(subst[1]) / \
                    (1.0 + config.date_comp_max_subst[1])
      subst_month_weight = month_agree_weight - penalty_fac*month_weight_range

    else:  # Too many substitutions in month field
      subst_month_weight = month_disagree_weight  # Set to disag. weight

    if (subst[2] == 0):  # No substitutions in year
      subst_year_weight = year_agree_weight  # Set to agreement weight

    elif (subst[2] <= config.date_comp_max_subst[2]):
      penalty_fac = float(subst[2]) / \
                    (1.0 + config.date_comp_max_subst[2])
      subst_year_weight = year_agree_weight - penalty_fac * year_weight_range

    else:  # Too many substitutions in year field
      subst_year_weight = year_disagree_weight  # Set to disag. weight

    subst_weight = subst_day_weight + subst_month_weight + subst_year_weight

  else: # Don't consider month substitution
    subst_weight = ''  # Set to a non-numerical value

  # Check transpositions - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (config.date_comp_max_trans != ''):  # Consider transpositions
    if (trans[0] == 0):  # No transpositions in days
      trans_day_weight = day_agree_weight  # Set to agreement weight

    elif (trans[0] <= config.date_comp_max_trans[0]):
      penalty_fac = float(trans[0]) / \
                    (1.0 + config.date_comp_max_trans[0])

      trans_day_weight = day_agree_weight - penalty_fac * day_weight_range

    else:  # Too many transpositions in day field
      trans_day_weight = day_disagree_weight  # Set to disagreement weight

    if (trans[1] == 0):  # No transpositions in month
      trans_month_weight = month_agree_weight  # Set to agreement weight

    elif (trans[1] <= config.date_comp_max_trans[1]):
      penalty_fac = float(trans[1]) / \
                    (1.0 + config.date_comp_max_trans[1])
      trans_month_weight = month_agree_weight - penalty_fac*month_weight_range

    else:  # Too many transpositions in month field
      trans_month_weight = month_disagree_weight  # Set to disag. weight

    if (trans[2] == 0):  # No transpositions in year
      trans_year_weight = year_agree_weight  # Set to agreement weight

    elif (trans[2] <= config.date_comp_max_trans[2]):
      penalty_fac = float(trans[2]) / \
                    (1.0 + config.date_comp_max_trans[2])
      trans_year_weight = year_agree_weight - penalty_fac * year_weight_range

    else:  # Too many transpositions in year field
      trans_year_weight = year_disagree_weight  # Set to disag. weight

    trans_weight = trans_day_weight + trans_month_weight + trans_year_weight

  else: # Don't consider month transposition
    trans_weight = ''  # Set to a non-numerical value

  # Check absolute day differences - - - - - - - - - - - - - - - - - - - - - -
  #
  if (day_diff < 0) and (config.date_comp_max_day_before != ''):
    if (-day_diff <= config.date_comp_max_day_before):
      penalty_fac = float(-day_diff) / \
                    (1.0 + config.date_comp_max_day_before)
      day_diff_weight = total_agree_weight - penalty_fac * total_weight_range

    else:  # Day difference too large
      day_diff_weight = total_disagree_weight  # Set to disagr. weight
   
  elif (day_diff > 0) and (config.date_comp_max_day_after != ''):
    if (day_diff <= config.date_comp_max_day_after):
      penalty_fac = float(day_diff) / \
                    (1.0 + config.date_comp_max_day_after)
      day_diff_weight = total_agree_weight - penalty_fac * total_weight_range

    else:  # Day difference too large
      day_diff_weight = total_disagree_weight  # Set to disagr. weight

  else:
    day_diff_weight = ''  # Set to a non-numerical value

  # Check percentage day differences - - - - - - - - - - - - - - - - - - - - -
  #
  if (day_diff < 0) and (config.date_comp_max_perc_before != ''):
    if (perc_diff <= config.date_comp_max_perc_before):
      penalty_fac = perc_diff / (1.0 + config.date_comp_max_perc_before)
      perc_diff_weight = total_agree_weight - penalty_fac * total_weight_range

    else:  # Percentage difference too large
      perc_diff_weight = total_disagree_weight  # Set to disagr. weight
   
  elif (day_diff > 0) and (config.date_comp_max_perc_after != ''):
    if (perc_diff <= config.date_comp_max_perc_after):
      penalty_fac = perc_diff / (1.0 + config.date_comp_max_perc_after)
      perc_diff_weight = total_agree_weight - penalty_fac * total_weight_range

    else:  # Percentage difference too large
      perc_diff_weight = total_disagree_weight  # Set to disagr. weight

  else:
    perc_diff_weight = ''  # Set to a non-numerical value

  msg = ['  Substitution weight:         '+str(subst_weight), \
         '  Transposition weight:        '+str(trans_weight), \
         '  Day difference weight:       '+str(day_diff_weight), \
         '  Percentage difference weight:'+str(perc_diff_weight)]
  inout.log_message(msg,'v2')

  # Create list with all numerical weight values  - - - - - - - - - - - - - - -
  #
  final_weights = []
  if (subst_weight != ''):
    final_weights.append(subst_weight)
  if (trans_weight != ''):
    final_weights.append(trans_weight)
  if (day_diff_weight != ''):
    final_weights.append(day_diff_weight)
  if (perc_diff_weight != ''):
    final_weights.append(perc_diff_weight)

  # Compute final date weight according to setting from config.py - - - - - - -
  #
  if (config.date_linkage_weight_comp == 'min'):
    date_weight = min(final_weights)

  elif (config.date_linkage_weight_comp == 'max'):
    date_weight = max(final_weights)

  elif (type(config.date_linkage_weight_comp) == types.ListType) and \
       (len(config.date_linkage_weight_comp) == 4):

    # Compute weighted final date weight
    #
    date_weight = 0.0
    if (subst_weight != ''):
      date_weight += subst_weight * config.date_linkage_weight_comp[0]
    if (trans_weight != ''):
      date_weight += trans_weight * config.date_linkage_weight_comp[1]
    if (day_diff_weight != ''):
      date_weight += day_diff_weight * config.date_linkage_weight_comp[2]
    if (perc_diff_weight != ''):
      date_weight += perc_diff_weight * config.date_linkage_weight_comp[3]

  else:
    inout.log_message('Illegal value or type for '+
                      '"config.date_linkage_weight_comp": '+ \
                      config.date_linkage_weight_comp,'err')
    raise TypeError(date_linkage_weights)

  return date_weight

# -----------------------------------------------------------------------------

def epoch_to_date(daynum):
  """Convert a Unix epoch day number into a date [year, month, day].

  USAGE:
    [year, month, day] = epoch_to_date(daynum)

  ARGUMENTS:
    daynum  A integer giving the Unix epoch day (0 = 1970-01-01)

  DESCRIPTION:
    Function for converting a number of days since Unix epoch time (integer
    value) into a date tuple [day, month, year].

  EXAMPLES:
    [year, month, day] = epoch_to_date(0)       # 1970-01-01
    [year, month, day] = epoch_to_date(11736)   # 2002-02-18
  """

  import time
  import types

  if (type(daynum) != types.IntType) and (type(daynum) != types.LongType):
    inout.log_message('Input argument is not of integer type: '+str(daynum), \
                      'err')
    raise TypeError(epoch_to_date)

  datetuple  = time.gmtime(daynum*24*3600)

  return [datetuple[0], datetuple[1], datetuple[2]]

# -----------------------------------------------------------------------------

def date_to_epoch(day, month, year):
  """ Convert a date [day, month, year] into a Unix epoch day number.

  USAGE:
    daynum = date_to_epoch(year, month, day)

  ARGUMENTS:
    day    Day value (integer number)
    month  Month value (integer number)
    year   Year value (integer number)

  DESCRIPTION:
    Function for converting a date into a Unix epoch day number
    (integer value).

    Based on a Perl script.. source unknown

  EXAMPLES:
    day1 = date_to_epoch(18, 02, 2002)  # 11736
    day2 = date_to_epoch(01, 01, 1970)  # 0
  """

  import types

  if (type(year) != types.IntType) and (type(year) != types.LongType) and \
     (type(month) != types.IntType) and (type(month) != types.LongType) and \
     (type(day) != types.IntType) and (type(day) != types.LongType):
    inout.log_message('Input argument is of illegal type', 'err')
    raise TypeError(date_to_epoch)

  # What is all of this about?
  if (month < 3):
    year = year - 1
  if (month > 2):
    month = month - 3
  else:
    month = month + 9

  c = year / 100.0
  ya = year - ( 100 * c )

  return int(((146097*c)/4) + ((1461*ya)/4) + \
             (((153*month)+2)/5) + day - 719469)

# -----------------------------------------------------------------------------

def date_to_age(day, month, year):
  """Convert a date into an age (relative to a fix date given in config.py)

  USAGE:
    age = date_to_age(day, month, year)

  ARGUMENTS:
    day    Day value (integer number)
    month  Month value (integer number)
    year   Year value (integer number)

  DESCRIPTION:
    Returns the age in years as a positive floating-point number, or a negative
    error number. 

    The fix date is defined in config.py as: 'date_age_fix_date'
    (either set to a specific date or to 'today').
  """

  # Get fix date from config.py or make fix date 'today'  - - - - - - - - - - -
  #
  if (not config.date_age_fix_date) or (config.date_age_fix_date == 'today'):
    sys_time = time.localtime(time.time())  # Get current system date
    fix_date = [sys_time[2], sys_time[1], sys_time[0]]
  else:
    parsed_date = parse_datestr(config.date_age_fix_date) # config.py fix date
    fix_date = parsed_date[0:3]

  inout.log_message('Date age computation fix date: '+str(fix_date),'v2')

  # Get epoch number for input date and fix date  - - - - - - - - - - - - - - -
  #
  date_epoch = date_to_epoch(day, month, year)
  fix_epoch  = date_to_epoch(fix_date[0], fix_date[1], fix_date[2])

  day_diff = fix_epoch - date_epoch  # Get day difference

  # Check for valid day difference and compute age  - - - - - - - - - - - - - -
  #
  if (day_diff < 0):
    inout.log_message(['Age computation not possible:', \
                       '  Fix date:  '+str(fix_date), \
                       '  This date: '+str([day,month,year])],'warn')
    age = -1.0  # Given date is before fix date -> Age computation not possible

  else:  # Approximation of age - - - - - - - - - - - - - - - - - - - - - - - -
    age = float(day_diff) / 365.25

  return age

# -----------------------------------------------------------------------------

def str2date(datestr, formatstr):
  """A routine that converts a string into a date using a format string.

  USAGE:
    [year,month,day] = str2date(datestr, formatstr)

  ARGUMENTS:
    datestr    Input date as a string
    formatstr  A format string, must be made of three directives, which can
               either be written one after the other or being separated by
               a space between them (e.g. "%d%b%Y" or "%m %d %y")
               Possible format directives are (similar to Python 'strptime'
               format directives):
                 %b  Abbreviated month name (Jan, Feb, Mar, etc.)
                 %B  Full month name (January, February, etc.)
                 %d  Day of the month as a decimal number [01,31]
                 %m  Month as a decimal number [01,12]
                 %y  Year without century as a decimal number [00,99]
                 %Y  Year with century as a decimal number

  DESCRIPTION:
    This routine parses the input date string according to the given format
    and extracts a [day,month,year] triplet if possible.

    If the routine can't parse the date string successfully an empty
    list is returned.

    Valid four digit year values must be in the interval [1850,2005].
  """

  formatstr = formatstr.strip()
  datestr   = datestr.strip().lower()

  if (datestr =='') or (formatstr ==''):
    return []  # No date or no format string given

  elif (' ' in formatstr) and (' ' in datestr):
    date_list   = datestr.split()
    format_list = formatstr.split()

  elif (' ' not in formatstr) and (' ' not in datestr):
    date_list   = []
    format_list = []

    workformatstr = formatstr
    workdatestr   = datestr
    while (workformatstr != ''):
      if (workformatstr[:2] == '%Y'):  # Four digit year
        format_list.append(workformatstr[:2])
        workformatstr = workformatstr[2:]
        date_list.append(workdatestr[:4])
        workdatestr = workdatestr[4:]

      elif (workformatstr[:2] in ['%m','%d','%y']):  # 2-digit year/month/day 
        format_list.append(workformatstr[:2])
        workformatstr = workformatstr[2:]
        date_list.append(workdatestr[:2])
        workdatestr = workdatestr[2:]

      else:
        inout.log_message('Illegal format string without spaces:'+formatstr, \
                          'err')

  else:  # A space in either date or format string but not in both
    return []  # Illegal combination of format string and date string

  day   = -1  # Set initially to illegal values
  month = -1
  year  = -1

  if (len(format_list) != 3) or (len(date_list) != 3):
    return []  # No valid format or date

  while (format_list != []):  # process all elements in the format list - - - -

    directive = format_list[0]

    if (directive in ['%b','%B']):  # Month abbreviation or name
      if (len(date_list[0]) >= 3):
        # Get month number (between 1..12) or -1 if not a valid month
        month = config.month_abbrev_dict.get(date_list[0][:3],-1)

    elif (directive == '%d'):  # Day of the month number
      try:
        day = int(date_list[0])  # Convert into an integer number
      except:
        day = -1  # Day is no integer

    elif (directive == '%m'):  # Month number (between 1..12)
      try:
        month = int(date_list[0])  # Convert into an integer number
      except:
        month = -1  # Month is no integer

    elif (directive == '%y'):  # Two digit year without century
      if (len(date_list[0]) == 2):  # Parse only if length is 2 digits
        try:
          year = int(date_list[0])  # Convert into an integer number
          if (year < 0) or (year > 99):
            year = -1  # Illegal year value in two digit
        except:
          year = -1  # year is no integer

        # Convert year into four digit value according to user defined pivot
        # year from 'project.py' module (wrapped into 'config.py')
        #
        if (year >= 0) and (year < config.date_pivot_year):
          year = year+2000
        elif (year >= config.date_pivot_year) and (year < 100):
          year = year+1900

    elif (directive == '%Y'):  # Four digit year with century
      if (len(date_list[0]) == 4):  # Parse only if length is 4 digits
        try:
          year = int(date_list[0])  # Convert into an integer number
          if (year < 1850) or (year > 2005):
            year = -1  # Illegal year value in four digit
        except:
          year = -1  # year is no integer

    else:
      inout.log_message('Illegal format directive: '+directive, 'err')

    date_list = date_list[1:]  # Remove processed first element
    format_list = format_list[1:]

  if (day == -1) or (month == -1) or (year == -1):  # No date parsed
    return []

  # Now do some test on the range of the values - - - - - - - - - - - - - - - -
  #
  if ((year % 4) == 0) and ((year % 100) == 0) and ((year % 400) == 0):
    leap_year = 1
  else:
    leap_year = 0

  valid = 1  # Set date to valid

  if (month == 2):
    if (leap_year == 1) and (day > 29):
      valid = -1  # Illegal day number in February in a leap year
    elif  (day > 28):
      valid = -1  # Illegal day number in February in a normal year
  elif (month in [1,3,5,7,8,10,12]) and (day > 31):
      valid = -1  # Illegal day number in 31-day months
  elif (day > 30):
      valid = -1  # Illegal day number in 30-day months

  if (valid == -1):
     return []
  else:
    return[day,month,year]

# -----------------------------------------------------------------------------

def test():

  test_dates = [['Sep 1, 68',      '18 Jan 2002'], \
                ['17:2:2002',      '2002-02-25'], \
                ['18,03,2001',     '21.12.1999'], \
                ['February 18,19', '23\\July\\1968'], \
                ['18-02-2002',     '5/03/01'], \
                ['19680429',       '600810'],
                ['3:05:2000',      '30.11.1989'], \
                ['01011970',       "1. January '70"], \
                ['10019170',       '01011970'], \
                ['10011790',       '01011970'], \
                ['13 Feb 1945',    'Feb 13, \'45'], \
                ['April 29 1968',  '29-4=68'], \
                ['11-01-1972',     'January 10. 1972'], \
                ['23 May 1934',    '28 March 34'], \
                ['11 Jun 1899',    '11 Jul 1989'], \
                ['12111968', '      21111969'], \
                ]

  for datepair in test_dates:
    date1 = parse_datestr(datepair[0])
    date2 = parse_datestr(datepair[1])
    print 'Date1:', datepair[0], '->', date1
    print 'Date2:', datepair[1], '->', date2

    [day_diff, perc_diff] = date_diff(date1,date2)
    [tr,sub] = date_comp(date1,date2)
    age1 = date_to_age(date1[0],date1[1],date1[2])
    age2 = date_to_age(date2[0],date2[1],date2[2])
    
    print '  Day difference:       ', day_diff
    print '  Percentage difference:', perc_diff
    print '  Transpositions:       ', tr
    print '  Substitutions:        ', sub
    print '  Age 1:                ', age1
    print '  Age 2:                ', age2
    #print

    fw = date_linkage_weight(date1,date2)
    print 'Linkage weight:', fw
    print
    print

# =============================================================================
