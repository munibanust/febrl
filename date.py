# =============================================================================
# date.py - Routines for date conversions and parsing.
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
# The Original Software is "date.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
# - Peter Viechnicki, Research Department, Vredenburg Corp. Lanham, MD, USA
#   E-Mail: pviechnicki@vredenburg.com
#   Fixed bug in str_to_date, added '%U' directive
#
# =============================================================================

"""Module date.py - Routines for date conversions and parsing.

   PUBLIC FUNCTIONS:
     epoch_to_date        Convert a epoch day number into a date
     date_to_epoch        Convert a date into a epoch day integer
     date_to_age          Convert a date into an age (relative to a fix date)
     str_to_date          A routine that converts a string into a date using
                          a format string
     get_today            Return today's date as a [day,month,year] tuple.

   The date conversion routines are based on the 'normalDate.py' module by
   Jeff Bauer, see:

      http://starship.python.net/crew/jbauer/normalDate/
   
   Note that the epoch date used is NOT the UNIX epoch date (1 January 1970),
   but instead the 1 January 1900.

   Note that dates are returned as a tuple of strings, with day and month being
   of length 2 (i.e. '01' etc.), and year being of length 4 (e.g. '2003').

   See doc strings of individual functions for detailed documentation.
"""

# =============================================================================
# Imports go here

import math
import string
import time

# =============================================================================

# A dictionary of month name abbreviations, used in date.str_to_date() routine

month_abbrev_dict = {'jan':'01', 'feb':'02', 'mar':'03', 'apr':'04', \
                     'may':'05', 'jun':'06', 'jul':'07', 'aug':'08', \
                     'sep':'09', 'oct':'10', 'nov':'11', 'dec':'12'}

# Define a character replace table for data strings - - - - - - - - - - - -
#
string_replace = ["'.,:-=_/\\", \
                  "         "]

# Characters in the first list are replaced by the corresponding character in
# the second list

replace_table = string.maketrans(string_replace[0], string_replace[1])

# =============================================================================
# Some simple funcions used for date conversions follow
# (based on functions from the 'normalDate.py' module by Jeff Bauer, see:
# http://starship.python.net/crew/jbauer/normalDate/)

days_in_month = [[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], \
                 [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]]

def first_day_of_year(year):
  """Calculate the day number (relative to 1 january 1900) of the first day in
     the given year.
  """

  if (year == 0):
    print 'error:A year value of 0 is not possible'
    raise Exception

  elif (year < 0):
    first_day = (year * 365) + int((year - 1) / 4) - 693596
  else:  # Positive year
    leap_adj = int ((year + 3) / 4)
    if (year > 1600):
      leap_adj = leap_adj - int((year + 99 - 1600) / 100) + \
                 int((year + 399 - 1600) / 400)

    first_day = year * 365 + leap_adj - 693963

    if (year > 1582):
      first_day -= 10

  return first_day

# -----------------------------------------------------------------------------

def is_leap_year(year):
  """Determine if the given year is a leap year. Returns 0 (no) or 1 (yes).
  """

  if (year < 1600):
    if ((year % 4) != 0):
      return 0
    else:
      return 1

  elif ((year % 4) != 0):
    return 0

  elif ((year % 100) != 0):
    return 1

  elif ((year % 400) != 0):
    return 0

  else:
    return 1

# =============================================================================

def epoch_to_date(daynum):
  """Convert an epoch day number into a date [day, month, year], with
     day, month and year being strings of length 2, 2, and 4, respectively.
     (based on a function from the 'normalDate.py' module by Jeff Bauer, see:
     http://starship.python.net/crew/jbauer/normalDate/)

  USAGE:
    [year, month, day] = epoch_to_date(daynum)

  ARGUMENTS:
    daynum  A integer giving the epoch day (0 = 1 January 1900)

  DESCRIPTION:
    Function for converting a number of days (integer value) since epoch time
    1 January 1900 (integer value) into a date tuple [day, month, year].

  EXAMPLES:
    [day, month, year] = epoch_to_date(0)      # returns ['01','01','1900']
    [day, month, year] = epoch_to_date(37734)  # returns ['25','04','2003']
  """

  if (not (isinstance(daynum, int) or isinstance(daynum, long))):
    print 'error:Input value for "daynum" is not of integer type: %s' % \
          (str(daynum))
    raise Exception

  if (daynum >= -115860):
    year = 1600 + int(math.floor((daynum + 109573) / 365.2425))
  elif (daynum >= -693597):
    year = 4 + int(math.floor((daynum + 692502) / 365.2425))
  else:
    year = -4 + int(math.floor((daynum+695058) / 365.2425))

  days = daynum - first_day_of_year(year) + 1

  if (days <= 0):
    year -= 1  
    days = daynum - first_day_of_year(year) + 1

  days_in_year = 365 + is_leap_year(year)  # Adjust for a leap year

  if (days > days_in_year):
    year += 1
    days = daynum - first_day_of_year(year) + 1

  # Add 10 days for dates between 15 October 1582 and 31 December 1582
  #
  if (daynum >= -115860) and (daynum <= -115783):
    days += 10

  day_count = 0
  month = 12
  leap_year_flag = is_leap_year(year)

  for m in range(12):
    day_count += days_in_month[leap_year_flag][m]
    if (day_count >= days):
      month = m + 1
      break

  # Add up the days in the prior months
  #
  prior_month_days = 0
  for m in range(month-1):
    prior_month_days += days_in_month[leap_year_flag][m]

  day = days - prior_month_days

  day_str =   string.zfill(str(day),2)  # Add '0' if necessary
  month_str = string.zfill(str(month),2)  # Add '0' if necessary
  year_str =  str(year)  # Is always four digits long

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Epoch: %i -> Day:%s, month:%s, year:%s' % \
        (daynum, day_str, month_str, year_str)

  return [day_str, month_str, year_str]

# =============================================================================

def date_to_epoch(day, month, year):
  """ Convert a date [day, month, year] into an epoch day number.
     (based on a function from the 'normalDate.py' module by Jeff Bauer, see:
     http://starship.python.net/crew/jbauer/normalDate/)

  USAGE:
    daynum = date_to_epoch(year, month, day)

  ARGUMENTS:
    day    Day value (string or integer number)
    month  Month value (string or integer number)
    year   Year value (string or integer number)

  DESCRIPTION:
    Function for converting a date into a epoch day number (integer value)
    since 1 january 1900.

  EXAMPLES:
    day = date_to_epoch('01', '01', '1900')  # returns 0
    day = date_to_epoch('25', '04', '2003')  # returns 37734
  """

  # Convert into integer values
  #
  try:
    day_int = int(day)
  except:
    print 'error:"day" value is not an integer'
    raise Exception
  try:
    month_int = int(month)
  except:
    print 'error:"month" value is not an integer'
    raise Exception
  try:
    year_int = int(year)
  except:
    print 'error:"year" value is not an integer'
    raise Exception

  # Test if values are within range
  #
  if (year_int <= 1000):
    print 'error:Input value for "year" is not a positive integer ' + \
          'number: %i' % (year)
    raise Exception

  leap_year_flag = is_leap_year(year_int)

  if (month_int <= 0) or (month_int > 12):
    print 'error:Input value for "month" is not a possible day number: %i' % \
          (month)
    raise Exception

  if (day_int <= 0) or (day_int > days_in_month[leap_year_flag][month_int-1]):
    print 'error:Input value for "day" is not a possible day number: %i' % \
          (day)
    raise Exception

  days = first_day_of_year(year_int) + day_int - 1

  for m in range(month_int-1):
    days += days_in_month[leap_year_flag][m]

  if (year_int == 1582):
    if (month_int > 10) or ((month_int == 10) and (day_int > 4)):
      days -= 10

  return days

# =============================================================================

def date_to_age(day, month, year, fix_date='today'):
  """Convert a date into an age (relative to a fix date)

  USAGE:
    age = date_to_age(day, month, year)
    age = date_to_age(day, month, year, fix_date)

  ARGUMENTS:
    day       Day value (integer number)
    month     Month value (integer number)
    year      Year value (integer number)
    fix_date  The date relative for which the age is computed. Can be a date
              tuple, the string 'today' (which is the default), or an integer
              (epoch day number)

  DESCRIPTION:
    Returns the age in years as a positive floating-point number.
    If the date is after the fix date a negative age is returned.
  """

  # Check if fix date is given, otherwise calculate it  - - - - - - - - - - - -
  #
  if (fix_date == 'today'):
    sys_time = time.localtime(time.time())  # Get current system date
    fix_day =   string.zfill(str(sys_time[2]),2)
    fix_month = string.zfill(str(sys_time[1]),2)
    fix_year =  str(sys_time[0])

  elif (isinstance(fix_date, list) and (len(fix_date) == 3)):
    fix_day =   string.zfill(str(fix_date[0]),2)
    fix_month = string.zfill(str(fix_date[1]),2)
    fix_year =  str(fix_date[2])

  elif (isinstance(fix_date, int)):
    fix_epoch = fix_date

  else:
    print 'error:"fix_date" is not in a valid form: %s' % (str(fix_date))
    raise Exception

  # Get epoch number for input date and fix date  - - - - - - - - - - - - - - -
  #
  date_epoch = date_to_epoch(day, month, year)

  if (not isinstance(fix_date, int)):
    fix_epoch  = date_to_epoch(fix_day, fix_month, fix_year)

  day_diff = fix_epoch - date_epoch  # Get day difference

  # Compute approximate age - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  age = float(day_diff) / 365.25  # Can be positive or negative

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Date: %s with fix date: %s -> Age: %.3f' % \
        (str([day,month,year]), str(fix_date), age)

  return age

# =============================================================================

def str_to_date(date_str, format_str, pivot_year):
  """A routine that converts a string into a date using a format string.

  USAGE:
    [day,month,year] = str_to_date(datestr, formatstr)

  ARGUMENTS:
    datestr     Input date as a string
    formatstr   A format string, must be made of three directives, which can
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
                  %U  'UNK' or 'UNKNOWN' (for day)
    pivot_year  If a two-digit year is given, the pivot year is used to
                detemine if it is expanded into 19XX or 20XX. Two-digits years
                smaller than the pivot year will be expanded into 20XX, years
                larger and equal than the pivot year will be expanded into 19xx
                Example: pivot_year = 03:  68 -> 1968
                                           03 -> 1903
                                           02 -> 2002

  DESCRIPTION:
    This routine parses the input date string according to the given format
    and extracts a [day,month,year] triplet if possible. The output is a list
    of three strings, with both day and month having length 2, and year having
    length 4. Example: ['01','02','2003']

    If the routine can't parse the date string successfully an empty
    list is returned.

    Valid four digit year values must be in the interval [1850,2005].
  """

  # Apply character replace table - - - - - - - - - - - - - - - - - - - - - - -
  #
  date_str = date_str.translate(replace_table)

  # Remove leading and trailing whitespaces and make lower case - - - - - - - -
  #
  format_str = format_str.strip()
  date_str   = date_str.strip().lower()

  # Replace triple and double spaces with one space only  - - - - - - - - - - -
  #
  date_str = date_str.replace('   ',' ')
  date_str = date_str.replace('  ',' ')
  date_str = date_str.replace('  ',' ')

  # Now check the date and format strings - - - - - - - - - - - - - - - - - - -
  #
  if (date_str =='') or (format_str ==''):
    return []  # No date or no format string given

  elif (' ' in format_str) and (' ' in date_str):
    date_list   = date_str.split()
    format_list = format_str.split()

  elif (' ' not in format_str) and (' ' not in date_str):
    date_list   = []
    format_list = []

    work_format_str = format_str
    work_date_str   = date_str
    while (work_format_str != ''):
      if (work_format_str[:2] == '%Y'):  # Four digit year
        format_list.append(work_format_str[:2])
        work_format_str = work_format_str[2:]
        date_list.append(work_date_str[:4])
        work_date_str = work_date_str[4:]

      elif (work_format_str[:2] in ['%m','%d','%y']):  # 2-digit year/month/day
        format_list.append(work_format_str[:2])
        work_format_str = work_format_str[2:]
        date_list.append(work_date_str[:2])
        work_date_str = work_date_str[2:]

      else:
        print 'error:Illegal format string without spaces: "%s"' % (formatstr)
        raise Exception

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
        month = month_abbrev_dict.get(date_list[0][:3],-1)

    elif (directive == '%d'):  # Day of the month number
      try:
        day = int(date_list[0])  # Convert into an integer number
        if (day < 1) or (day > 31):
          day = -1
      except:
        day = -1  # Day is no integer

    elif (directive == '%m'):  # Month number (between 1..12)
      try:
        month = int(date_list[0])  # Convert into an integer number
        if (month < 1) or (month > 12):
          month = -1
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

        # Convert year into four digit value according to value of pivot_year
        #
        if (year >= 0) and (year < pivot_year):
          year = year+2000
        elif (year >= pivot_year) and (year < 100):
          year = year+1900

    elif (directive == '%Y'):  # Four digit year with century
      if (len(date_list[0]) == 4):  # Parse only if length is 4 digits
        try:
          year = int(date_list[0])  # Convert into an integer number
          if (year < 1850) or (year > 2005):
            year = -1  # Illegal year value in four digit
        except:
          year = -1  # year is no integer

    elif (directive == '%U'): # New date expression for UNKNOWN
                              # (Peter Viechnicki)
      # print "New %U: date_list[0] = ", date_list[0]
      if (date_list[0] == 'unk' or date_list[0] == 'unknown'):
        day = 'UNK'
      else:
        day = -1

    else:
      print 'error:Illegal format directive: "%s"' % (directive)
      raise Exception

    date_list =   date_list[1:]    # Remove processed first element
    format_list = format_list[1:]

  if (day == -1) or (month == -1) or (year == -1):  # No date parsed
    return []

  # Now do some test on the range of the values - - - - - - - - - - - - - - - -
  #
  if ((year % 4) == 0):
    if ((year % 100) == 0):
      if ((year % 400) == 0):
        leap_year = True
      else:
        leap_year = False
    else:
      leap_year = True
  else:
    leap_year = False

  valid = True

  if (day != 'UNK'):  # Only test day value if known
    if (month == 2):
      if (leap_year == True) and (day > 29):
        valid = False  # Illegal day number in February in a leap year
      elif (leap_year == False) and (day > 28):
        valid = False  # Illegal day number in February in a normal year
    elif (month in [1,3,5,7,8,10,12]) and (day > 31):
      valid = False  # Illegal day number in 31-day months
    elif (month in [4,6,9,11]) and (day > 30):
      valid = False  # Illegal day number in 30-day months

  if (valid == False):
    return []

  else:
    day_str =   string.zfill(str(day),2)
    month_str = string.zfill(str(month),2)
    year_str =  str(year)

    return [day_str,month_str,year_str]

# =============================================================================

def get_today():
  """Return today's date as a [day,month,year] tuple, with three string (both
     day and month having length 2, and year having length 4).
  """

  sys_time = time.localtime(time.time())  # Get current system time and date
  today = [string.zfill(str(sys_time[2]),2), \
           string.zfill(str(sys_time[1]),2), \
           str(sys_time[0])]

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:    Today is %s' % (str(today))

  return today

# =============================================================================
