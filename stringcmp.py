# =============================================================================
# stringcmp.py - Several approximate string comparison routines.
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
# The Original Software is "stringcmp.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module stringcmp.py - Several approximate string comparison routines.

   Provides routines for several approximate string comparisons. All return
   a value between 1.0 (the strings are the same) and 0.0 (strings are totally
   different).

   ROUTINES
     jaro      Jaro
     winkler   Winkler (based on Jaro)
     bigram    Bigram based
     editdist  Edit-distance (or Levenshtein distance) 
     seqmatch  Uses Python's standard library 'difflib'

   See doc strings of individual functions for detailed documentation.

   *** Note that for the 'jaro' and 'winkler' routines the charcter '*' must
   *** not be in either strings, as it is used as as a marker.

   If called from command line, a test routine is run which prints example
   approximate string comparisons for various string pairs.
"""

# =============================================================================
# Imports go here

import difflib

# =============================================================================

def jaro(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = jaro(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string
 
  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)
  halflen = max(len1, len2) / 2 + 1

  ass1 = ''  # Characters assigned in str1
  ass2 = ''  # Characters assigned in str2

  workstr1 = str1  # Copy of original string
  workstr2 = str2

  common1 = 0  # Number of common characters
  common2 = 0

  # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len1):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1 += 1
      ass1 = ass1+str1[i]
      workstr2 = workstr2[:index]+'*'+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2 += 1
      ass2 = ass2 + str2[i]
      workstr1 = workstr1[:index]+'*'+workstr1[index+1:]

  if (common1 != common2):
    print 'error:Jaro: Something is wrong. String 1: "%s", string2: "%s"' % \
          (str1, str2) + ', common1: %i, common2: %i' % (common1, common2) + \
          ', common should be the same.'
    common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition = 0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition += 1
  transposition = transposition / 2.0

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Jaro comparator string 1: "%s", string 2: "%s"' % (str1, str2)
  print '3:    Common: %i' % (common1)
  print '3:    Assigned 1: %s, assigned 2: %s' % (ass1, ass2)
  print '3:    Transpositions: %i' % (transposition)
  print '3:    Final approximate string weight: %f' % (w)

  return w

# =============================================================================

def winkler(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = winkler(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string
 
  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.

    Based on the 'jaro' string comparator, but modifies it according to wether
    the first few characters are the same or not.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)
  halflen = max(len1,len2) / 2 + 1

  ass1 = ''  # Characters assigned in str1
  ass2 = ''  # Characters assigned in str2
  workstr1 = str1
  workstr2 = str2

  common1 = 0  # Number of common characters
  common2 = 0

  # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len1):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1 += 1
      ass1 = ass1 + str1[i]
      workstr2 = workstr2[:index]+'*'+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2 += 1
      ass2 = ass2 + str2[i]
      workstr1 = workstr1[:index]+'*'+workstr1[index+1:]

  if (common1 != common2):
    print 'error:Winkler: Something is wrong. String 1: "%s"' % (str1) + \
          ', string2: "%s", common1: %i, common2: %i' % \
          (str1, common1, common2) + ', common should be the same.'
    common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition = 0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition += 1
  transposition = transposition / 2.0

  # Now compute how many characters are common at beginning - - - - - - - - - -
  #
  minlen = min(len1,len2)
  for same in range(minlen+1):
    if (str1[:same] != str2[:same]):
      break
  same -= 1
  if (same > 4):
    same = 4

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  wn = w + same*0.1 * (1.0 - w)

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Winkler comparator string 1: "%s", string 2: "%s"' % (str1, str2)
  print '3:    Common: %i' % (common1)
  print '3:    Assigned 1: %s, assigned 2: %s' % (ass1, ass2)
  print '3:    Transpositions: %i' % (transposition)
  print '3:    Same at beginning: %i' % (same)
  print '3:    Jaro weight: %f ' % (w)
  print '3:    Final approximate string weight: %f' % (wn)

  return wn

# =============================================================================

def bigram(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using bigrams.

  USAGE:
    score = bigram(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string
 
  DESCRIPTION:
    Bigrams are two-character sub-strings contained in a string. For example,
    'peter' contains the bigrams: pe,et,te,er.

    This routine counts the number of common bigrams and divides by the
    average number of bigrams. The resulting number is returned.
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  bigr1 = []
  bigr2 = []

  # Make a list of bigrams for both strings - - - - - - - - - - - - - - - - - -
  #
  for i in range(1,len(str1)):
    bigr1.append(str1[i-1:i+1])

  for i in range(1,len(str2)):
    bigr2.append(str2[i-1:i+1])

  # Compute average number of bigrams - - - - - - - - - - - - - - - - - - - - -
  #
  average = (len(bigr1)+len(bigr2)) / 2.0
  if (average == 0.0):
    return 0.0

  # Get common bigrams  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  common = 0.0

  if (len(bigr1) < len(bigr2)):  # Count using the shorter bigram list
    short_bigr = bigr1
    long_bigr  = bigr2
  else:
    short_bigr = bigr2
    long_bigr  = bigr1

  for b in short_bigr:
    if (b in long_bigr):
      common += 1.0
      long_bigr[long_bigr.index(b)] = []  # Mark this bigram as counted

  w = common / average

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Bigram comparator string 1: "%s", string 2: "%s"' % (str1, str2)
  print '3:    Number of bigrams 1: %i' % (len(bigr1))
  print '3:    Number of bigrams 2: %i' % (len(bigr2))
  print '3:    Average: %i' % (average)
  print '3:    Common: %i' % (common)
  print '3:    Final approximate string weight: %f' % (w)

  return w

# =============================================================================

def editdist(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the edit (or Levenshtein) distance.

  USAGE:
    score = editdist(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string
 
  DESCRIPTION:
    The edit distance is the minimal number of insertions, deletions and
    substitutions needed to make two strings equal.

    For more information on the modified Soundex see:
    - http://www.nist.gov/dads/HTML/editdistance.html
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  if (n == 0) or (m == 0):  # Check if strings are of length zero
    return 0.0

  d = range(n+1)  # Create matrix
  for i in range(n+1):
    d[i] = range(m+1)
    for j in range(m+1):
      d[i][j] = 0

  for i in range(n+1):  # Set initial values
    d[i][0] = i
  for j in range(m+1):
    d[0][j] = j

  for i in range(1,n+1):
    s = str1[i-1]

    for j in range(1,m+1):
      t = str2[j-1]
      if (s == t):
        cost = 0
      else:
        cost = 1

      d[i][j] = min(d[i-1][j]+1, d[i][j-1]+1,d[i-1][j-1]+cost)

  w = float(max(n,m) - d[n][m]) / float(max(n,m))

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Edit Distance comparator string 1: "%s", string 2: "%s"' % \
        (str1, str2)
  print '3:    n: %i' % (n)
  print '3:    m: %i' % (m)
  print '3:    Final approximate string weight: %f' % (w)

  return w

# =============================================================================

def seqmatch(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the Python standard library 'difflib' sequence matcher.

     Because the matches is not commutative, the pair and the swapped pair are
     compared and the average is taken.

  USAGE:
    score = seqmatch(str1, str2)

  ARGUMENTS:
    str1  The first string
    str2  The second string
 
  DESCRIPTION:
    For more information on Python's 'difflib' library see:

      http://www.python.org/doc/current/lib/module-difflib.html
  """

  # Quick check if the strings are the same - - - - - - - - - - - - - - - - - -
  #
  if (str1 == str2):
    return 1.0

  seq_matcher_1 = difflib.SequenceMatcher(None, str1, str2)
  seq_matcher_2 = difflib.SequenceMatcher(None, str2, str1)

  w = (seq_matcher_1.ratio()+seq_matcher_2.ratio()) / 2.0  # Return average

  # A log message for high volume log output (level 3) - - - - - - - - - - - -
  #
  print '3:  Seq Match comparator string 1: "%s", string 2: "%s"' % \
        (str1, str2)
  print '3:    Ratio 1: %f' % (seq_matcher_1.ratio())
  print '3:    Ratio 2: %f' % (seq_matcher_2.ratio())
  print '3:    Final approximate string weight: %f' % (w)

  return w

# =============================================================================
#
# Do some tests if called from command line
#
# Most test strings are taken from:
#   Approximate String Comparison and its Effect on an Advanced Record
#   Linkage System, Edward H. Porter and William W. Winkler, Bureau of
#   Census, 1997. Research report RR97/02.
#

if (__name__ == '__main__'):

  msg = []

  msg.append('Febrl module "stringcmp.py"')
  msg.append('---------------------------')
  msg.append('')

  strings = [['shackleford','dunningham','nichleson','jones','massey', \
              'abroms','hardin','itman','jeraldine','marhta','michelle', \
              'julies','tanya','dwayne','sean','jon','jon','brookhaven', \
              'brook hallow','decatur','fitzrureiter','higbee','higbee', \
              'lacura','iowa','1st','peter','abcde','yz','cunningham', \
              'campell','galloway','frederick','michele','jesse', \
              'jonathon','julies','yvette','dickson','dixon','peter', \
              'gondiwindi'], \
             ['shackelford','cunnigham','nichulson','johnson','massie', \
              'abrams','martinez','smith','geraldine','martha','michael', \
              'julius','tonya','duane','susan','john','jan','brrokhaven', \
              'brook hllw','decatir','fitzenreiter','highee','higvee', \
              'locura','iona','ist','peter','fghij','abcdef', \
              'cunnigham','campbell','calloway','fredrick','michelle', \
              'jessie','jonathan','juluis','yevett','dixon','dickson', \
              'ole', 'gondiwindiro']]

  msg.append('       String 1            String 2       Jaro   Winkler'+ \
             '    Bigram    Edit distance  Sequence matcher')

  for i in range(len(strings[0])):
    str1 = strings[0][i]
    str2 = strings[1][i]

    s = '%15s     %15s      %.3f     %.3f' %(str1, str2, jaro(str1, str2), \
                                             winkler(str1, str2))+ \
        '     %.3f     %.3f           %.3f'% (bigram(str1, str2), \
                            editdist(str1, str2),seqmatch(str1, str2))
    msg.append(s)
  for m in msg:
    print m

# =============================================================================
