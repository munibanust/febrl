# =============================================================================
# stringcmp.py - Several approximate string comparison routines.
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
# The Original Software is "stringcmp.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module stringcmp.py - Several approximate string comparison routines.

   Provides routines for several approximate string comparisons.

   PUBLIC FUNCTIONS:
     compare  Approximate string comparison
     test     Sime test routine.

   Possible comparison modes (input argument for compare routine) are:
     jaro      Jaro
     winkler   Winkler (based on Jaro)
     bigram    Bigram based
     editdist  Edit-distance (or Levenshtein distance) 

   See doc strings of individual functions for detailed documentation.
"""

# -----------------------------------------------------------------------------

import config
import inout

# -----------------------------------------------------------------------------

def compare(str1, str2, mode):
  """Approximate string comparison.

  USAGE:
    score = compare(str1, str2, mode)

  ARGUMENTS:
    str1   The first string
    str2   The second string
    mode   One of the possible comparison modes:
           'jaro', 'winkler', 'editdist', or 'bigram'

  DESCRIPTION:
    This routine computes and returnes a scaled comparison value between 0.0
    (strings differ totally) and 1.0 (both strings are the same).

  EXAMPLES:
    score1 = compare('meier', 'meyer', 'bigram')
    score2 = compare('meier', 'meyer', 'winkler')
  """

  if (mode == 'jaro'):
    v = jaro(str1, str2)
  elif (mode == 'winkler'):
    v = winkler(str1, str2)
  elif (mode == 'editdist'):
    v = editdist(str1, str2)
  elif (mode == 'bigram'):
    v = bigram(str1, str2)
  else:
    inout.log_message('Illegal string comparison method: '+mode,'err')
    raise Exception(compare)  # Illegal mode value

  return v

# -----------------------------------------------------------------------------

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

  len1 = len(str1)
  len2 = len(str2)
  if (len1 > len2):
    halflen = len1/2 +1 #-1
  else:
    halflen = len2/2 +1 #-1

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
    end   = min(i+halflen,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1=common1+1
      ass1=ass1+str1[i]
      workstr2 = workstr2[:index]+'*'+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2=common2+1
      ass2=ass2+str2[i]
      workstr1 = workstr1[:index]+'*'+workstr1[index+1:]

  if (common1 != common2):
    return -1.0

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition=0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition=transposition+1
  transposition=transposition/2.0

  msg = ['Jaro string comparator:', \
         '  str1: '+str1+', len1='+str(len1)+', assigned1='+ass1, \
         '  str1: '+str2+', len2='+str(len2)+', assigned2='+ass2, \
         '  halflen='+str(halflen)+', common='+str(common1)]
  inout.log_message(msg,'v2')

  w = 1./3.*(float(common1)/float(len1) + float(common1)/float(len2) + \
           (float(common1)-transposition)/float(common1))
  return w

# -----------------------------------------------------------------------------

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

  len1 = len(str1)
  len2 = len(str2)
  #halflen = max(len1,len2)/2 + 1

  if (len1 > len2):
    halflen = len1/2 +1 #-1
  else:
    halflen = len2/2 +1 #-1

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
    end   = min(i+halflen,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1=common1+1
      ass1=ass1+str1[i]
      workstr2 = workstr2[:index]+'*'+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2=common2+1
      ass2=ass2+str2[i]
      workstr1 = workstr1[:index]+'*'+workstr1[index+1:]

  if (common1 != common2):
    return -1.0

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition=0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition=transposition+1
  transposition=transposition/2.0

  # Now compute how many characters are common at beginning - - - - - - - - - -
  #
  minlen = min(len1,len2)
  for same in range(minlen+1):
    if (str1[:same] != str2[:same]):
      break
  same = same-1
  if (same > 4):
    same = 4

  msg = ['Winkler string comparator:', \
         '  str1: '+str1+', len1='+str(len1)+', assigned1='+ass1, \
         '  str1: '+str2+', len2='+str(len2)+', assigned2='+ass2, \
         '  halflen='+str(halflen)+', common='+str(common1), \
         '  number of same characters at beginning='+str(same)]
  inout.log_message(msg,'v2')

  w = 1./3.*(float(common1)/float(len1) + float(common1)/float(len2) + \
           (float(common1)-transposition)/float(common1))

  wn = w + same*0.1 * (1.0 - w)

  return wn

# -----------------------------------------------------------------------------

#def orls(str1, str2):
#  """Return approximate string comparator measure (between 0.0 and 1.0)
#
#  USAGE:
#    score = orls(str1, str2)
#
#  ARGUMENTS:
#    str1  The first string
#    str2  The second string
# 
#  DESCRIPTION:
#   As desribed in 'Methods for Automatic Record Matching and Linking and their
#    use in national Statistics", Leicester Gill, National Statistics,
#    London 2001, page 74, box 7.11
#
#    Modified to return a value between 0.0 and 1.0.
#  """
#
#  len_min = min(len(str1),len(str2))  # Length of shorter string
#
#  aggree = 0  # Number of characters that aggree
#
#  ... to be finished, PC
#
#  return w

# -----------------------------------------------------------------------------

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

  bigr1 = []
  bigr2 = []

  # Make a list of bigrams for both strings - - - - - - - - - - - - - - - - - -
  #
  for i in range(1,len(str1)):
    bigr1.append(str1[i-1:i+1])
  for i in range(1,len(str2)):
    bigr2.append(str2[i-1:i+1])

  # Get common bigrams  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  common=0.0
  for b in bigr1:
    if (b in bigr2):
      common=common+1.0

  # Compute average number of bigrams - - - - - - - - - - - - - - - - - - - - -
  #
  average = (len(bigr1)+len(bigr2)) / 2.0
  if (average == 0):
    return 0.0

  msg = ['Bigram string comparator:', \
         '  str1: '+str1+' with bigrams: '+str(bigr1), \
         '  str2: '+str2+' with bigrams: '+str(bigr2), \
         '  common='+str(common)+', average='+str(average)]
  inout.log_message(msg,'v2')

  w = common / average
  return w

# -----------------------------------------------------------------------------

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

  n = len(str1)
  if (n == 0):
    return 0.0
  m = len(str2)
  if (m == 0):
    return 0.0

  d = range(n+1)
  for i in range(n+1):
    d[i] = range(m+1)
    for j in range(m+1):
      d[i][j] = 0

  for i in range(n+1):
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

  msg = ['Edit distance string comparator:', \
         '  str1: '+str1+', len1='+str(n), \
         '  str2: '+str2+', len2='+str(m), \
         '  edit distance='+str(d[n][m])]
  inout.log_message(msg,'v2')

  return float(max(n,m) - d[n][m]) / float(max(n,m))

# -----------------------------------------------------------------------------

def test():
  """Simple test routine.

  Most test strings are taken from:
    Approximate String Comparison and its Effect on an Advanced Record
    Linkage System, Edward H. Porter and William W. Winkler, Bureau of
    Census, 1997. Research report RR97/02.
  """

  strings = [['shackleford','dunningham','nichleson','jones','massey', \
              'abroms','hardin','itman','jeraldine','marhta','michelle', \
              'julies','tanya','dwayne','sean','jon','jon','brookhaven', \
              'brook hallow','decatur','fitzrureiter','higbee','higbee', \
              'lacura','iowa','1st','peter','abcde','yz','cunningham', \
              'campell','galloway','frederick','michele','jesse', \
              'jonathon','julies','yvette','dickson','dixon','peter'], \
             ['shackelford','cunnigham','nichulson','johnson','massie', \
              'abrams','martinez','smith','geraldine','martha','michael', \
              'julius','tonya','duane','susan','john','jan','brrokhaven', \
              'brook hllw','decatir','fitzenreiter','highee','higvee', \
              'locura','iona','ist','peter','fghij','abcdef', \
              'cunnigham','campbell','calloway','fredrick','michelle', \
              'jessie','jonathan','juluis','yevett','dixon','dickson', \
              'ole']]

  print '       String 1            String 2       Jaro   Winkler    Bigram',
  print '    Edit distance'

  for i in range(len(strings[0])):
    str1 = strings[0][i]
    str2 = strings[1][i]
    print '%15s     %15s      %.3f     %.3f     %.3f     %.3f' % \
          (str1, str2, compare(str1, str2, 'jaro'), \
           compare(str1, str2, 'winkler'), \
           compare(str1, str2, 'bigram'), \
           compare(str1, str2, 'editdist'))

# -----------------------------------------------------------------------------
