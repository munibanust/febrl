# =============================================================================
# tcsv.py - Module for parsing malformed CSV and other delimited files
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
# The Original Software is "tcsv.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module tcsv.py - Module for parsing malformed CSV and other delimited files.
   Author: Tim Churches
"""

# -----------------------------------------------------------------------------

import re  # Regular expression module

# -----------------------------------------------------------------------------

class delimited_parser:
  """A parser for comma-separated value (CSV) and other delimited strings.

  # Maybe change 'null_value' to 'missing' ?? (PC, 22/8/2002)
  """

  def __init__(self, delimiter_chars=',', quote_chars='"', escape_chars='/', \
               comment_chars='#', null_chars='.', \
               null_strings =['null','None','none','missing'], null_value='', \
               as_strings=0):
    """Creates an instance of a delimited string parser, which parses a 
       delimited string into a into a sequence of fields.

    ATTRIBUTES:
      delimiter_chars  A string of characters each of which may act as a field
                       delimiter
      quote_chars      A string of characters, each of which may act as a
                       quoting character enclosing a field. The first
                       character in this string is substituted for all the
                       others in the output.
      escape_chars     A string of characters, each of which may act as an
                       escaping character for delimiter chacters embedded in
                       fields.
      comment_chars    A string of characters, each of which signifies a
                       comment line if it appears as the first character on
                       that line (that is, as the first character in line).
      null_chars       A string of characters, each of which is interpreted
                       as signifying a null or missing value in line.
      null_strings     A sequence of strings, each of which is interpreted as
                       signifying a null or missing value in line.
      null_value       The value to return when an element of null_chars or
                       null_strings is encountered.
      as_strings       If true, return each field as a string, not as an
                       appropriate data type.
    """

    self.delimiter_chars=delimiter_chars
    self.quote_chars=quote_chars
    self.escape_chars=escape_chars
    self.comment_chars=comment_chars
    self.null_chars=null_chars
    self.null_strings=null_strings
    self.null_value=null_value
    self.as_strings=as_strings

  # ---------------------------------------------------------------------------

  def typefld(self,fld):
    """Returns an appropriate type of value depending on contents of fld,
       which should be a string. Utility function used by parse_delimited().

    ARGUMENTS:
      fld  A string to be returned as the appropriate type
    """

    try:
      return int(fld) 
    except:
      try:
        return long(fld)
      except:
        try:
          return float(fld)
        except:
          fld = fld.strip()
          # print "fld:", fld
          if (fld == ''):
            return self.null_value
          elif (fld in self.null_strings):
            return self.null_value
          elif (len(fld) == 1) and (fld in self.null_chars):
            return self.null_value
          else: 
            return fld

  # ---------------------------------------------------------------------------

  def parse(self,line):
    """Parses the argument line into a sequence of fields.
    """
    
    delimiter_pattern = re.compile('[' + self.delimiter_chars + ']')
    quote_pattern = re.compile('[' + self.quote_chars + ']')
    embedded_quote_pattern = re.compile('^(.*)(?:[' + self.quote_chars + \
                             '])(.*)(?:[' + self.quote_chars + '])(.*)$')

    single_embedded_quote_pattern = re.compile('^(.*)(?:[' + \
                                               self.quote_chars + '])(.*)$')

    even_dequote_pattern = re.compile('(?:^[' + self.quote_chars + \
                           '])(.*)(?:[' + self.quote_chars + ']$)')

    odd_opening_dequote_pattern = re.compile('(?:^[' + self.quote_chars + \
                                  '])(.*)')
    odd_closing_dequote_pattern = re.compile('(.*)(?:[' + self.quote_chars + \
                                  ']$)')

    escape_pattern = '[' + self.escape_chars + ']$'
    comment_pattern = '^[' + self.comment_chars + ']'

    if (re.search(comment_pattern, line) == None):
      n_quotes_line = len(re.findall(quote_pattern,line))
      if ((n_quotes_line % 2) == 0):  # Balanced quotes - good!
        balanced_quotes = 1
      else:
        balanced_quotes = 0

      raw_parsed_line = re.split(delimiter_pattern,line)
      parsed_line = []
      quoted  = 0
      escaped = 0

      for fld in raw_parsed_line:
        n_quotes = len(re.findall(quote_pattern,fld))

        if (n_quotes % 2 == 0) and (n_quotes > 1):
          fld = fld.strip()
          if (len(re.findall(even_dequote_pattern,fld)) > 0):
            fld = re.findall(even_dequote_pattern,fld)[0]
            parsed_line.append(fld)
          elif (len(re.findall(embedded_quote_pattern,fld)) > 0):
            pfld = re.findall(embedded_quote_pattern,fld)[0]
            if (len(pfld) == 3):
              fld = pfld[0] + '"' + pfld[1] + '" ' + pfld[2] 
              parsed_line.append(fld.strip())
          else:
            parsed_line.append(self.typefld(fld))

        elif (n_quotes % 2 == 1) and (quoted == 0):
          fld = fld.lstrip()
          if (len(re.findall(odd_opening_dequote_pattern,fld)) > 0):
            fld = re.findall(odd_opening_dequote_pattern,fld)[0]
            parsed_line.append(fld)
            quoted = 1
          elif (len(re.findall(single_embedded_quote_pattern,fld)) > 0):
            fld = re.findall(single_embedded_quote_pattern,fld)[0]
            parsed_line.append(fld[0] + self.quote_chars[0] + fld[1])

        elif (n_quotes % 2 == 1) and (quoted == 1):
          fld = fld.rstrip()
          parsed_line[-1] = parsed_line[-1] + "," + \
                            re.findall(odd_closing_dequote_pattern,fld)[0]
          quoted = 0

        elif (n_quotes == 0) and (quoted == 1):
          if (balanced_quotes == 1):
            parsed_line[-1] = parsed_line[-1] + "," + fld
          else:
            parsed_line.append(fld)
            quoted = 0

        elif (re.search(escape_pattern, fld) <> None) and (escaped == 0):
          parsed_line.append(fld[:-1])
          escaped = 1

        elif (escaped == 1):
          if re.search(escape_pattern, fld) == None:
            escaped = 0
            parsed_line[-1] = parsed_line[-1] + "," + fld
          else:
            parsed_line[-1] = parsed_line[-1] + "," + fld[:-1]

        else: 
          parsed_line.append(self.typefld(fld))
        # print fld, quoted

    else:
      parsed_line = None

    if (self.as_strings):
      stringified = []
      for fld in parsed_line:
        stringified.append(str(fld))
      return stringified

    else:
      return parsed_line

  # ---------------------------------------------------------------------------

  def testparse(self,line):
    print "Input:", line
    print "Output:", self.parse(line)
    print

# ---------------------------------------------------------------------------
# Various test examples, call module as command line argument: 'python tcsv.py'

if __name__ == "__main__":

  p = delimited_parser()  # Get a parsing object

  print "With defaults"  #  - - - - - - - - - - - - - - - - - - - - - - - - - -
  p.testparse('1,  2,"3",4,"five, cinque or lima",6," seven ",8,8.5,9,0')

  print "With as_strings=1"  #  - - - - - - - - - - - - - - - - - - - - - - - -
  p.as_strings=1
  p.testparse('1,  2,"3",4,"five, cinque or lima",6," seven ",8,8.5,9,0')

  print "With defaults"  #  - - - - - - - - - - - - - - - - - - - - - - - - - -
  p.as_strings=0
  p.testparse('1,  2,"3  ",  4   , "five, cinque or lima",6," seven ",8,9,0')
  p.testparse('1,  2,"3  ",  4   , "five, cinque or lima,6," seven ",8,9,0')
  p.testparse('1,  2,"3  ",  4   , five, cinque or lima",6," seven ",8,9,0')

  print "With delimiter_chars='\t'"  #  - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars='\t'
  p.testparse('1\t2\t3\tHello world\t4')
  p.delimiter_chars='\t'
  p.testparse('1,  2,"3  ",  4   , five/, cinque or lima,6," seven ",8,9,0')

  print "With delimiter_chars=',;'"   # - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars=',;'
  p.testparse('1,2;3,4;5;6,7,8,9;0')
  p.testparse('1,2;.,4;5;6,null,8,9;0')
  p.testparse('1,2;.,4;5;6,null,,98765432198765432198721987654321987654321;0')

  print "With defaults"  #  - - - - - - - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars=','
  p.testparse('one, two,three, "four,five,six,,,,seven", eight,"nine , ten",eleven')
  p.testparse('one, two,three, four/,five/,six/,/,/,/,seven, eight,"nine , ten",eleven')
  p.testparse('#one, two,three, four/,five/,six/,/,/,/,seven, eight,"nine , ten",eleven')
  p.testparse('1,  2,"3  ",  4   , "five cinque"" or lima",6," seven ",8,9,0')
  p.testparse('1,  2,"3  ",  4   , five, cinque"" or lima,6," seven ",8,9,0')

  print "With delimiter_chars=' '"  # - - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars=' '
  p.testparse('1 2 3 4 5 6 7 8 9 0')
  p.testparse('words  delimited  by  double  spaces')

  print "With delimiter_chars='  '"  #  - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars='  '
  p.testparse('words  delimited  by  double  spaces')

  print "With defaults" # - - - - - - - - - - - - - - - - - - - - - - - - - - -
  p.delimiter_chars=','
  p.testparse('1988,"ayr" gobi rd,ulan bator,mongolia,3456')
  p.testparse('1988,""ayr"" gobi rd,ulan bator,mongolia,3456')
  p.testparse(' 1988, 13th flr "tumbi-umbi st,  ulan bator  ,mongolia, 3456 ')
  p.testparse('1996,"old mill culkinny rd,brookhome,Missing,2333,')

  print "With quote_chars set to double and single quotes"  # - - - - - - - - -
  p.quote_chars='"' + "'"
  p.testparse('1,  2,"3",4,' + "'five, cinque or lima',6,' seven ',8,8.5,9,0")
  p.testparse('1988,"ayr'+ "' gobi rd,ulan bator,mongolia,3456")
