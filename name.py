# =============================================================================
# name.py - Routines for name parsing, comparing and linkage.
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
# The Original Software is "name.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module name.py - Routines for name parsing, comparing and linkage.

   PUBLIC FUNCTIONS:
     clean_name_component  Clean a name component input string
     tag_name_component    Tag a name input component string and make a list
     get_gender_guess      Extract a gender guess from tags and words
     get_title             Extract the title component of a name
     get_name_component    Parse the input word list and extract names and
                           alternative names for either givenname or surnames
     get_names_rules       Parse the input using rules and extract given- and
                           surnames
     get_names_hmm         Process the input word and tag lists using a Hidden
                           Markov Model (HMM) to extract name output fields
     test                  Simple test routine with example inputs

   See doc strings of individual functions for detailed documentation.

   See also the relevant section in the 'config.py' module.
"""

# -----------------------------------------------------------------------------

import string
import types

import config
import inout
import mymath

# -----------------------------------------------------------------------------

def clean_name_component(name_str):
  """Clean a name component input string.

  USAGE:
    cleaned_str = clean__name_component(name_str)

  ARGUMENTS:
    name_str  A string containing the name component (or parts of it)

  DESCRIPTION:
    This routines cleans the input string by using the 'name_corr_list'. It
    also strips off all leading and trailing spaces. A cleaned string is
    returned.
  """

  # First add a trailing and leading space  - - - - - - - - - - - - - - - - - -
  # (this is to make sure replacement strings do match at beginning and end)
  #
  name_str = ' '+name_str+' '

  # Check for strings from the name correction-list - - - - - - - - - - - - - -
  #
  for (org,repl) in config.name_corr_list:
    name_str = name_str.replace(org,repl)

  # Make sure commas are separated from words so they become list elements  - -
  #
  name_str = name_str.replace(',', ' , ')

  return name_str.strip()

# -----------------------------------------------------------------------------

def tag_name_component(name_str):
  """Tag a name input component string and make a list.

  USAGE:
    [word_list, tag_list] = tag_name_component(name_str)

  ARGUMENTS:
    name_str  A string containing the name component

  DESCRIPTION:
    This routine cleans the input string and extracts words, numbers and
    separators into a list. Each element of this list is assigned one or more
    tags. A 'greedy tagger' is applied, which cheques sequences of list
    elements in the name lookup table (longer sequences first) and replaces
    them with the string and tag from the lookup-table if found.

    The routine returns two lists: words and their tags
  """

  # First, split input string into elements at spaces - - - - - - - - - - - - -
  #
  org_list = name_str.split()  # The original list from the input string
  inout.log_message('  Initial word list: '+str(org_list),'v2')

  tag_list  = []  # The initially empty list of tags
  word_list = []  # The initially empty list of words

  while (org_list != []):  # As long as not all elements have been processed
    tmp_list = org_list[:config.name_dict_seq_len]  # Extract longest sub-list
    tmp_val = []  # Start with empty value
    tmp_key = tuple(tmp_list)

    while (tmp_key != ()):  # As long as key not empty and not found in lookup
      if (config.name_lookup_dict.has_key(tmp_key)):
        tmp_val = config.name_lookup_dict[tmp_key]
        break
      tmp_key = tmp_key[:-1]  # Remove last element in key

    if (tmp_val != []):  # A value has been found in the dictionary
      tmp_len = len(tmp_key)  # Length of found sequence

      if (tmp_val[0] != ''):  # it's not an empty value
        word_list.append(tmp_val[0])  # Append corrected word (or sequence)
        tag_list.append(tmp_val[1])   # Append tag or tags

    else:  # No value has been found in the lookup dictionary, try other tags

      tmp_val = org_list[0]  # Value is first element in the original list
      tmp_len = 1

      if (len(tmp_val) == 1) and (tmp_val.isalpha()):  # A 1-letter word
        word_list.append(tmp_val)
        tag_list.append('II')

      elif (tmp_val.isdigit()):  # Element is a number
        word_list.append(tmp_val)
        tag_list.append('NU')

      elif (not tmp_val.isalpha()) and tmp_val.isalnum():  # Alpha-numeric
        word_list.append(tmp_val)
        tag_list.append('AN')

      elif (tmp_val == '-'):  # Element is a hyphen
        if (tag_list != []) and (tag_list[-1] not in ['BO','SP']):
                              # Don't append hyphen at beginning or after BO/SP
          word_list.append(tmp_val)
          tag_list.append('HY')

      elif (tmp_val == ','):  # Element is a comma
        word_list.append(tmp_val)
        tag_list.append('CO')

      elif (tmp_val == '|'):  # Element is a vertical bar
        word_list.append(tmp_val)
        tag_list.append('VB')

      else:  # An unknown element
        word_list.append(tmp_val)
        tag_list.append('UN')

    # Finally remove the processed elements from the original element list
    #
    org_list = org_list[tmp_len:]  # Remove processed elements

  return [word_list, tag_list]

# -----------------------------------------------------------------------------

def get_gender_guess(word_list, tag_list):
  """Extract a gender guess from tags and words.

  USAGE:
    gender_guess = get_gender_guess(word_list, tag_list)

  ARGUMENTS:
    word_list  List of words as produces with clean_tag_names()
    tag_list   Corresponding list of tags as produces with clean_tag_names()

  DESCRIPTION:
    First, titles are checked and if a word is found in a gender list of
    titles the gender is set.

    If no gender could be extracted from title words, the remaining list of
    tags is checked, and if a givenname gender is found, it is used.
    A final gender value is only returned if it is un-ambiguous (i.e. if
    both male and female givennames are found, no gender information is
    returned).

    The returned string is either 'female', 'male' or '' (if no gender has been
    found).
  """

  male_count   = 0
  female_count = 0

  # Check title words - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len(tag_list)):
    if (tag_list[i] == 'TI'):
      if (word_list[i] in config.name_male_title):
        male_count += 1
      elif (word_list[i] in config.name_female_title):
        female_count += 1

  # Check givenname tags only if no gender in title words found - - - - - - - -
  #
  if (male_count == 0) and (female_count == 0):
    for t in tag_list:
      if (t.find('GM') >= 0):
        male_count += 1
      elif (t.find('GF') >= 0):
        female_count += 1

  # Check if only male or only female gender found - - - - - - - - - - - - - -
  #
  if (male_count > 0) and (female_count == 0):
    gender = 'male'

  elif (male_count == 0) and (female_count > 0):
    gender = 'female'

  else:
    gender = ''

  inout.log_message('  Gender guess: '+gender,'v1')

  return gender

# -----------------------------------------------------------------------------

def get_title(word_list, tag_list):
  """Extract the title component of a name.

  USAGE:
    [word_list, tag_list, title_list] = get_title(word_list, tag_list)

  ARGUMENTS:
    word_list  List of words as produces with clean_tag_names()
    tag_list   Corresponding list of tags as produces with clean_tag_names()

  DESCRIPTION:
    This routine extracts all title words (and removes them from both the word
    and tag list).
    It returns a list with title words that were at the beginning of the word
    list, but not after a non-title word.

    Returns the modified word and tag list with all title words removed, and
    the list of extracted title words sorted alphabetically.
  """

  title_list = []
  tmp_list = []
  tmp_tags = []
  list_len = len(word_list)

  # Extract all title words (list elements with a TI tag) - - - - - - - - - - -
  #
  for i in range(list_len):
    if (tag_list[i] == 'TI'):
      title_list.append(word_list[i])
    else:
      title_list.append('')
      tmp_list.append(word_list[i])
      tmp_tags.append(tag_list[i])

  # Make a dictionary of title words found at beggining of word list  - - - - -
  #
  title_dict = {}
  i = 0
  while (i < list_len) and (title_list[i] == ''):
    i+=1
  while (i < list_len) and (title_list[i] != ''):
    title_dict.update({title_list[i]:i})
    i+=1

  title_list = title_dict.keys()
  title_list.sort()

  inout.log_message('  Extracted title list: '+str(title_list),'v1')
  inout.log_message(['  Modified word list: '+str(tmp_list), \
                     '  Modified tag list:  '+str(tmp_tags)],'v2')

  return [tmp_list, tmp_tags, title_list]

# -----------------------------------------------------------------------------

def get_name_component(word_list, tag_list, name_mode):
  """Parse the input word list and extracts names and alternative names.

  USAGE:
    [names, alt_names] = get_name_component(word_list, tag_list, name_mode)

  ARGUMENTS:
    word_list  List of words as produces with clean_tag_names()
    tag_list   Corresponding list of tags as produces with clean_tag_names()
    name_mode  Either 'gname' if the name component contains givennames, or
               'sname' if it contains surnames.

  DESCRIPTION:
    This routine can be used to parse and process a name list that contains
    either givennames or surnames (but not both of them). It uses the tag list
    to distinguish between names and alternative names (which occur after a
    separator element or within brackets).

    Two lists are returned, the first one with the names and the second one
    with the alternative names.
  """

  if (name_mode not in ['gname','sname']):
    inout.log_message('Illegal name mode: '+name_mode,'err')
    raise Exception(get_name_component)

  names_list=[[],[]]  # Output list with primary and alternative names

  # Flags (indices into 'names_list') for different name modes
  #
  name_mode_prim = 0  # Flag for 'primary' name mode
  name_mode_alt  = 1  # Flag for 'alternative' name mode
  curr_name_mode = name_mode_prim  # Set first name mode

  # Flags for vertical bars mode
  #
  vb_mode_outside = 0  # 'outside' vertical bars
  vb_mode_inside  = 1  # 'inside' vertical bars
  curr_vb_mode    = vb_mode_outside

  sep_index = -1  # Gives the index of the last separator element in input list

  list_len = len(word_list)

  # Loop over all elements in the input list  - - - - - - - - - - - - - - - - -
  #
  i = 0
  while i < list_len:
    w = word_list[i]  # Process this word
    t = tag_list[i]   # Corresponding tag

    if (t == 'ST'):  # Process a saint word or name - - - - - - - - - - - - - -
      if ('_' not in w):  # A saint word only, like saint, brother, holy, etc.
        names_list[curr_name_mode].append(w)
      else:  # A 'saint' plus name sequence, append name to alternative names
        names_list[curr_name_mode].append(w)
        unders_ind = w.find('_')
        names_list[name_mode_alt].append(w[unders_ind+1:])  # Name word only

    elif (t in ['RU','CO']):  # A 'rubbuish' word or a comma  - - - - - - - - -
      pass  # Just pass over it, don't append it to the names list

    elif (t == 'BO'):  # Process 'baby of' and similar sequences  - - - - - - -
      names_list[curr_name_mode].append(w)

    elif (t == 'PR'):  # Process name prefix  - - - - - - - - - - - - - - - - -
      names_list[curr_name_mode].append(w)

    elif (t in ['NU','AN']):  # A Number or alphanumeric  - - - - - - - - - - -
      names_list[curr_name_mode].append(w)
      inout.log_message('Number or alphanumeric word: '+str(word_list),'warn')

    elif (t == 'HY'):  # Process hyphen - - - - - - - - - - - - - - - - - - - -
      if (i > 0):
        if (tag_list[i-1] not in ['HY','BO','CO','SP','VB']):
          names_list[curr_name_mode].append(w)  # Don't always append a hyphen
        else:
          inout.log_message('Strange hyphen situation: '+str(word_list),'warn')

    elif (t == 'VB'):  # Vertical bar, switch vertical bar mode - - - - - - - -
      if (i == 0):
        inout.log_message('Strange situation: Vertical bar at beginning: '+ \
                          str(word_list),'warn')

      elif (tag_list[i-1] in ['HY','BO','CO','SP']) or \
           ((tag_list[i-1] == 'ST') and ('_' not in word_list[i-1])):
        inout.log_message('Strange vertical bar situation: '+str(word_list), \
                          'warn')
      else:  # Switch vertical bar (from outside to inside or vice versa)
        if (curr_vb_mode == vb_mode_outside):
          curr_vb_mode   = vb_mode_inside
          curr_name_mode = name_mode_alt  # Switch to alternative name mode
          # All names between vertival bars (brackets) are in alternative mode
        else:
          curr_vb_mode = vb_mode_outside
          curr_name_mode = name_mode_prim  # Switch to primary name mode

    elif (t == 'SP'):  # A separator like 'known_as', 'and' or 'or' - - - - - -
      sep_index = i  # Store index of separator list element

      # Switch to alternative name if it's 'known_as' not at beginning
      #
      if (i > 0) and (w == 'known_as'):
        curr_name_mode = name_mode_alt  # Following names are now alt. names

    elif (t == 'NE'):  # The word 'nee' (name or separator) - - - - - - - - - -
      # 'nee' is a separator if it is not at the beginning and the previous
      # word is neither a name prefix, saint name, nor a 'baby_of' sequence
      #
      if (i > 0) and ((tag_list[i-1] not in ['PR','BO']) or \
                      ((tag_list[i-1] == 'ST') and \
                        ('_' not in word_list[i-1]))):
        sep_index = i  # Get index of separator word
        curr_name_mode = name_mode_alt  # Switch to alternative name mode

    elif (t == 'II'):  # A single character (e.g. initial)  - - - - - - - - - -
      names_list[name_mode_alt].append(w)  # Append single letter to alt. name

    # Word must be a name word, tagged with: SN, GF, GM or UN - - - - - - - - -
    #
    else:

      # Check alternative mode: If we're outside vertical bars ('curr_vb_mode'
      # is 'vb_mode_outside') and the 'curr_name_mode' is 'name_mode_alt'
      # check previous list element
      #
      if (curr_vb_mode == vb_mode_outside) and \
         (curr_name_mode == name_mode_alt):
        if (tag_list[i-1] not in ['SP', 'HY', 'PR']):

          # Switch back to primary name mode if list element before is a
          # 'normal' word (assuming we're behind a 'known_as' sequence)
          #
          curr_name_mode == name_mode_prim

      names_list[curr_name_mode].append(w)  # Append word to current name mode

    i+=1

  # Now check if words are stored in both name lists, and if so remove them
  # from the alternative name list
  tmp_name_list1 = []
  for w in names_list[1]:
    if (w not in names_list[0]):
      tmp_name_list1.append(w)

  # Now check if a word appears twice in the name lists, if so remove one
  #
  tmp_name_list0 = []
  for w in names_list[0]:
    if (w not in tmp_name_list0):
      tmp_name_list0.append(w)
  tmp_name_list11 = []
  for w in tmp_name_list1:
    if (w not in tmp_name_list11):
      tmp_name_list11.append(w)

  return (tmp_name_list0, tmp_name_list11)

# -----------------------------------------------------------------------------

def get_names_rules(word_list, tag_list, first_name_comp):
  """Parse the input word list using rules and extracts given- and surnames.

  USAGE:
    name_dict = get_names(word_list, tag_list, first_name_comp)

  ARGUMENTS:
    word_list        List of words as produces with clean_tag_names()
    tag_list         Corresponding list of tags as produces with
                     clean_tag_names()
    first_name_comp  Set to 'gname' if the input is most likely to start
                     with givennames, or to 'sname' if it most likely
                     starts with surnames.

  DESCRIPTION:
    This routine can be used to parse and process a name list that contains
    both given- and surnames. It uses the tag list to distinguish between
    given- and surnames, and between names and alternative names (which occur
    after a separator element or within brackets).

    The dictionary returned can contain the following key words:
    - givenname
    - alt_givenname
    - surname
    - alt_surname
  """

  if (first_name_comp not in ['gname','sname']):
    inout.log_message('Illegal value for first_name_comp: '+first_name_comp, \
                      'err')
    raise Exception(get_names)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Phase one: Split input word and tag list into up to five sub-lists  - - - -
  # (for words and for tags) at separator elements
  #
  # Sub-list 1: All words/tags until first separator
  # Sub-list 2: All words/tags in first alternative given word area
  # Sub-list 3: All words/tags until next separator or end
  # Sub-list 4: All words/tags in second alternative given word area
  # Sub-list 5: All words/tags after second alternative word area
  #
  word_sub_list = [[],[],[],[],[]]
  tag_sub_list  = [[],[],[],[],[]]

  # Check if whe word 'nee' (with tag 'NE') is in the list, and decide if - - -
  # it is a separator or a name word
  #
  if ('NE' in tag_list):
    i = 0
    for i in range(len(tag_list)):
      if (tag_list[i] == 'NE'):
        # 'nee' is a separator if it is not at the beginning and the word
        # before 'nee' is neither a name prefix, a saint name, nor a 'baby_of'
        # sequence
        #
        if (i > 0) and ((tag_list[i-1] not in ['PR','BO']) or \
           ((tag_list[i-1] == 'ST') and ('_' not in word_list[i-1]))):
          tag_list[i] = 'SP'
        else:
          tag_list[i] == 'UN'  # Make it an unknown word

  # If no separator or vertical bar in the input then no splitting into - - - -
  # sub-lists is needed
  #
  if ('SP' not in tag_list) and ('VB' not in tag_list):
    word_sub_list[0] = word_list
    tag_sub_list[0]  = tag_list

  else:  # Perform splitting into sub-lists

    list_len = len(word_list)
    list_ptr = 0  # Pointer into five sub-lists
    i = 0  # Loop index

    while (i < list_len):

      # Handle separators ('known_as' or 'and' or 'or') - - - - - - - - - - - -
      #
      if (tag_list[i] == 'SP'):

        # Only switch to next sub-list if it's a 'known_as' or 'nee' separator
        # not at the beginning, and the element before was not a vertical bar
        #
        if (i > 0) and (word_list[i] in ['known_as','nee']) and \
           (tag_list[i-1] != 'VB'):
          if (list_ptr == 4):
            inout.log_message('To many separators/vertical bars in input: '+ \
                              str(word_list),'warn')
          else:
            list_ptr += 1  # Switch to next sub-list

        if (i == list_len-1):
            inout.log_message('Last element in input is a separator: '+ \
                              str(word_list),'warn')
        else:
          i += 1  # Go to next word/tag, i.e. skip over separator

          if (tag_list[i] not in ['SP','VB','RU','CO']):
            word_sub_list[list_ptr].append(word_list[i])  # Append next element
            tag_sub_list[list_ptr].append(tag_list[i])
          else:  # A strange situation
            inout.log_message('Strange separator situation: '+str(word_list), \
                              'warn')

          # If current word is a name prefix append following word(s) as well
          #
          if (tag_list[i] == 'PR') and (i < list_len-1):
            i += 1
            word_sub_list[list_ptr].append(word_list[i])  # Append name word
            tag_sub_list[list_ptr].append(tag_list[i])

            if (tag_list[i] == 'PR') and (i < list_len-1):  # Check 2. prefix
              i += 1
              word_sub_list[list_ptr].append(word_list[i])  # Append name word
              tag_sub_list[list_ptr].append(tag_list[i])

          # Check if word is hyphened with another word, if so append
          #
          elif (i < list_len-2) and (tag_list[i+1] == 'HY'):
            word_sub_list[list_ptr].append(word_list[i+1])  # Append hyphen
            tag_sub_list[list_ptr].append(tag_list[i+1])
            word_sub_list[list_ptr].append(word_list[i+2])  # Append word
            tag_sub_list[list_ptr].append(tag_list[i+2])
            i += 2  # Increase pointer in list to next element

          # If the following element is not a separator switch to next sub-list
          #
          if (i < list_len-1) and (tag_list[i+1] != 'SP'):
            if (list_ptr == 4):
              inout.log_message('To many separators/vertical bars in input: '+\
                                str(word_list),'warn')
            else:
              list_ptr += 1  # Switch to next sub-list

          # If the following word is again 'known_as' change it to 'and' so in
          # the next iteration no sub-list switching occurs
          #
          elif (i < list_len-1) and (word_list[i+1] == 'known_as'):
            word_list[i+1] = 'and'

      # Handle vertical bars  - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      elif (tag_list[i] == 'VB'):

        # Only switch to next sub-list if vertical bar is not at the beginning
        #
        if (i > 0):
          if (list_ptr == 4):
            inout.log_message('To many separators/vertical bars in input: '+\
                              str(word_list),'warn')
          else:
            list_ptr += 1  # Switch to next sub-list

        if (i != list_len-1):
          i += 1  # Go to next word/tag, i.e. skip vertical bar

        # Process all elements until next vertical bar
        #
        while (i < list_len-1) and (tag_list[i] != 'VB'):

          if (tag_list[i] != 'SP'):  # Jump over separator elements
            word_sub_list[list_ptr].append(word_list[i])
            tag_sub_list[list_ptr].append(tag_list[i])
          i += 1

        # Only switch to next sub-list if the element following the vertical
        # bar is not a separator 'known_as'
        #
        if (i < list_len-1) and (word_list[i+1] != 'known_as'):
          if (list_ptr == 4):
            inout.log_message('To many separators/vertical bars in input: '+\
                              str(word_list),'warn')
          else:
            list_ptr += 1  # Switch to next sub-list

      # Append all other words/tags to the current sub-lists  - - - - - - - - -
      #
      else:
        word_sub_list[list_ptr].append(word_list[i])
        tag_sub_list[list_ptr].append(tag_list[i])

      if (i < list_len):
        i += 1  # Go to next word/tag

  inout.log_message(['  Extracted sub-lists:', \
    '    Potential givennames:                 '+str(word_sub_list[0]), \
    '                                          '+str(tag_sub_list[0]), \
    '    Potential alternative givennames:     '+str(word_sub_list[1]), \
    '                                          '+str(tag_sub_list[1]), \
    '    Potential surnames:                   '+str(word_sub_list[2]), \
    '                                          '+str(tag_sub_list[2]), \
    '    Potential alternative surnames:       '+str(word_sub_list[3]), \
    '                                          '+str(tag_sub_list[3]), \
    '    Potential second alternative surnames:'+str(word_sub_list[4]), \
    '                                          '+str(tag_sub_list[4])],'v2')

  # Make sure there are no empty sub-lists between filled sub-lists
  #
  if (((word_sub_list[3] == []) and (word_sub_list[4] != [])) or \
      ((word_sub_list[2] == []) and (word_sub_list[3] != [])) or \
      ((word_sub_list[1] == []) and (word_sub_list[2] != [])) or \
      ((word_sub_list[0] == []) and (word_sub_list[1] != []))):
    inout.log_message(['Empty sub-lists between filled sub-lists:', \
                         str(word_sub_list[0]),str(word_sub_list[1]), 
                         str(word_sub_list[2]),str(word_sub_list[3]), 
                         str(word_sub_list[4])],'warn')

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Phase two: Parse sub-lists and assign into four output name lists - - - - -
  #
  givenname_list     = []
  alt_givenname_list = []
  surname_list       = []
  alt_surname_list   = []

  # The second name component is the opposite of the first  - - - - - - - - - -
  #  (givenname <-> surname)
  #
  if (first_name_comp == 'gname'):
    second_name_comp = 'sname'
  else:
    second_name_comp = 'gname'

  # Sub-list 4 not empty: Potential second alternative names  - - - - - - - - -
  #
  if (word_sub_list[4] != []):
    [n, an] = get_name_component(word_sub_list[4], tag_sub_list[4], \
                                 second_name_comp)
    if (second_name_comp == 'sname'):
      alt_surname_list += n + an
    else:
      alt_givenname_list += n + an

  # Sub-list 3 not empty: Potential second alternative names  - - - - - - - - -
  #
  if (word_sub_list[3] != []):
    [n, an] = get_name_component(word_sub_list[3], tag_sub_list[3], \
                                 second_name_comp)
    if (second_name_comp == 'sname'):
      alt_surname_list = n + an + alt_surname_list
    else:
      alt_givenname_list = n + an + alt_givenname_list

  # Sub-list 2 not empty: Potential second primary names  - - - - - - - - - - -
  #
  if (word_sub_list[2] != []):
    [n, an] = get_name_component(word_sub_list[2], tag_sub_list[2], \
                                 second_name_comp)
    if (second_name_comp == 'sname'):
      surname_list += n
      alt_surname_list = an + alt_surname_list
    else:
      givenname_list += n
      alt_givenname_list = an + alt_givenname_list

  # Sub-list 1 not empty: Potential first or second alternative names - - - - -
  # This list can now contain either first or second alternative names.
  # - If the first name component (e.g. givenname) doesn't have an alternative
  #   name, but the second name component (e.g. surname) does, then alternative
  #   surnames are stored in this list.
  # - But if the first name component has an alternative name, then sub-list 2
  #   was not empty, it contained the second name component.
  #
  if (word_sub_list[1] != []):

    if (word_sub_list[2] == []):
      # Sub-list 2 is empty, so this sub-list contain alternative names for the
      # second name component
      #
      [n, an] = get_name_component(word_sub_list[1], tag_sub_list[1], \
                                   second_name_comp)
      if (second_name_comp == 'sname'):
        alt_surname_list = n + an + alt_surname_list
      else:
        alt_givenname_list = n + an + alt_givenname_list

    else:
      # Sub-list 2 not empty, so this sub-list contains alternative names for
      # the first name component
      #
      [n, an] = get_name_component(word_sub_list[1], tag_sub_list[1], \
                                   first_name_comp)
      if (first_name_comp == 'sname'):
        alt_surname_list = n + an + alt_surname_list
      else:
        alt_givenname_list = n + an + alt_givenname_list

  # Sub-list 0 should be non-empty all the time - - - - - - - - - - - - - - - -
  # - If both sub-list 1 and 2 are not empty, then sub-list 0 (this list)
  #   contains only the first primary names.
  # - But if sub-list two is empty, then for both cases sub-list 1 empty or not
  #   this sub-list 0 contains both primary name components.
  #
  if (word_sub_list[0] == []):  # This should never happen!
    inout.log_message('Empty first sub-list: '+str(word_sub_list),'err')

  else:

    if (word_sub_list[1] != []) and (word_sub_list[2] != []):

      # Sub-list 1 and 2 not empty, so this first sub-list contains primary
      # names for the first name component
      #
      [n, an] = get_name_component(word_sub_list[0], tag_sub_list[0], \
                                   first_name_comp)
      if (first_name_comp == 'sname'):
        surname_list += n
        alt_surname_list = an + alt_surname_list
      else:
        givenname_list += n
        alt_givenname_list = an + alt_givenname_list

    else:

      # Sub-list 0 contains both primary names for first and second name
      # components. This is the most common case of names.

      # Loop over all elements in sub-lists 0 - - - - - - - - - - - - - - - - -
      #
      list_len = len(word_sub_list[0])

      # Set current name component where words will be assigned to first
      #
      name_comp_assign = first_name_comp

      gname_ass_count = 0  # Counter for number of words assigned to givenname
      sname_ass_count = 0  # Counter for number of words assigned to surname

      i = 0
      while i < list_len:
        w = word_sub_list[0][i]  # Process this word
        t = tag_sub_list[0][i]   # Corresponding tag

        if (t == 'ST'):  # Process a saint word or name - - - - - - - - - - - -
          if ('_' not in w):  # A saint word only, like saint, brother, etc.
            if (name_comp_assign == 'gname'):
              givenname_list.append(w)
              gname_ass_count += 1
            else:
              surname_list.append(w)
              sname_ass_count += 1
          else:  # A 'saint' plus name sequence, append name to altern. names
            unders_ind = w.find('_')
            if (name_comp_assign == 'gname'):
              givenname_list.append(w)
              alt_givenname_list.append(w[unders_ind+1:])
              gname_ass_count += 1
            else:
              surname_list.append(w)
              alt_surname_list.append(w[unders_ind+1:])
              sname_ass_count += 1

        elif (t in ['RU','CO']):  # A 'rubbuish' word or a comma  - - - - - - -
          pass  # Just pass over it, don't append it to the names list

        elif (t == 'BO'):  # Process 'baby of' and similar sequences  - - - - -
          if (name_comp_assign == 'gname'):
            givenname_list.append(w)
          else:
            surname_list.append(w)

        elif (t == 'PR'):  # Process name prefix  - - - - - - - - - - - - - - -
          if (name_comp_assign == 'gname'):
            givenname_list.append(w)
            gname_ass_count += 1
          else:
            surname_list.append(w)
            sname_ass_count += 1

        elif (t in ['NU','AN']):  # A Number or alphanumeric  - - - - - - - - -
          if (name_comp_assign == 'gname'):
            givenname_list.append(w)
            gname_ass_count += 1
          else:
            surname_list.append(w)
            sname_ass_count += 1
          inout.log_message('Number or alphanumeric word: '+str(word_list), \
                            'warn')

        elif (t == 'HY'):  # Process hyphen - - - - - - - - - - - - - - - - - -
          if (i > 0):
            if (tag_list[i-1] not in ['HY','BO','CO','SP','VB']):
              if (name_comp_assign == 'gname'):
                givenname_list.append(w)
              else:
                surname_list.append(w)
            else:
              inout.log_message('Strange hyphen situation: '+str(word_list), \
                                'warn')

        elif (t == 'II'):  # A single character (initial) - - - - - - - - - - -
            # Append single letter to alternative givennames and switch
            # component to surnames (assuming surname follows after an initial)
            alt_givenname_list += [w]
            name_comp_assign = 'sname'

        # Word must be a name word, tagged with: SN, GF, GM, or UN  - - - - - -
        #
        else:

          # First check if a hyphen is following followed by another name - - -
          # (in which case make sure they get all assigned to the same
          # component)
          #
          if (i < list_len-2):
            if (tag_list[i+1] == 'HY') and \
               (tag_list[i+2] not in ['SP','HY','CO','VB','RU']):
              tmp_w = word_list[i:i+3]
              tmp_len = 3
            else:
              tmp_w = [word_list[i]]
              tmp_len = 1
          else:
            tmp_w = [word_list[i]]
            tmp_len = 1

          # If the current word (or hyphened words) is the last and no
          # givenname has been assigned yet, then assign it to givenname
          #
          if (gname_ass_count == 0) and \
              (((tmp_len == 1) and (i == list_len-1)) or \
               ((tmp_len == 3) and (i == list_len-3))):
            givenname_list += tmp_w
            gname_ass_count += 1

          # If the current word (or hyphened words) is the last and no surname
          # has been assigned yet, then assign it to surname
          #
          elif (sname_ass_count == 0) and \
              (((tmp_len == 1) and (i == list_len-1)) or \
               ((tmp_len == 3) and (i == list_len-3))):
            surname_list += tmp_w
            sname_ass_count += 1

          # If the current name component is givenname and no givenname has
          # been assigned yet
          #
          elif (name_comp_assign == 'gname') and (gname_ass_count == 0):
            givenname_list += tmp_w
            gname_ass_count += tmp_len

          # If the current name component is surname and no surname has
          # been assigned yet
          #
          elif (name_comp_assign == 'sname') and (sname_ass_count == 0):
            surname_list += tmp_w
            sname_ass_count += tmp_len

          # If it no hyphened name and the tag of the current word is
          # givenname (only) then assign to givenname
          #
          elif (tmp_len == 1) and ((t == 'GF') or (t == 'GM')):
            givenname_list += tmp_w
            gname_ass_count += 1

          # If it no hyphened name and the tag of the current word is
          # surname (only) then assign to surname
          #
          elif (tmp_len == 1) and (t == 'SN'):
            surname_list += tmp_w
            sname_ass_count += 1

          # If the current name component is givenname and the current word is
          # (not only) a surname, then assign current name to givenname
          #
          elif (name_comp_assign == 'gname') and (tmp_len == 1) and \
               (t != 'SN'):
            givenname_list += tmp_w
            gname_ass_count += 1

          # If the current name component is surname and the current word is
          # (not only) a givenname, then assign current name to surname
          #
          elif (name_comp_assign == 'sname')  and (tmp_len == 1) and \
                (t != 'GF') and (t != 'GM'):
            surname_list += tmp_w
            sname_ass_count += 1

          else:  # Append to current name component
            if (name_comp_assign == 'gname'):
              givenname_list += tmp_w
              gname_ass_count += 1
            else:
              surname_list += tmp_w
              sname_ass_count += 1

          if (tmp_len == 3):
            i += 2

        i += 1

  # Now do some post-processing steps - - - - - - - - - - - - - - - - - - - - -
  #
  # First, if a hyphened name is in a primary name, or a name containing an
  # underscore, add the components (if they are not name prefixes) into the
  # corresponding alternative name.
  #
  #for n in givenname_list:
  #  if ('-' in n):
  #    tmp_gname_list = n.split('-')
  #  elif ('_' in n):
  #    tmp_gname_list = n.split('_')
  #  else:
  #    tmp_gname_list = []

    # Check each word in the list if they are in a givenname dictionary, not
    # a name prefix, not in ('baby', 'son', 'daughter', 'of') or 'saint'
    #
  #  for gn in tmp_gname_list:  # Now check each component in this list
  #    if (not config.nameprefix_dict.has_key(gn)) and \
  #       (gn not in ['baby','son','daughter','of']):
  #      if (config.givenname_f_dict.has_key(gn)):
  #        alt_givenname_list.append(config.givenname_f_dict[gn])
  #      if (config.givenname_m_dict.has_key(gn)):
  #        alt_givenname_list.append(config.givenname_m_dict[gn])
  #      if (not config.givenname_f_dict.has_key(gn)) and \
  #         (not config.givenname_m_dict.has_key(gn)):
  #        alt_givenname_list.append(gn)

  # Now the same for surnames
  #
  #for n in surname_list:
  #  if ('-' in n):
  #    tmp_sname_list = n.split('-')
  #  elif ('_' in n):
  #    tmp_sname_list = n.split('_')
  #  else:
  #    tmp_sname_list = []

    # Check each word in the list if they are in surname dictionary, not
    # a name prefix, not in ('baby', 'son', 'daughter', 'of') or 'saint'
    # 
  #  for sn in tmp_sname_list:  # Now check each component in this list
  #    if (not config.nameprefix_dict.has_key(sn)) and \
  #       (sn not in ['baby','son','daughter','of']):
  #      if (config.surname_dict.has_key(sn)):
  #        alt_surname_list.append(config.surname_dict[sn])
  #      else:
  #        alt_surname_list.append(sn)

  # Now check for duplicate names and names in alternative names that are 
  # already in primary names
  tmp_givenname_list = []
  for n in givenname_list:
    if (n not in tmp_givenname_list):
      tmp_givenname_list.append(n)
  tmp_alt_givenname_list = []
  for n in alt_givenname_list:
    if (n not in tmp_givenname_list) and (n not in tmp_alt_givenname_list):
      tmp_alt_givenname_list.append(n)

  tmp_surname_list = []
  for n in surname_list:
    if (n not in tmp_surname_list):
      tmp_surname_list.append(n)
  tmp_alt_surname_list = []
  for n in alt_surname_list:
    if (n not in tmp_surname_list) and (n not in tmp_alt_surname_list):
      tmp_alt_surname_list.append(n)

  inout.log_message(['  Extracted output name lists:', \
            '    Givennames:            '+str(tmp_givenname_list), \
            '    Alternative givennames:'+str(tmp_alt_givenname_list), \
            '    Surnames:              '+str(tmp_surname_list), \
            '    Alternative surnames:  '+str(tmp_alt_surname_list)],'v1')

  name_dict = {}
  if (tmp_givenname_list != []):
    name_dict.update({'givenname':tmp_givenname_list})
  if (tmp_alt_givenname_list != []):
    name_dict.update({'alt_givenname':tmp_alt_givenname_list})
  if (tmp_surname_list != []):
    name_dict.update({'surname':tmp_surname_list})
  if (tmp_alt_surname_list != []):
    name_dict.update({'alt_surname':tmp_alt_surname_list})

  return name_dict

# -----------------------------------------------------------------------------

def get_names_hmm(word_list, tag_list):
  """Process the input using HMM to extract name output fields.

  USAGE:
    name_dict = get_names_hmm(word_list, tag_list)

  ARGUMENTS:
    word_list        List of words as produces with clean_tag_names()
    tag_list         Corresponding list of tags as produces with
                     clean_tag_names()

  DESCRIPTION:
    The routine returns a dictionary with the parsed and extracted output
    fields for title and both given- and surnames. A Hidden Markov Model (HMM)
    is used for this task.

    The dictionary returned can contain the following key words:
    - title
    - givenname
    - alt_givenname
    - surname
    - alt_surname
    - name_hmm_proba (the probability returned by the Viterbi algorithm for the
                      most likely HMM state seqence)
  """

  # First, create all permutations of the input tag sequence
  #
  tag_list_seq = mymath.perm_tag_sequence(tag_list)

  msg = ['  Input tag sequence: '+str(tag_list), '  Output tag sequences:']
  for t in tag_list_seq:
    msg.append('    '+str(t))
  inout.log_message(msg,'v2')

  # Now give all tag sequences to the HMM - - - - - - - - - - - - - - - - - - -
  # and keep the one with highest probability
  #
  max_prob = -1.0
  best_obs_seq   = []
  best_tag_list  = []

  for t in tag_list_seq:
    [obs_seq, prob] = config.name_hmm.viterbi(t)
    if (prob > max_prob):
       best_obs_seq  = obs_seq
       best_tag_list = t
       max_prob = prob

    inout.log_message('  Probability '+str(prob)+'  for sequence '+str(t),'v2')

  inout.log_message(['  Best observation sequence: '+str(best_obs_seq),
                     '          with tag sequence: '+str(best_tag_list)],'v2')

  # Now process the observation sequence and add elements into dictionary - - -
  #
  tag_list_len = len(tag_list)
  norm_max_prob = max_prob / float(tag_list_len)  # Normalise max. probability
  name_dict = {'name_hmm_proba':[str(norm_max_prob)]}

  i = 0
  list_len = len(word_list)
  while (i < list_len):  # Loop over words and states
    w = word_list[i]
    s = best_obs_seq[i]

    if (w in ['|',',']):  # Do not output commas and vertical bars  - - - - - -
      pass

    elif (s in ['knwn','andor']):  # Skip over' known as' or 'and/or' - - - - -
      pass  # These should only be separators between names and altern. names

    elif (s == 'titl'):  # Title  - - - - - - - - - - - - - - - - - - - - - - -
      v = name_dict.get('title',[])
      v.append(w)
      name_dict.update({'title':v})

    elif (s == 'baby'):  # 'Baby of' sequence - - - - - - - - - - - - - - - - -

      if (i < list_len-1):
        if (best_obs_seq[i+1] in ['gname1','gname2']):  # Append to givenname
          v = name_dict.get('givenname',[])
          v.append(w)
          name_dict.update({'givenname':v})

        elif (best_obs_seq[i+1] in ['agname1','agname2']):  # Append to agname
          v = name_dict.get('alt_givenname',[])
          v.append(w)
          name_dict.update({'alt_givenname':v})

        elif (best_obs_seq[i+1] in ['sname1','sname2']):  # Append to surname
          v = name_dict.get('surname',[])
          v.append(w)
          name_dict.update({'surname':v})

        elif (best_obs_seq[i+1] in ['asname1','asname2']):  # Append to asname
          v = name_dict.get('alt_surname',[])
          v.append(w)
          name_dict.update({'alt_surname':v})

        else:
          inout.log_message('Strange situation with "baby of": '+ \
                            str(word_list),'warn')
      else:
        inout.log_message('Strange situation with "baby of": '+ \
                          str(word_list),'warn')

    elif (s in ['gname1','gname2']):  # Givennames  - - - - - - - - - - - - - -
      v = name_dict.get('givenname',[])
      v.append(w)
      name_dict.update({'givenname':v})

    elif (s in ['agname1','agname2']):  # Alternative givennames  - - - - - - -
      v = name_dict.get('alt_givenname',[])
      v.append(w)
      name_dict.update({'alt_givenname':v})

    elif (s == 'ghyph'):  # Givenname hyphen  - - - - - - - - - - - - - - - - -
      if (i > 0) and (best_obs_seq[i-1] in ['gname1','gname2']):  # Givenname
        if (name_dict.has_key('givenname')):  # Only append hyphen after gname
          v = name_dict['givenname'] + [w]
          name_dict.update({'givenname':v})
        else:
          inout.log_message('Strange hyphen situation in givenname: '+ \
                            str(word_list),'warn')

      elif (i > 0) and (best_obs_seq[i-1] in ['agname1','agname2']): # agname
        if (name_dict.has_key('alt_givenname')):
          v = name_dict['alt_givenname'] + [w] # Only append hyphen after agn.
        else:
          inout.log_message('Strange hyphen situation in alternative '+ \
                            'givenname: '+str(word_list),'warn')
        name_dict.update({'alt_givenname':v})

      else:
        inout.log_message('Strange givenname hyphen situation: '+ \
                          str(word_list),'warn')

    elif (s in ['sname1','sname2']):  # Surnames  - - - - - - - - - - - - - - -
      v = name_dict.get('surname',[])
      v.append(w)
      name_dict.update({'surname':v})

    elif (s in ['asname1','asname2']):  # Alternative surnames  - - - - - - - -
      v = name_dict.get('alt_surname',[])
      v.append(w)
      name_dict.update({'alt_surname':v})

    elif (s == 'shyph'):  # Surname hyphen  - - - - - - - - - - - - - - - - - -
      if (i > 0) and (best_obs_seq[i-1] in ['sname1','sname2']):  # Surname
        if (name_dict.has_key('surname')):  # Only append hyphen after sname
          v = name_dict['surname'] + [w]
        else:
          inout.log_message('Strange hyphen situation in surname: '+ \
                            str(word_list),'warn')
        name_dict.update({'surname':v})

      elif (i > 0) and (best_obs_seq[i-1] in ['asname1','asname2']):  # asname
        if (name_dict.has_key('alt_surname')):
          v = name_dict['alt_surname'] + [w]  # Only append hyphen after asname
        else:
          inout.log_message('Strange hyphen situation in alternative '+ \
                            'surname: '+str(word_list),'warn')
        name_dict.update({'alt_surname':v})

      else:
        inout.log_message('Strange surname hyphen situation: '+ \
                          str(word_list),'warn')

    elif (s in ['pref1','pref2']):  # Name prefix - - - - - - - - - - - - - - -
      if (i < list_len-1) and (best_obs_seq[i+1] in ['pref1','pref2']):
                                        # Followed by another prefix
        w = w+' '+w[i+1]  # Concatenate
        i = i+1
        s = best_obs_seq[i]

      if (i < list_len-1) and (best_obs_seq[i+1] in ['gname1','gname2']):
        if (name_dict.has_key('givenname')):
          v = name_dict['givenname'] + [w]
          inout.log_message('Strange name prefix situation in givenname: '+ \
                            str(word_list),'warn')
        else:
          v = [w]
        name_dict.update({'givenname':v})

      elif (i < list_len-1) and (best_obs_seq[i+1] in ['agname1','agname2']):
        if (name_dict.has_key('alt_givenname')):
          v = name_dict['alt_givenname'] + [w]
          inout.log_message('Strange name prefix situation in alternative '+ \
                            'givenname: '+str(word_list),'warn')
        else:
          v = [w]
        name_dict.update({'alt_givenname':v})

      elif (i < list_len-1) and (best_obs_seq[i+1] in ['sname1','sname2']):
        if (name_dict.has_key('surname')):
          v = name_dict['surname'] + [w]
          inout.log_message('Strange name prefix situation in surname: '+ \
                            str(word_list),'warn')
        else:
          v = [w]
        name_dict.update({'surname':v})

      elif (i < list_len-1) and (best_obs_seq[i+1] in ['asname1','asname2']):
        if (name_dict.has_key('alt_surname')):
          v = name_dict['alt_surname'] + [w]
          inout.log_message('Strange name prefix situation in alternative '+ \
                            'surname: '+str(word_list),'warn')
        else:
          v = [w]
        name_dict.update({'alt_surname':v})

      else:
        inout.log_message('Strange name prefix situation: '+ \
                          str(word_list),'warn')

    else:  # Should never happen
      msg = ['This should never happen!', '  Tag: '+str(s), '  Word: '+w, \
             '  Word list: '+str(word_list), \
             '  tag list:  '+str(tag_list)]
      inout.log_message(msg,'warn')

    i +=1

  # Finally do some tests on the output fields  - - - - - - - - - - - - - - - -
  #
  name_items = name_dict.items()

  # Check if a value list has more than three elements, if so print out
  #
  for i in name_items:
    if (len(i[1]) > 3):
      inout.log_message('Name output field '+ i[0]+ \
              ' contains more than three elements: '+str(i[1]),'warn')

  # Check if a name component is not allocated but its alternative name is
  #
  if (name_dict.get('givenname',[]) == []) and \
     (name_dict.get('alt_givenname',[]) != []):
    inout.log_message('No givenname but an alternative givenname: '+ \
                      str(name_dict['givenname']),'warn')

  if (name_dict.get('surname',[]) == []) and \
     (name_dict.get('alt_surname',[]) != []):
    inout.log_message('No surname but an alternative surname: '+ \
                       str(name_dict['surname']),'warn')

  # Check if a givenname is given but no surname  - - - - - - - - - - - - - - -
  #
  if (name_dict.get('givenname',[]) != []) and \
     (name_dict.get('surname',[]) == []):
    inout.log_message('No surname but one or more givennames: '+ \
                       str(name_dict['givenname']),'warn')

  # Check if title words are known from lookup-table  - - - - - - - - - - - - -
  #
  tmp_title = name_dict.get('title',[])
  for t in tmp_title:
    if (not config.name_lookup_dict.has_key((t,))):
      inout.log_message('Title word: '+t+' is not a known title','warn')

  return name_dict

# -----------------------------------------------------------------------------

def test():
  """Simple test routine with example inputs.
  """

  # Records with givenname and surname in one string
  #
  test_records_one = [ \
    'dr. steve roberts                                 ', \
    'miss monica moneypenny                            ', \
    'john (parrot) cleese                              ', \
    'john {parrot} cleese                              ', \
    'john "parrot" cleese                              ', \
    'mr. dr. christen, peter                           ', \
    'peter, dr. christen                               ', \
    'mr. doctor christen, mr peter                     ', \
    'peter, dr. christen, Mr. PhD                      ', \
    'doctor john                                       ', \
    'mr. rev rees                                      ', \
    'dr jean-pierre (alfred) echo                      ', \
    'paola del rio, miss                               ', \
    'paola del rio, miss also known as the queen       ', \
    'paola del-rio (queen)                             ', \
    'henrik van de felde                               ', \
    'peter a.k.a. pete christen, mr., PHD              ', \
    'peter k.a. pete christen, mr., PhD                ', \
    'peter known as pete christen-miller (jimmy)       ', \
    'henrik van de felde also known as henriko         ', \
    'henrik van-de-felde also known as henriko         ', \
    'henrik van de-felde also known as henriko         ', \
    'henrik van-de felde also known as henriko         ', \
    'jo ling                                           ', \
    'jo-ling hu                                        ', \
    'del-ponte-monte ricotta                           ', \
    'baby of gayle west (guapa)                        ', \
    'babe molly miller-richie                          ', \
    'babe of veronica-louise (also known as sue) west  ', \
    'b/o veronica (lousie or lou) miller, miss         ', \
    'peter a.k.a. pete also known as sharky christen   ', \
    'peter j. (sharky) known as pete christen, Phd     ', \
    'peter del rio a.k.a. pet or josef christen-jutz   ', \
    'peter (pete) christen-del rio a.k.a. jumbo        ', \
    'del (van a.k.a. peter) miller, phd                ', \
    'peter-josef christen-jutz                         ', \
    'peter-josef christen-jutz, baby                   ', \
    'saint peter christen                              ', \
    'saint peter-christen josef                        ', \
    'peter st christen (josef)                         ', \
    'peter-christen                                    ', \
    'peter-christen-jutz                               ', \
    'peter-christen jutz                               ', \
    'peter christen-jutz                               ', \
    'AASIA FARHAT UL AIN YEUNG AKA YU                  ', \
    'TE ARAKOREKO AILLIE KNOWN AS AILEEN GIVSON        ', \
    'FRANCIS KN AS BOB GIVSON (NEE BATES)              ', \
    'KNOWN AS ROY YEUNG AKA YU                         ', \
    'van-den pete                                      ', \
    'vanden peter                                      ', \
   ]

  # Records with givenname and surname separated
  #
  test_records_two = [ \
    ['peter                              ','christen                     '], \
    ['ole m.                             ','nielsen, phd                 '], \
    ['peter josef                        ','christen                     '], \
    ['doctor churches,                   ','tim                          '], \
    ['AASIA FARHAT UL AIN                ','GIVSON (NEE BATES)           '], \
    ['TE ARAKOREKO AILLIE KNOWN AS AILEEN','YEUNG AKA YU                 '], \
    ['FRANCIS KN AS BOB                  ','GIVSON NEE                   '], \
    ['KNOWN AS ROY                       ','nee bates                    '], \
    ['jean-pierre (alfredo)              ','miller-meyer                 '], \
    ['son of peter also known as piedro  ','christen-jutz                '], \
    ['paola (queen)                      ','del-rio                      '], \
    ['paola "queen"                      ','del rio                      '], \
    ['sonof gayle                      ','west-miller (nee mueller)    '], \
    ['babe molly                         ','van der miller-richie        '], \
    ['babe-of veronica-louise (also known as sue)','best west            '], \
    ['b/o veronica (lousie or lou)       ','miller, miss                 '], \
    ['peter a.k.a. pete also known as sharky','del jean-pierre           '], \
    ['peter j. (sharky) known as pete','christen, mister Phd             '], \
    ['peter del rio a.k.a. pet or josef','christen-jutz                  '], \
    ['peter (pete)','christen-del rio a.k.a. jumbo                       '], \
    ['del (van a.k.a. peter)             ','miller, phd                  '], \
    ['st john                            ','st-john a.k.a. miller/meyer  '], \
    ['st (john)                          ','van der (known as peter) mike'], \
    ['(ko) lee                           ','st clair-telferd             '], \
    ['allan-(squire)                     ','watson (harrison)            '], \
    ['baby of mary miller                ','baby-of-mary miller          '], \
    ['baby of-mary miller                ','babyof-mary-miller           '], \
    ['baby-mary miller                   ','babyof mary-miller           '], \
    ['known mary miller                  ','knownas-mary-miller          '], \
   ]

  for relem in test_records_two:

    for r in relem:
      print
      print 'INPUT: "%s"'%(r)
      clean_str = clean_name_component(r.lower())
      print 'CLEAN: "%s"'%(clean_str)
      [tmp_list,tag_list] = tag_name_component(clean_str)
      print '  Words:',tmp_list 
      print '  Tags: ',tag_list 

      g = get_gender_guess(tmp_list,tag_list)
      if (g != ''):
        print '  Gender guess:', g
      [tmp_list,tag_list,title_list] = get_title(tmp_list,tag_list)
      if title_list != []:
        print '  Title(s):', title_list
 
      #[l1, l2] = get_name_component(tmp_list, tag_list, 'gname')
      #if (l1 != []):
      #  print '  G: Names:', l1 
      #if (l2 != []):
      #  print '  G: Altnames:', l2
      #[l1, l2] = get_name_component(tmp_list, tag_list, 'sname')
      #if (l1 != []):
      #  print '  S: Names:', l1 
      #if (l2 != []):
      #  print '  S: Altnames:', l2
      print

  print
  print
  print
  print
  print
  print

  for r in test_records_one:
    print
    print 'INPUT: "%s"'%(r)
    clean_str = clean_name_component(r.lower())
    print 'CLEAN: "%s"'%(clean_str)
    [tmp_list,tag_list] = tag_name_component(clean_str)
    print '  Words:',tmp_list 
    print '  Tags: ',tag_list 

    g = get_gender_guess(tmp_list,tag_list)
    if (g != ''):
      print '  Gender guess:', g
    [tmp_list,tag_list,title_list] = get_title(tmp_list,tag_list)
    if title_list != []:
      print '  Title(s):', title_list

#      [l1,l2,l3,l4] = get_names_rules(tmp_list,tag_list, 'gname')
#      if (l1 != []):
#        print '  G: Names:', l1 
#      if (l2 != []):
#        print '  G: Altnames:', l2
#      if (l3 != []):
#        print '  S: Names:', l3 
#      if (l4 != []):
#        print '  S: Altnames:', l4

    line_dict = get_names_hmm(tmp_list,tag_list)
    print line_dict

    print

# =============================================================================
