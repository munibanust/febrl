# =============================================================================
# address.py - Routines for address cleaning and standardisation.
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
# The Original Software is "address.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module address.py - Routines for address cleaning and standardisation.

   PUBLIC FUNCTIONS:
     tag_address_component    Tag an address component input string and make a
                              list
     get_address_hmm          Process the input word and tag lists using a
                              Hidden Markov Model (HMM) to extract address
                              output fields

   See doc strings of individual functions for detailed documentation.

   TODO:
   - PC 10/12/2002: Try to correct warnings in get_address_hmm
"""

# =============================================================================
# Imports go here

import string
import mymath

# =============================================================================

def tag_address_component(address_str, tag_lookup_table, record_id):
  """Tag an address input component string and make a list.

  USAGE:
    [word_list, tag_list] = tag_address_component(address_str,
                                                  tag_lookup_table)
  ARGUMENTS:
    address_str        A string containing the address component
    tag_lookup_table  A tagging look-up table as defined in 'lookup.py'
    record_id         A string identifying the current record

  DESCRIPTION:
    This routines cleans the input string and extracts words, numbers and
    separators into a list. Each element of this list is assigned one or more
    tags. A 'greedy tagger' is applied, which cheques sequences of list
    elements in the given lookup table (longer sequences first) and replaces
    them with the string and tag from the lookup-table if found.

    The routine returns two lists: words and their tags
  """

  # First, split input string into elements at spaces - - - - - - - - - - - - -
  #
  org_list = address_str.split()  # The original list from the input string

  tag_list  = []  # The initially empty list of tags
  word_list = []  # The initially empty list of words

  while (org_list != []):  # As long as not all elements have been processed
    tmp_list = org_list[:tag_lookup_table.max_key_length]
                                                  # Extract longest sub-list
    tmp_val = []  # Start with empty value
    tmp_key = tuple(tmp_list)

    while (tmp_key != ()):  # As long as key not empty and not found in lookup
      if (tag_lookup_table.has_key(tmp_key)):
        tmp_val = tag_lookup_table[tmp_key]
        break
      tmp_key = tmp_key[:-1]  # Remove last element in key

    if (tmp_val != []):  # A value has been found in the dictionary
      tmp_len = len(tmp_key)  # Length of found sequence

      if (tmp_val[0] != ''):  # It's not an empty value
        word_list.append(tmp_val[0])  # Append corrected word (or sequence)
        tag_list.append(tmp_val[1])   # Append tag or tags

    else:  # No value has been found in the lookup dictionary, try other tags

      tmp_val = org_list[0]  # Value is first element in the original list
      tmp_len = 1

      if (tmp_val.isdigit()):  # Element is a number
        word_list.append(tmp_val)
        if (len(tmp_val) == 4):
          tag_list.append('N4')
        else:
          tag_list.append('NU')

      elif (not tmp_val.isalpha()) and tmp_val.isalnum():  # Alpha-numeric
        word_list.append(tmp_val)
        tag_list.append('AN')

      elif (tmp_val == '-'):  # Element is a hyphen
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

  # A log message for high volume log output (level 3)  - - - - - - - - - - - -
  #
  print '3:%s  Address string "%s"' % (record_id, address_str) 
  print '3:%s    Split into word list: %s' % (record_id, word_list)
  print '3:%s            and tag list: %s' % (record_id, tag_list)

  return [word_list, tag_list]

# =============================================================================

def get_address_hmm(word_list, tag_list, address_hmm, tag_lookup_table,
                    record_id, fields_str):
  """Process the input using a HMM to extract address output fields.

  USAGE:
    address_dict = get_address_hmm(word_list, tag_list, address_hmm,
                                   tag_lookup_table)

  ARGUMENTS:
    word_list         List of words as produces with tag_address_component()
    tag_list          Corresponding list of tags as produces with
                      tag_address_component()
    address_hmm       A reference to the address hidden Markov model
    tag_lookup_table  A tagging look-up table as defined in 'lookup.py'
    record_id         A string identifying the current record
    fields_str        A string representation of the input fields

  DESCRIPTION:
    The routine returns a dictionary with the parsed and extracted output
    fields for the address component. A Hidden Markov Model (HMM) is used for
    this task.

    The dictionary returned can contain the following key words:
      wayfare_number
      wayfare_name
      wayfare_qualifier
      wayfare_type
      unit_number
      unit_type
      property_name
      institution_name
      institution_type    
      postaddress_number
      postaddress_type
      locality_name
      locality_qualifier
      postcode
      territory
      country
      address_hmm_prob (the probability returned by the Viterbi algorithm for
                        the most likely HMM state seqence)
  """

  # First, create all permutations of the input tag sequence
  #
  tag_list_seq = mymath.perm_tag_sequence(tag_list)

  # Now give all tag sequences to the HMM - - - - - - - - - - - - - - - - - - -
  # and keep the one with highest probability
  #
  max_prob = -1.0
  best_obs_seq   = []
  best_tag_list  = []

  for t in tag_list_seq:
    [obs_seq, prob] = address_hmm.viterbi(t)
    if (prob > max_prob):
       best_obs_seq  = obs_seq
       best_tag_list = t
       max_prob = prob

    print '3:%s  Sequence: %s has Viterbi probability: %f' % \
          (record_id, str(t), prob)

  print '2:%s  Best observation sequence: %s with tag sequence: %s' % \
        (record_id, str(best_obs_seq), str(best_tag_list))

  # Now process the observation sequence and add elements into dictionary - - -
  #
  if (len(tag_list) != len(word_list)):
    print 'error:%s Length of word list and tag list differs: %s, %s%s' % \
          (record_id, str(word_list), str(tag_list), fields_str)
    raise Exception

  list_len = len(tag_list)

  if (list_len == 0):
    print 'warning:%s Empty tag list returned from HMM %s' % \
          (record_id, fields_str)
    return {}  # Return an empty dictionary if not output fields given

  # norm_max_prob = max_prob / float(list_len)  # Normalise max. probability
  address_dict = {'address_hmm_prob':[str(max_prob)]}

  for i in range(list_len):  # Loop over words and states
    w = word_list[i]
    s = best_obs_seq[i]

    #  Do not output commas, vertical bars and hyphens  - - - - - - - - - - - -
    #
    if (w in ['|', ',', '-','/']):
      pass

    elif (s == 'wfnu'):  # Wayfare number - - - - - - - - - - - - - - - - - - -
      v = address_dict.get('wayfare_number',[])
      v.append(w)
      address_dict.update({'wayfare_number':v})

    elif (s in ['wfna1','wfna2','wfna3']):  # Wayfare name  - - - - - - - - - -
      v = address_dict.get('wayfare_name',[])
      v.append(w)
      address_dict.update({'wayfare_name':v})

    elif (s == 'wfql'):  # Wayfare qualifier  - - - - - - - - - - - - - - - - -
      v = address_dict.get('wayfare_qualifier',[])
      v.append(w)
      address_dict.update({'wayfare_qualifier':v})

    elif (s == 'wfty'):  # Wayfare type - - - - - - - - - - - - - - - - - - - -
      v = address_dict.get('wayfare_type',[])
      v.append(w)
      address_dict.update({'wayfare_type':v})

    elif (s == 'unnu'):  # Unit number  - - - - - - - - - - - - - - - - - - - -
      v = address_dict.get('unit_number',[])
      v.append(w)
      address_dict.update({'unit_number':v})

    elif (s == 'unty'):  # Unit type  - - - - - - - - - - - - - - - - - - - - -
      v = address_dict.get('unit_type',[])
      v.append(w)
      address_dict.update({'unit_type':v})

    elif (s in ['prna1','prna2']):  # Property name - - - - - - - - - - - - - -
      v = address_dict.get('property_name',[])
      v.append(w)
      address_dict.update({'property_name':v})

    elif (s in ['inna1','inna2']):  # Institution name  - - - - - - - - - - - -
      v = address_dict.get('institution_name',[])
      v.append(w)
      address_dict.update({'institution_name':v})

    elif (s == 'inty'):  # Institution type - - - - - - - - - - - - - - - - - -
      v = address_dict.get('institution_type',[])
      v.append(w)
      address_dict.update({'institution_type':v})

    elif (s == 'panu'):  # Postal address number  - - - - - - - - - - - - - - -
      v = address_dict.get('postaddress_number',[])
      v.append(w)
      address_dict.update({'postaddress_number':v})

    elif (s == 'paty'):  # Postal address type  - - - - - - - - - - - - - - - -
      v = address_dict.get('postaddress_type',[])
      v.append(w)
      address_dict.update({'postaddress_type':v})

    elif (s in ['loc1','loc2']):  # Locality name - - - - - - - - - - - - - - -
      v = address_dict.get('locality_name',[])
      v.append(w)
      address_dict.update({'locality_name':v})

    elif (s == 'locql'):  # Locality qualifier  - - - - - - - - - - - - - - - -
      v = address_dict.get('locality_qualifier',[])
      v.append(w)
      address_dict.update({'locality_qualifier':v})

    elif (s == 'pc'):  # Postcode - - - - - - - - - - - - - - - - - - - - - - -
      v = address_dict.get('postcode',[])
      v.append(w)
      address_dict.update({'postcode':v})

    elif (s in ['ter1','ter2']):  # Territory - - - - - - - - - - - - - - - - -
      v = address_dict.get('territory',[])
      v.append(w)
      address_dict.update({'territory':v})

    elif (s in ['cntr1','cntr2']):  # Country - - - - - - - - - - - - - - - - -
      v = address_dict.get('country',[])
      v.append(w)
      address_dict.update({'country':v})

    else:  # Should never happen
      print 'warning:%s This should never happen! ' % (record_id) + \
            ' Tag: %s, word: %s, word list: %s, tag list: %s%s' % \
            (str(s), w, str(word_list), str(tag_list),fields_str)

  # Check if concatenated locality and territory words are in lookup-table  - -
  #
  if (address_dict.has_key('locality_name')):
    loc = address_dict['locality_name']
    if (len(loc) > 1):  # Locality contains more than one word
      loc_tuple = tuple(loc)  # Make it a tuple
      if (tag_lookup_table.has_key(loc_tuple)):
         new_loc = tag_lookup_table[loc_tuple][0]
         address_dict.update({'locality_name':[new_loc]})

  if (address_dict.has_key('territory')):
    terr = address_dict['territory']
    if (len(terr) > 1):  # Territory contains more than one word
      terr_tuple = tuple(terr)  # Make it a tuple
      if (tag_lookup_table.has_key(terr_tuple)):
         new_terr = tag_lookup_table[terr_tuple][0]
         address_dict.update({'territory':[new_terr]})

  if (address_dict.has_key('country')):
    cntr = address_dict['country']
    if (len(cntr) > 1):  # Country contains more than one word
      cntr_tuple = tuple(cntr)  # Make it a tuple
      if (tag_lookup_table.has_key(cntr_tuple)):
         new_cntr = tag_lookup_table[cntr_tuple][0]
         address_dict.update({'country':[new_cntr]})

  # Finally do some tests on the output fields  - - - - - - - - - - - - - - - -
  #
  address_items = address_dict.items()

  # Check if a value list has more than three elements, if so print out
  #
  for i in address_items:
    if (len(i[1]) > 3):
      print 'warning:%s Output field "%s" contains' % (record_id, str(i[0]))+ \
            ' more than three elements: %s%s' % (str(i[1]), fields_str)

  # Check if 'number' elements only contain (alpha-) numerical values - - - - -
  # and also check how many numbers in an element
  #
  if (address_dict.has_key('wayfare_number')): # Check how many wayfare numbers
    v = address_dict['wayfare_number']
    if (len(v) > 2):
      print 'warning:%s More than two wayfare numbers: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        print 'warning:%s Wayfare number contains no ' % (record_id) + \
              'digits: %s%s' % (str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('unit_number')):  # Check how many unit numbers
    v = address_dict['unit_number']
    if (len(v) > 1):
      print 'warning:%s More than one unit number: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        print 'warning:%s Unit number contains no ' % (record_id) + \
              'digits: %s%s' % (str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('postaddress_number')): # Check postaddress numbers
    v = address_dict['postaddress_number']
    if (len(v) > 1):
      print 'warning:%s More than one post-address number: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        print 'warning:%s Post-address number contains no ' % (record_id) + \
              'digits: %s%s' % (str(v), fields_str)
        break  # Exit for loop

  # Check if 'type' elements contain one word only  - - - - - - - - - - - - - -
  # if it's a known type word
  #
  if (address_dict.has_key('wayfare_type')):  # Check wayfare type
    v = address_dict['wayfare_type']
    if (len(v) > 1):
      print 'warning:%s More than one wayfare type: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not tag_lookup_table.has_key((i))) or \
         (tag_lookup_table.has_key((i)) and 
          (tag_lookup_table[(i)][1].find('WT') < 0)):
        print 'warning:%s Wayfare type word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('unit_type')):  # Check unit type
    v = address_dict['unit_type']
    if (len(v) > 1):
      print 'warning:%s More than one unit type: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not tag_lookup_table.has_key((i))) or \
         (tag_lookup_table.has_key((i)) and \
          (tag_lookup_table[(i)][1].find('UT') < 0)):
        print 'warning:%s Unit type word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('institution_type')):  # Check institution type
    v = address_dict['institution_type']
    if (len(v) > 1):
      print 'warning:%s More than one institution type: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not tag_lookup_table.has_key((i))) or \
         (tag_lookup_table.has_key((i)) and \
          (tag_lookup_table[(i)][1].find('IT') < 0)):
        print 'warning:%s Institution type word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('postaddress_type')):  # Check postaddress type
    v = address_dict['postaddress_type']
    if (len(v) > 2):
      print 'warning:%s More than two post-address type: %s%s' % \
            (record_id, str(v), fields_str)
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not tag_lookup_table.has_key((i))) or \
         (tag_lookup_table.has_key((i)) and \
          (tag_lookup_table[(i)][1].find('PA') < 0)):
        print 'warning:%s Post-address type word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  # Check if 'qualifier' elements only contain known qualifier words  - - - - -
  #
  if (address_dict.has_key('wayfare_qualifier')):  # Check wayfare qualifier
    v = address_dict['wayfare_qualifier']
    for i in v:
      if (not tag_lookup_table.has_key((i,))) or \
         (tag_lookup_table.has_key((i,)) and \
          (tag_lookup_table[(i,)][1].find('LQ') < 0)):
        print 'warning:%s Wayfare qualifier word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  if (address_dict.has_key('locality_qualifier')):  # Check locality qualifier
    v = address_dict['locality_qualifier']
    for i in v:
      if (not tag_lookup_table.has_key((i,))) or \
         (tag_lookup_table.has_key((i,)) and \
          (tag_lookup_table[(i,)][1].find('LQ') < 0)):
        print 'warning:%s Locality qualifier word is not known: %s%s' % \
              (record_id, str(v), fields_str)
        break  # Exit for loop

  return address_dict

# =============================================================================
