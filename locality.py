# =============================================================================
# locality.py - Routines for locality standardisation and linkage.
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
# The Original Software is "locality.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module locality.py - Routines for locality standardisation and linkage.

   PUBLIC FUNCTIONS:
     clean_geoloc_component  Clean a geocode or locality input component string
     tag_geoloc_component    Tag a geocode or locality component input string
                             and make a list
     get_geoloc_hmm          Process the input word and tag lists using a
                             Hidden Markov Model (HMM) to extract geocode and
                             locality output fields
     test                    Simple test routine with example inputs

   See doc strings of individual functions for detailed documentation.

   See also the relevant section in the config.py module.
"""

# -----------------------------------------------------------------------------

import string
import types

import config
import inout
import mymath

# -----------------------------------------------------------------------------

def clean_geoloc_component(geoloc_str):
  """Clean a geocode or locality component input string.

  USAGE:
    cleaned_str = clean_geoloc_component(geoloc_str)

  ARGUMENTS:
    geoloc_str  A string containing the geocode and/or locality components
                (or parts of it/them)

  DESCRIPTION:
    This routine cleans the input string by using the 'geoloc_corr_list'. It
    also strips off all leading and trailing spaces. A cleaned string is
    returned.
  """
  # First add a trailing and leading space  - - - - - - - - - - - - - - - - - -
  # (this is to make sure replacement strings do match at beginning and end)
  #
  geoloc_str = ' '+geoloc_str+' '

  # Check for strings from the geocode/locality correction-list - - - - - - - -
  #
  for (org,repl) in config.geoloc_corr_list:
    geoloc_str = geoloc_str.replace(org,repl)

  # Make sure commas are separated from words so they become list elements  - -
  #
  geoloc_str = geoloc_str.replace(',', ' , ')

  return geoloc_str.strip()

# -----------------------------------------------------------------------------

def tag_geoloc_component(geoloc_str):
  """Tag a geocode locality input component string and make a list.

  USAGE:
    [word_list, tag_list] = tag_geoloc_component(loc_str)

  ARGUMENTS:
    geoloc_str  A string containing the geocode and/or locality components

  DESCRIPTION:
    This routines cleans the input string and extracts words, numbers and
    separators into a list. Each element of this list is assigned one or more
    tags. A 'greedy tagger' is applied, which cheques sequences of list
    elements in the name lookup table (longer sequences first) and replaces
    them with the string and tag from the lookup-table if found.

    The routine returns two lists: words and their tags
  """

  # First, split input string into elements at spaces - - - - - - - - - - - - -
  #
  org_list = geoloc_str.split()  # The original list from the input string
  inout.log_message('  Initial word list: '+str(org_list),'v2')

  tag_list  = []  # The initially empty list of tags
  word_list = []  # The initially empty list of words

  while (org_list != []):  # As long as not all elements have been processed
    tmp_list = org_list[:config.geoloc_dict_seq_len] # Extract longest sub-list
    tmp_val = []  # Start with empty value
    tmp_key = tuple(tmp_list)

    while (tmp_key != ()):  # As long as key not empty and not found in lookup
      if (config.geoloc_lookup_dict.has_key(tmp_key)):
        tmp_val = config.geoloc_lookup_dict[tmp_key]
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

  return [word_list, tag_list]

# -----------------------------------------------------------------------------

def get_geoloc_hmm(word_list, tag_list):
  """Process input using a HMM to extract geocode and locality output fields.

  USAGE:
    geoloc_dict = get_geoloc_hmm(word_list, tag_list)

  ARGUMENTS:
    word_list  List of words as produces with clean_tag_locality()
    tag_list   Corresponding list of tags as produces with
               clean_tag_locality()

  DESCRIPTION:
    The routine returns a dictionary with the parsed and extracted output
    fields for both the locality and geocode components. A Hidden Markov Model
    (HMM) is used for this task.

    The dictionary returned can contain the following key words:
    - wayfare_number
    - wayfare_name
    - wayfare_qualifier
    - wayfare_type
    - unit_number
    - unit_type
    - property_name
    - institution_name
    - institution_type    
    - postaddress_number
    - postaddress_type
    - locality_name
    - locality_qualifier
    - postcode
    - territory
    - country
    - geoloc_hmm_proba (the probability returned by the Viterbi algorithm for
                        the most likely HMM state seqence)
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
    [obs_seq, prob] = config.geoloc_hmm.viterbi(t)
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
  geoloc_dict = {'geoloc_hmm_proba':[str(norm_max_prob)]}

  list_len = len(word_list)
  for i in range(list_len):  # Loop over words and states
    w = word_list[i]
    s = best_obs_seq[i]

    #  Do not output commas, vertical bars and hyphens  - - - - - - - - - - - -
    #
    if (w in ['|', ',', '-','/']):
      pass

    elif (s == 'wfnu'):  # Wayfare number - - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('wayfare_number',[])
      v.append(w)
      geoloc_dict.update({'wayfare_number':v})

    elif (s in ['wfna1','wfna2','wfna3']):  # Wayfare name  - - - - - - - - - -
      v = geoloc_dict.get('wayfare_name',[])
      v.append(w)
      geoloc_dict.update({'wayfare_name':v})

    elif (s == 'wfql'):  # Wayfare qualifier  - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('wayfare_qualifier',[])
      v.append(w)
      geoloc_dict.update({'wayfare_qualifier':v})

    elif (s == 'wfty'):  # Wayfare type - - - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('wayfare_type',[])
      v.append(w)
      geoloc_dict.update({'wayfare_type':v})

    elif (s == 'unnu'):  # Unit number  - - - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('unit_number',[])
      v.append(w)
      geoloc_dict.update({'unit_number':v})

    elif (s == 'unty'):  # Unit type  - - - - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('unit_type',[])
      v.append(w)
      geoloc_dict.update({'unit_type':v})

    elif (s in ['prna1','prna2']):  # Property name - - - - - - - - - - - - - -
      v = geoloc_dict.get('property_name',[])
      v.append(w)
      geoloc_dict.update({'property_name':v})

    elif (s in ['inna1','inna2']):  # Institution name  - - - - - - - - - - - -
      v = geoloc_dict.get('institution_name',[])
      v.append(w)
      geoloc_dict.update({'institution_name':v})

    elif (s == 'inty'):  # Institution type - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('institution_type',[])
      v.append(w)
      geoloc_dict.update({'institution_type':v})

    elif (s == 'panu'):  # Postal address number  - - - - - - - - - - - - - - -
      v = geoloc_dict.get('postaddress_number',[])
      v.append(w)
      geoloc_dict.update({'postaddress_number':v})

    elif (s == 'paty'):  # Postal address type  - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('postaddress_type',[])
      v.append(w)
      geoloc_dict.update({'postaddress_type':v})

    elif (s in ['loc1','loc2']):  # Locality name - - - - - - - - - - - - - - -
      v = geoloc_dict.get('locality_name',[])
      v.append(w)
      geoloc_dict.update({'locality_name':v})

    elif (s == 'locql'):  # Locality qualifier  - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('locality_qualifier',[])
      v.append(w)
      geoloc_dict.update({'locality_qualifier':v})

    elif (s == 'pc'):  # Postcode - - - - - - - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('postcode',[])
      v.append(w)
      geoloc_dict.update({'postcode':v})

    elif (s in ['ter1','ter2']):  # Territory - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('territory',[])
      v.append(w)
      geoloc_dict.update({'territory':v})

    elif (s in ['cntr1','cntr2']):  # Country - - - - - - - - - - - - - - - - -
      v = geoloc_dict.get('country',[])
      v.append(w)
      geoloc_dict.update({'country':v})

    else:  # Should never happen
      msg = ['This should never happen!', '  Tag: '+str(s), '  Word: '+w, \
             '  Word list: '+str(word_list), \
             '  tag list:  '+str(tag_list)]
      inout.log_message(msg,'warn')

  # Check if concatenated locality and territory words are in lookup-table  - -
  #
  if (geoloc_dict.has_key('locality_name')):
    loc = geoloc_dict['locality_name']
    if (len(loc) > 1):  # Locality contains more than one word
      loc_tuple = tuple(loc)  # Make it a tuple
      if (config.geoloc_lookup_dict.has_key(loc_tuple)):
         new_loc = config.geoloc_lookup_dict[loc_tuple][0]
         geoloc_dict.update({'locality_name':[new_loc]})

  if (geoloc_dict.has_key('territory')):
    terr = geoloc_dict['territory']
    if (len(terr) > 1):  # Territory contains more than one word
      terr_tuple = tuple(terr)  # Make it a tuple
      if (config.geoloc_lookup_dict.has_key(terr_tuple)):
         new_terr = config.geoloc_lookup_dict[terr_tuple][0]
         geoloc_dict.update({'territory':[new_terr]})

  if (geoloc_dict.has_key('country')):
    cntr = geoloc_dict['country']
    if (len(cntr) > 1):  # Country contains more than one word
      cntr_tuple = tuple(cntr)  # Make it a tuple
      if (config.geoloc_lookup_dict.has_key(cntr_tuple)):
         new_cntr = config.geoloc_lookup_dict[cntr_tuple][0]
         geoloc_dict.update({'country':[new_cntr]})

  # Finally do some tests on the output fields  - - - - - - - - - - - - - - - -
  #
  geoloc_items = geoloc_dict.items()

  # Check if a value list has more than three elements, if so print out
  #
  for i in geoloc_items:
    if (len(i[1]) > 3):
      inout.log_message('Geocode/locality output field '+ str(i[0])+ \
              ' contains more than three elements: '+str(i[1]),'warn')

  # Check if 'number' elements only contain (alpha-) numerical values - - - - -
  # and also check how many numbers in an element
  #
  if (geoloc_dict.has_key('wayfare_number')):  # Check how many wayfare numbers
    v = geoloc_dict['wayfare_number']
    if (len(v) > 2):
      inout.log_message('More than two wayfare numbers: '+str(v),'warn')
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        inout.log_message('Wayfare number element contains no digits: '+ \
                          str(v),'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('unit_number')):  # Check how many unit numbers
    v = geoloc_dict['unit_number']
    if (len(v) > 1):
      inout.log_message('More than one unit numbers: '+str(v),'warn')
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        inout.log_message('Unit number element contains no digits: '+str(v),\
                          'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('postaddress_number')): # Check postaddress numbers
    v = geoloc_dict['postaddress_number']
    if (len(v) > 1):
      inout.log_message('More than one postaddress numbers: '+str(v),'warn')
    for i in v:
      if (i.isalpha()):  # Element contains only letters
        inout.log_message('Postaddress number element contains no digits: '+ \
                          str(v),'warn')
        break  # Exit for loop

  # Check if 'type' elements contain one word only  - - - - - - - - - - - - - -
  # if it's a known type word
  #
  if (geoloc_dict.has_key('wayfare_type')):  # Check wayfare type
    v = geoloc_dict['wayfare_type']
    if (len(v) > 1):
      inout.log_message('More than one wayfare type: '+str(v),'warn')
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not config.geoloc_lookup_dict.has_key((i))) or \
         (config.geoloc_lookup_dict.has_key((i)) and \
          (config.geoloc_lookup_dict[(i)][1].find('WT') < 0)):
        inout.log_message('Wayfare type word is not known: '+str(v),'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('unit_type')):  # Check unit type
    v = geoloc_dict['unit_type']
    if (len(v) > 1):
      inout.log_message('More than one unit type: '+str(v),'warn')
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not config.geoloc_lookup_dict.has_key((i))) or \
         (config.geoloc_lookup_dict.has_key((i)) and \
          (config.geoloc_lookup_dict[(i)][1].find('UT') < 0)):
        inout.log_message('Unit type word is not known: '+str(v),'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('institution_type')):  # Check institution type
    v = geoloc_dict['institution_type']
    if (len(v) > 1):
      inout.log_message('More than one institution type: '+str(v),'warn')
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not config.geoloc_lookup_dict.has_key((i))) or \
         (config.geoloc_lookup_dict.has_key((i)) and \
          (config.geoloc_lookup_dict[(i)][1].find('IT') < 0)):
        inout.log_message('Institution type word is not known: '+str(v),'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('postaddress_type')):  # Check postaddress type
    v = geoloc_dict['postaddress_type']
    if (len(v) > 2):
      inout.log_message('More than two postaddress type: '+str(v),'warn')
    for i in v:
      i = i.split('_')
      i = tuple(i)  # Make it a tuple
      if (not config.geoloc_lookup_dict.has_key((i))) or \
         (config.geoloc_lookup_dict.has_key((i)) and \
          (config.geoloc_lookup_dict[(i)][1].find('PA') < 0)):
        inout.log_message('Postaddress type word is not known: '+str(v),'warn')
        break  # Exit for loop

  # Check if 'qualifier' elements only contain known qualifier words  - - - - -
  #
  if (geoloc_dict.has_key('wayfare_qualifier')):  # Check wayfare qualifier
    v = geoloc_dict['wayfare_qualifier']
    for i in v:
      if (not config.geoloc_lookup_dict.has_key((i,))) or \
         (config.geoloc_lookup_dict.has_key((i,)) and \
          (config.geoloc_lookup_dict[(i,)][1].find('LQ') < 0)):
        inout.log_message('Wayfare qualifier word is not known: '+str(v), \
                          'warn')
        break  # Exit for loop

  if (geoloc_dict.has_key('locality_qualifier')):  # Check locality qualifier
    v = geoloc_dict['locality_qualifier']
    for i in v:
      if (not config.geoloc_lookup_dict.has_key((i,))) or \
         (config.geoloc_lookup_dict.has_key((i,)) and \
          (config.geoloc_lookup_dict[(i,)][1].find('LQ') < 0)):
        inout.log_message('Locality qualifier word is not known: '+str(v), \
                          'warn')
        break  # Exit for loop

  return geoloc_dict

# -----------------------------------------------------------------------------

def test():
  """Simple test routine with example inputs.
  """

  test_localities = ["2602 o'connor a.c.t", \
                     "27o2 o-connor austr. capit. ter.", \
                     "dickson 2602 a.c.t", \
                     "sydney nsw 2000", \
                     "haymarket 2000 new s wales", \
                     "sydney nsw 2001", \
                     "huskinson n-s-w 2407", \
                     "2600 Custance Street estern suburbs ACT 2607", \
                     "Zincke Cl. 15, ACT2907 byron bay north east", \
                     "4/8 Biddell Place east upper the entrance q.l.d. 2913", \
                     "16 Balonne Street 2617 Kalleen, Austr. Capit. Terr.", \
                     "42 victoria street, 2602 bawley point, s. a.", \
                     "42 victoria street, 2602 east boyd town, vic", \
                     "new south wales", \
                     "upper norah head west n-s-w", \
                     "bawley point heights n ter", \
                     "upper north sydney southern austr", \
                     "goldburn lower downs queen lnd" , \
                     " new victoria street 42, south australia", \
                    ]

  ilp_geocode_test_loc = ["100  cottenham ave kingsford 2032", \
                          "126  jamison rd penrith 2750", \
                          "29  violet st south bathurst 2795", \
                          "110  cawley st corrimal 2518", \
                          "39r  old mendooran rd dubbo 2830", \
                          "64  kanoona lane whittingham 2330", \
                          "36  campbell pde mannering park 2259", \
                          "504  collombatti rd collombatti 2440", \
                          "22  libra cl elermore vale 2287", \
                          "298a 298c mayfield rd pyree 2540", \
                          "127  kallaroo rd terrey hills 2084", \
                          "41  west lanitza rd lanitza 2460", \
                          "1115  beaties lane barellan 2665", \
                          "42  banks st monterey 2217", \
                          "7r  darling's rd dubbo 2830", \
                          "2b  kensington rd kensington 2033", \
                          "154  denison rd dulwich hill 2203", \
                          "21  brook rd glenbrook 2773", \
                          "10  attunga st seven hills 2147", \
                          "18  depot rd mortdale 2223", \
                          "139  bestic st brighton/le/sands 2216", \
                          "12  glenview rd mount kuring-gai 2080", \
                          "  deniliquin-barham rd bunnaloo 2731", \
                          "17 the king george gardents 2-8 hurstville 2220", \
                          "46  the grand pde brighton-le-sands 2216", \
                          "54  reading rd brighton le sands 2216", \
                          "9095  armidale-kempsey rd carrai 2440", \
                          "  bendemeer-watsons creek rd bendemeer 2355", \
                          "11  bilmark pl brighton-le-sands 2216", \
                          "286  rydal-hampton rd hampton 2790", \
                          "17  mundowi rd mount kuring'gai heights 2080", \
                          "  ghinni-ghi rd lower east geneva 2474", \
                          "  military land - holsworthy  holsworthy 2173", \
                          "1080  firth-heinz rd pillar valley 2462",\
                          "1080 firth-heinz rd un 1234 mnt vic 2462",\
                          "'the trees' newell h'way", \
                          "'bridley stud' via the bland", \
                          "'cloonawillin' duffy's lane", \
                          "'kunnawah' shephard's siding via", \
                          "'oak valley' rmb 142b", \
                          "rmb 122 rose valley creek", \
                          "conv home - 1359 pacific h", \
                          "diamond warers/diamond hd rd", \
                          "'winburn park' 11 mile park", \
                          "3/25- 35a park road", \
                          "three mile flats c/- po", \
                          "nsw academy of sport wakehurst parkway", \
                          "'shepherds lodge' lost river", \
                          "c/- 16-20 bardwell rd", \
                          "c/- 606/10 martin place", \
                          "c/- post office lower duck creek rd", \
                          "4 teatree road (po box400)", \
                          "c1/17 sunnyside ave", \
                          "lot 3 'marala' abbington park /barry way", \
                          "lot 1 pygmy possum pl woodlands east", \
                          "lot 1 wallagoot lane jellat jellat via", \
                          "lot 6 off old maitland rd", \
                          "lot 193 north dorrigo via", \
                          "2/lot 6 frank cooper st", \
                          "lot 1\\20 illawarra hwy", \
                          "lot 275 eight ave", \
                          "holbrook stud widden valley", \
                          "13/lot 3 equity pl", \
                          "lot 4, the northern road", \
                          "lot 4 the heights", \
                          "lot 3 lostock dam via", \
                          "portion 93 old great north rd", \
                          "3/30 cape three points rd", \
                         ]

  test_loc_all = ["corner Wellington Ave & Jackson Roads, Mulgrave", \
                  "123 smith street, smithfield NSW 2203", \
                  "7/2 smith st, north sydney, new south wales, 2201", \
                  "lot 3 miller-meyer arcade", \
                  "7a/5 corner smith miller nambucca heads 2347 n.s.w.", \
                  "lot 7/4 smith parade, wollongong", \
                  "'cheddar' unit 17a, smith street", \
                  "po box 19 berry north", \
                  "rmb 23, kangoroo valley", \
                  "'cheddar' via smalltown, large town, austr. cap. ter.", \
                  "123 smith road, smalltown via large city, 2690", \
                  "st. xavier nursing home, st. george road, 2560", \
                  "villa 17 mowell village", \
                  "van 42 lazy park caravan park", \
                  "c;van 17 lazy hill north park, 2210 n s w", \
                  "the boulevard, smithville", \
                  "corner round street and square avenue, victoria", \
                  "'Glencary' Sams Corner Rd", \
                  "Park Glen Cottage Kia-Ora Ln", \
                  "12/blk-B-2 Herbert Street", \
                 ]

  all_tests = test_localities[:] + ilp_geocode_test_loc[:]+test_loc_all[:]

  for l in all_tests:
    print
    print 'INPUT: "%s"'%(l)
    clean_str = clean_geoloc_component(l.lower())
    print 'CLEAN: "%s"'%(clean_str)
    [tmp_list,tag_list] = tag_geoloc_component(clean_str)
    print '  Words:',tmp_list
    print '  Tags: ',tag_list

    print get_geoloc_hmm(tmp_list, tag_list)
    print
    print
