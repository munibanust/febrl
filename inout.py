# =============================================================================
# inout.py - Routines related to file input and output.
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
# The Original Software is "inout.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module inout.py - Routines related to input and output.

   Public functions:
     load_lookup_tables  Load lookup-table files with word corrections and tags
                         into one dictionary
     load_corr_list      Load one correction list file into a sorted list
     process_line        Process a text line (one record) according to settings
                         from config.py
     compose_line        Compose a text line (one record) for the output file
                         according to settings from config.py
     check_field_spill   Routine to check if a known word is spilling from
                         one string into another
     log_message         Print and/or append and a message to a log file

   See doc strings of individual functions for detailed documentation.
"""

# -----------------------------------------------------------------------------

import string
import types
import os

import config
import name
import locality

# -----------------------------------------------------------------------------

def load_lookup_tables(file_names):
  """Load lookup-table files with word corrections and tags into one dictionary

  USAGE:
    [lookup_dict, max_key_length] = load_lookup_tables(file_names)

  ARGUMENTS:
    file_names  Either one file name (a string) or a list of file names.

  DESCRIPTION:
    This routine loads one or more lookup-table files and converts them into
    one Python dictionary.

    If one file name is given, the lookup-table file is loaded and converted
    into a Python dictionary. If a list of file names is given, they are all
    loaded and merged into one single dictionary.

    A lookup-table file can contain one or more blocks of entires. A block
    starts with a line that contains a tag assignment with a tag in brackets:

      tag=<tag>

    This tag assignment must be written at the beginning of a line. It is
    possible to have a comment (starting with #) after the assignmnt.

    The following lines then contain the entries that will be tagged with the
    current tag. Each entry is of the form:

                key : values

    where values is either a string (space sbetwee nwords are allowed) or a
    comma separated list of strings. All values will become keys in the
    dictionary, and the key becomes the corresponding value. A value can
    consist of more than one word (spaces between words are possible).

    The list of values can be longer than one line. In such a case, simply
    write a text line with values without a colon. Example:

           northern : no, northern, nrth, nrthn, nrthrn, nth,
                      nthn, nthrn

    If a value occurs more than once in the lookup-table files and is expanded
    or corrected into the same string but with different tags, all tags are
    kept and stored in a list for this value.
    If a value occurs more than once in the lookup-table files and is expanded
    or corrected into different strings, an error message is printed and the
    program stops. The user then manually has to correct this error.

    Here is an example of a lookup-table file:

      tag=<WT>  # Tag for wayfare type words
               road : rd
             street : st, str

      tag=<CR>  # Tag for country words
        new zealand : new sealand, neu seeland, nz
  """

  # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (type(file_names) == types.StringType):
    file_names = [file_names]  # Make a list out of a single file name

  if (type(file_names) != types.ListType):
    msg = ['Input argument is of wrong type: '+file_names, \
           'Must be either of type "string" or "list"']
    log_message(msg,'err')
    raise Exception()

  dict = {}  # New empty dictionary
  max_len = 1  # Maximal length of the longest key sequence

  # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for fn in file_names:

    # Open file and read all lines into a list
    #
    try:
      f = open(fn,'r')
    except:
      log_message('Cannot open file: '+fn,'err')
      raise IOError()

    file_data = f.readlines()
    f.close()

    tag = ''  # Start with no tag

    # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for line in file_data:
      l = line.strip()  # Remove line separators
      if (len(l) > 0) and (l[0] != '#'):  # Not empty line and no comment line

        if (l[:5] == 'tag=<'):  # It's a line with a new tag
          tag = l[5:7]
          if (tag not in config.name_hmm_obser) and \
             (tag not in config.geoloc_hmm_obser):  # Make sure tag is valid
            log_message('Illegal tag: '+tag,'err')
            raise Exception()

        else:  # A line with an entry

          if (tag == ''):  # make sure a tag is set
            log_message('No tag set in file: '+fn,'err')
            raise Exception()

          ll = l.split(':')  # Separate key from values

          if (len(ll) == 2):  # Line contains a key - - - - - - - - - - - - - -
            k = ll[0].strip().lower()  # Get key, make lower and strip spaces

            k_list = k.split(' ')  # Make a list of key words
            if (len(k_list) > max_len):  # Update maximal key sequence length
              max_len = len(k_list)

            val = string.join(k_list,'_')
            key = tuple(k_list)
            this_tag = tag

            if (k != ''):  # If key is non-empty insert it into dictionary
              if (dict.has_key(key)):
                test_val = dict[key][0]  # Value without tag
                test_tag = dict[key][1]

                if (val == test_val):  # Same values
                  if (test_tag.find(this_tag) < 0):  # This tag is not in tags
                    this_tag = test_tag+'/'+this_tag
                else:
                  msg = ['Key already in dictionary with different value', \
                         'Key: "'+str(key)+'", old value: "'+ \
                         str(dict[key][0])+'", new value: "'+str(val)+'"']
                  log_message(msg,'err')
                  raise Exception()

              this_val = (val, this_tag)
              dict.update({key:this_val})  # Insert key itself into dicionary

            v = ll[1].lower() # Get values in a string

          elif (len(ll) == 1):  # Line contains only values - - - - - - - - - -
            v = ll[0].lower() # Get values in a string

          else:
            log_message('Illegal file format in file: '+fn+', line: '+l,'err')
            raise Exception()

          vv = v.split(',')  # Split values into a list

          for v in vv:  # Loop over all values  - - - - - - - - - - - - - - - -
            vs = v.strip()
            if (vs != ''):  # Only append non-empty values
              k_list = vs.split(' ')  # Make a list of key words
              if (len(k_list) > max_len):  # Update maximal key sequence length
                max_len = len(k_list)
              key = tuple(k_list)
              this_tag = tag

              if (dict.has_key(key)):
                test_val = dict[key][0]  # Value without tag
                test_tag = dict[key][1]

                if (val == test_val):  # Same values
                  if (test_tag.find(this_tag) < 0):  # This tag is not in tags
                    this_tag = test_tag+'/'+this_tag
                else:
                  msg = ['Key already in dictionary with different value', \
                         'Key: "'+str(key)+'", old value: "'+ \
                         str(dict[key][0])+'", new value: "'+str(val)+'"']
                  log_message(msg,'err')
                  raise Exception()

              this_val = (val, this_tag)
              dict.update({key:this_val})  # Insert key itself into dicionary

  return [dict, max_len]

# -----------------------------------------------------------------------------

def load_corr_list(file_name):
  """Load one correction list file into a sorted list

  USAGE:
    corr_list = load_corr_list((file_name)

  ARGUMENTS:
    file_name  The name of a correction list file ('.lst')

  DESCRIPTION:
    This routine loads one correction list file into a list sorted according to
    the length (longest first) of the orignal strings to be replaced.

    Entries in the file are of the form:

          replacement := values

    where 'values' can be one word or character or a comma separated list of
    words or characters. Each value will be replaced by the 'replacement'
    string on the left side.
    Both the replacement string and each of the value strings must be enclosed
    in either single or double quotes.
  """

  # Check input argument type and open file - - - - - - - - - - - - - - - - - -
  #
  if (type(file_name) != types.StringType):
    log_message('File name is not of type "string": '+file_name,'err')
    raise Exception()

  try:
    f = open(file_name,'r')
  except:
    log_message('Cannot open file: '+file_name,'err')
    raise IOError()

  file_data = f.readlines()
  f.close()

  org_list  = []  # List of original strings (the ones to be replaced)
  repl_list = []  # List of replacement strings
  len_list  = []  # List of original string lengths
  repl = ''  # Set inital replacement to nothing

  # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for line in file_data:
    l = line.strip()  # Remove line separators at the end

    if (len(l) > 0) and (l[0] != '#'):  # Not an empty line and no comment line
      ll = l.split(':=')  # Separate replacement from values

      if (len(ll) == 2):  # Line contains a replacement - - - - - - - - - - - -

        repl = ll[0].strip().lower()  # Make replacement lower and strip spaces

        if (not ((repl[0] == '"') and (repl[-1] == '"') or \
                 (repl[0] == "'") and (repl[-1] == "'"))):
          log_message(['Replacement string is not properly quoted: '+repl, \
                       'In file: '+file_name], 'err')
          raise Exception()

        repl = repl[1:-1]  # Remove quotes from replacement string

        v = ll[1].lower() # Get values in a string and make lowercase

      elif (len(ll) == 1):  # Line contains only values - - - - - - - - - - - -
        v = ll[0].lower() # Get values in a string and make lowercase

      else:  # More than one ':=' separator in the line - - - - - - - - - - - -
        log_message('Too many ":=" separators in this line: '+l, 'err')
        raise Exception()

      # Now process the values and append them to the list  - - - - - - - - - -

      vv = v.split(',')  # Split values into a list

      for v in vv:  # Loop over all values  - - - - - - - - - - - - - - - - -
        org = v.strip()  # Get the original string

        if (org != ''):  # Only process non-empty values
          if (not ((org[0] == '"') and (org[-1] == '"') or \
                   (org[0] == "'") and (org[-1] == "'"))):
            log_message(['Original string is not properly quoted: '+org, \
                         'In file: '+file_name], 'err')
            raise Exception()

          org = org[1:-1]  # Remove quotes from original string

          if (org != ''):  # Only append non-empty values
            org_list.append(org)
            repl_list.append(repl)
            len_list.append(len(org))

  tmp_list = map(None,len_list,org_list,repl_list)
  tmp_list.sort()
  tmp_list.reverse()

  corr_list = []
  for (i,org,repl) in tmp_list:
    corr_list.append((org,repl))

  msg = ['Loaded "'+file_name+'" into correction list"']  # Prepare log message
  for i in corr_list:
    msg.append('  '+str(i))
  log_message(msg,'v2')

  return corr_list

# -----------------------------------------------------------------------------

def load_dict(file_names):
  """Load one or more comma separated text files into a lookup-dictionary.

  USAGE:
    lookup_dict = load_dict(file_names)

  ARGUMENTS:
    file_names  Either one file name (a string) or a list of file names.

  DESCRIPTION:
    If one file name is given, the lookup-table file is loaded and converted
    into a Python dictionary. If a list of file names is given, they are all
    loaded and merged into one single dictionary.

    For more than one input files, warnings are printed out if keys occur more
    than once, and if they have different values.

    The files have to be text files in a two-column format (CSV files). The
    first column becomes the keys of the dictionary, the second the values.

  EXAMPLES:
    title_dict  = loaddict('title.csv')

    suburb_dict = load_dict(['suburb_act.csv', 'suburb_nsw.csv'])
  """

  # Check input argument  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (type(file_names) == types.StringType):
    file_names = [file_names]  # Make a list out of a single file name

  if (type(file_names) != types.ListType):
    msg = ['Input argument is of wrong type: '+file_names, \
           'Must be either of type string or list']
    log_message(msg,'err')
    raise TypeError(load_dict)

  dict = {}  # New empty dictionary

  # Loop over file names  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for fn in file_names:

    # Open file and read all lines into a list
    #
    try:
      f = open(fn,'r')
    except:
      log_message('Cannot open file: '+fn,'err')
      raise IOError(load_dict)

    file_data = f.readlines()
    f.close()

    # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for line in file_data:
      l = lline.strip()
      if (len(l) > 0) and (l[0] != '#'):  # Not empty line and no comment line
        ll = l.split(',')  # Get fields from a line
        if (len(ll) != 2):  # Not two columns
          msg = 'Illegal file format (not 2 columns) in file: '+fn+', line: '+l
          log_message(msg,'err')
          raise Exception(load_dict)

        k = ll[0].strip().lower()  # Make sure it's lower case
        v = ll[1].strip().lower()

        if (dict.has_key(k)):
          if (dict[k] != v):
            msg = ['Key '+v+' already in dictionary with different value.', \
                   '(old value will be over written)', \
                   'Old value: '+dict[v]+', new value: '+k]
            log_message(msg, 'warn')
            dict.update({k:v})

        else:  # New key
          dict.update({k:v})

  return dict

# -----------------------------------------------------------------------------

def process_line(line):
  """Process a text line (one record) according to settings from config.py

  USAGE:
    [name_comp, geocode_comp, locality_comp, date1_comp, date2_comp] = \
                          process_line(line)
 
  ARGUMENTS
    line         The input text line (a string) containing a record

  DESCRIPTION:
    This routine gets a raw text line as input, and splits it into five
    components (name, geocode, locality, date1 and date2). Splitting is
    done according to the file type, and field quote as set in config.py

    For comma and tabulator separated files the csv.py module (Object-Craft)
    is used.

    The routine also checks if words are 'spilling' from one field into the
    next (e.g. "miller r","d canberra") and tries to fix such errors by
    concatenating adjacent words and then looking them up in various lookup-
    tables. If a concatenated word is found in a lookup-table, the corrected
    word is inserted (instead of the two separate parts).
    This functionality can be activated in 'config.py'.

    Returns a list containing five components, each of them is normally a
    string. If 'givenname' and 'surname' are in separate fields (as specified
    in config.py), the name component is returned as a list with the two
    entries [givenname_component, surname_component].
  """

  name_comp_list       = []
  givenname_comp_list  = []
  surname_comp_list    = []
  geocode_comp_list    = []
  locality_comp_list   = []
  date1_comp_list      = []
  date2_comp_list      = []

  # Split the line into the basic fields  - - - - - - - - - - - - - - - - - - -
  #
  if (config.in_file_type in ['CSV','CSVQ','TAB','TABQ']):
                                                 # Comma or tabulator separated
    try:
      line_list = config.line_parser.parse(line)
    except:
      log_message('CSV line parsing failed with inout: '+line,'err')

    if (len(line_list) < config.input_len):
      log_message('Input line does not contain enough fields,' +\
                  'fill up with empty fields','warn')
      while (len(line_list) < config.input_len):
        line_list.append('')

    config.curr_line_list = line_list  # Save current line list

    # Extract fields into different component lists - - - - - - - - - - - - - -
    #
    if (config.input_component['name'] != []):  # Extract name fields
      for i in config.input_component['name']:
        name_comp_list.append(line_list[i])

    else:  # Extract givenname and surname into separate components - - - - - -
      if (config.input_component['givenname'] != []):  # Extract g-name fields
        for i in config.input_component['givenname']:
          givenname_comp_list.append(line_list[i])

      if (config.input_component['surname'] != []):  # Extract surname fields
        for i in config.input_component['surname']:
          surname_comp_list.append(line_list[i])

    if (config.input_component['geocode'] != []):  # Extract geocode fields
      for i in config.input_component['geocode']:
        geocode_comp_list.append(line_list[i])

    if (config.input_component['locality'] != []):  # Extract locality fields
      for i in config.input_component['locality']:
        locality_comp_list.append(line_list[i])

    if (config.input_component['date1'] != []):  # Extract date1 fields
      for i in config.input_component['date1']:
        date1_comp_list.append(line_list[i])

    if (config.input_component['date2'] != []):  # Extract date2 fields
      for i in config.input_component['date2']:
        date2_comp_list.append(line_list[i])

  elif (config.in_file_type == 'COL'):  # Column based input file - - - - - - -

    if (len(line) < config.input_len):
      log_message('Input line is not long enough, fill up with spaces','warn')
      line += ' '*(config.input_len-len(line))

    if (config.input_component['name'] != []):  # Extract name fields
      for (col_start,length) in config.input_component['name']:
        name_comp_list.append(line[col_start,col_start+length])

    else:  # Extract givenname and surname into separate components - - - - - -
      if (config.input_component['givenname'] != []):  # Extract g-name fields
        for (col_start,length) in config.input_component['givenname']:
          givenname_comp_list.append(line[col_start,col_start+length])

      if (config.input_component['surname'] != []):  # Extract surname fields
        for (col_start,length) in config.input_component['surname']:
          surname_comp_list.append(line[col_start,col_start+length])

    if (config.input_component['geocode'] != []):  # Extract geocode fields
      for (col_start,length) in config.input_component['geocode']:
        geocode_comp_list.append(line[col_start,col_start+length])

    if (config.input_component['locality'] != []):  # Extract locality fields
      for (col_start,length) in config.input_component['locality']:
        locality_comp_list.append(line[col_start,col_start+length])

    if (config.input_component['date1'] != []):  # Extract date1 fields
      for (col_start,length) in config.input_component['date1']:
        date1_comp_list.append(line[col_start,col_start+length])

    if (config.input_component['date2'] != []):  # Extract date2 fields
      for (col_start,length) in config.input_component['date2']:
        date2_comp_list.append(line[col_start,col_start+length])

  # elif (config.in_file_type == 'SQL'):  # - - - - - - - - - - - - - - - - - -

    ################################
    # Add later: SQL database access
    ################################

  msg = ['  Component basic field lists:', \
         '    Name:       '+str(name_comp_list), \
         '    Given name: '+str(givenname_comp_list), \
         '    Surname:    '+str(surname_comp_list), \
         '    Geocode:    '+str(geocode_comp_list), \
         '    Locality:   '+str(locality_comp_list), \
         '    Date1:      '+str(date1_comp_list), \
         '    Date2:      '+str(date2_comp_list)]
  log_message(msg,'v2')

  name_comp      = ''
  givenname_comp = ''
  surname_comp   = ''
  geocode_comp   = ''
  locality_comp  = ''
  date1_comp     = ''
  date2_comp     = ''

  # Now clean and then concatenate component lists into strings - - - - - - - -
  #
  if (name_comp_list != []):  # Name component
    name_comp = name_comp_list[0]  # Start with first field in list

    for f in name_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['name'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['name'] == 1):
          sep = check_field_spill(name_comp, f)

        name_comp = name_comp+sep+f  # Append separator and field

  if (givenname_comp_list != []):  # Givenname component  - - - - - - - - - - -
    givenname_comp = givenname_comp_list[0]  # Start with first field in list

    for f in givenname_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['givenname'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['givenname'] == 1):
          sep = check_field_spill(givenname_comp, f)

        givenname_comp = givenname_comp+sep+f  # Append separator and field

  if (surname_comp_list != []):  # Surname component  - - - - - - - - - - - - -
    surname_comp = surname_comp_list[0]  # Start with first field in list

    for f in surname_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['surname'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['surname'] == 1):
          sep = check_field_spill(surname_comp, f)

        surname_comp = surname_comp+sep+f  # Append separator and field

  if (geocode_comp_list != []):  # Geocode component  - - - - - - - - - - - - -
    geocode_comp = geocode_comp_list[0]  # Start with first field in list

    for f in geocode_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['geocode'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['geocode'] == 1):
          sep = check_field_spill(geocode_comp, f)

        geocode_comp = geocode_comp+sep+f  # Append separator and field

  if (locality_comp_list != []):  # Locality component  - - - - - - - - - - - -
    locality_comp = locality_comp_list[0]  # Start with first field in list

    for f in locality_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['locality'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['locality'] == 1):
          sep = check_field_spill(locality_comp, f)

        locality_comp = locality_comp+sep+f  # Append separator and field

  if (date1_comp_list != []):  # Date1 component  - - - - - - - - - - - - - - -
    date1_comp = date1_comp_list[0]  # Start with first field in list

    for f in date1_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['date1'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['date1'] == 1):
          if (date1_comp[-1] != ' ') and (f[0] != ' '):
            tmp_list0 = date1_comp.split()
            tmp_list1 = f.split()
            check_word = tmp_list0[-1]+tmp_list1[0]

          if (check_word in ['jan','feb','mar','apr','may','jun','jul','aug', \
              'sep','oct','nov','dec','january','february','march','april', \
              'may','june','july','august','september','october','november', \
              'december']):

              sep = ''  # Set separator to no space
              msg = '  Correct date1 word spilling: "'+date1_comp+'","'+f+'"'
              log_message(msg,'v1')

        date1_comp = date1_comp+sep+f  # Append separator and field

  if (date2_comp_list != []):  # Date2 component  - - - - - - - - - - - - - - -
    date2_comp = date2_comp_list[0]  # Start with first field in list

    for f in date2_comp_list[1:]:  # Loop over following fields (if any)
      if (f != ''):
        if (config.input_space_sep['date2'] == 1):
          sep = ' '  # Set separator to space between fields
        else:
          sep = ''  # No space between fields

        # Check field spilling only if space separator is set to ' ' 
        #
        if (sep == ' ') and (config.input_check_spilling['date2'] == 1):
          if (date2_comp[-1] != ' ') and (f[0] != ' '):
            tmp_list0 = date1_comp.split()
            tmp_list1 = f.split()
            check_word = tmp_list0[-1]+tmp_list1[0]

          if (check_word in ['jan','feb','mar','apr','may','jun','jul','aug', \
              'sep','oct','nov','dec','january','february','march','april', \
              'may','june','july','august','september','october','november', \
              'december']):

              sep = ''  # Set separator to no space
              msg = '  Correct date1 word spilling: "'+date1_comp+'","'+f+'"'
              log_message(msg,'v1')

        date2_comp = date2_comp+sep+f  # Append separator and field

  # Check if name component is given or givenname and surname separately - - -
  #
  if (config.input_component['givenname'] != []) or \
     (config.input_component['surname'] != []):
     name_comp = [givenname_comp, surname_comp]

  msg = ['  Components:', \
         '    Name:       "'+str(name_comp)+'"', \
         '    Geocode:    "'+geocode_comp+'"', \
         '    Locality:   "'+locality_comp+'"', \
         '    Date1:      "'+date1_comp+'"', \
         '    Date2:      "'+date2_comp+'"']
  log_message(msg,'v1')

  return [name_comp, geocode_comp, locality_comp, date1_comp, date2_comp]

# -----------------------------------------------------------------------------

def compose_line(output_dict,header=None):
  """Compose a text line (one record) for the output file.

  USAGE:
    line_text = compose_line(output_dict)
 
  ARGUMENTS
    output_dict  The dictionary containing the standardised output values
    header       A flag, set to '1' if the 'output_dict' contains a header
                 line and not a data line. Default is None (not a header line)

  DESCRIPTION:
    This routine gets as input a dictionary with the output field values which
    are then processed according to the output file type and output file
    settings from 'config.py' into one single string that is returned.
  """

  line = ''

  # Loop over all activated output fields - - - - - - - - - - - - - - - - - - -
  #

  for (field_pos,field_name) in config.output_field_list:

    if (header == 1):  # This is a header line, field values are field names
      field = output_dict[field_name]

    elif (field_name == 'record_id'):  # Add current record (line) number
      field = str(config.curr_line_no)

    elif (field_name[:14] == 'original_input'):  # Output original input
 
      if (len(field_name) == 14):
        field = config.curr_line  # The complete input line (without line sep.)
      else:  # There is an input field or input column range given
        v = field_name[14:].strip()  # Get the value
        v = v[1:-1]  # Remove brackets '[' and ']' from beginning and end
        if (v[0] == '(') and (v[-1] == ')'):  # It's a tuple
          v = v[1:-1]  # Remove tuple brackets
        v = v.split(',')  # Make a list
        for i in range(len(v)):
          v[i] = int(v[i])  # Make integers
        if (len(v) == 1):  # A input field number
          field = config.curr_line_list[v[0]]
        else:  # Length must be two, for column start and end
          field = config.curr_line(v[0],v[1]+1)

    elif (output_dict.has_key(field_name)):
      field = output_dict[field_name]  # Get current output field

      if (field == []):
        field = ''  # Make all empty fields to empty strings
      elif (type(field) == types.ListType):
        field = string.join(field,' ')  # Convert into string (space separator)

      if (type(field) != types.StringType):  # Illegal field type
        log_message('Output field is of illegal type: '+str(field),'err')
        raise Exception()

    else:  # No current output field
      field = ''  # Create an empty output field

    # Now append current field to output line - - - - - - - - - - - - - - - - -
    #
    if (config.out_file_type in ['CSV','TAB']):  # Comma or tabulator separated
      if (',' in field):  # The field contains a comma, quote it
        field = config.output_quote_character + field + \
                config.output_quote_character
      line = line + field + config.out_field_sep

    if (config.out_file_type in ['CSVQ','TABQ']):  # Quoted comma or tabulator
      line = line + config.output_quote_character + field + \
             config.output_quote_character + config.out_field_sep

    elif (config.out_file_type == 'COL'):
      [(start_col, field_len)] = field_pos
      if (len(field) > field_len):  # Field content is too long
        log_message('Field content too long: '+str(field)+' (max '+ \
                    str(field_len)+' characters)', 'warn')
      field = field.ljust(field_len)  # Expand field with spaces

      line = line + field[:field_len]  # Append field to line, cut trail off if
                                       # too long

    # elif (config.out_file_type == 'SQL'):

      ################################
      # Add later: SQL database access
      ################################

  if (config.out_file_type in ['CSV','CSVQ','TAB','TABQ']):
    line = line[:-len(config.out_field_sep)]  # Remove last field separator

  log_message('  Output line: |'+line+'|','v2')

  return line

# -----------------------------------------------------------------------------

def check_field_spill(str1, str2):
  """Routine to check if a known word is spilling from one string into another

  Returns a separator, either '' (if a word spilling has been found) or ' '
  otherwise.
  """

  org_str1 = str1  # Keep copies of the original input string
  org_str2 = str2

  alphanum = string.letters+string.digits  # String with all letters and digits
  sep = ' '  # Default return value for separator

  if (str1 == '') or (str1 == ''):  # One input string is empty
    pass  # Make sure the strings are not empty

  # Check if the characters at the 'boundary' are either letters or digits
  #
  elif (str1[-1] in alphanum) and (str2[0] in alphanum):
    test_str1 = ''
    while (str1 != '') and (str1[-1] in alphanum): # Get last char. from str1
      test_str1 =  str1[-1] + test_str1
      str1 = str1[:-1]
    test_str2 = ''
    while (str2 != '') and (str2[0] in alphanum): # Get first char. from str2
      test_str2 = test_str2 + str2[0]
      str2 = str2[1:]

    check_word = test_str1 + test_str2

    if (config.name_lookup_dict.has_key((check_word,)) or \
        config.geoloc_lookup_dict.has_key((check_word,))):
      sep = ''  # Set separator to no space
      msg = '  Correct name word spilling: "'+org_str1.strip()+ \
            '","'+org_str2.strip()+'" -> "'+check_word+'"'
      log_message(msg,'v2')
      config.num_word_spills += 1  # Increment number of corrected word spills

  return sep

# -----------------------------------------------------------------------------

def log_message(message, msg_type):
  """Print and/or append and a message to a log file.

  USAGE:
    log_message(message, msg_type)

  ARGUMENTS:
    message   Either a string containing the message (a one liner) or a list of
              strings (several lines).
    msg_type  The type of the message. Possible message types:
                'warn' Warning information
                'err'  Error information
                'v1'   Information message of verbose level 1
                'v2'   Information message of verbose level 2

  DESCRIPTION:
    Depending on the values of 'config.verbose' and 'config.logging' the
    message will be printed and/or logged to file.

    Warning and error messages will always be printed and/or logged. They are
    printed and/or logged with a '***** Warning' or '***** Error **********'
    before the message, and a '***** ' at the beginning of each line of the
    message.

    The difference between warning and error messages is that error messages
    usually are printed before a program is stopped due to a non-recoverable
    exception, while after warning messages the program continues.

    Messages of type 'v1' and 'v2' will be printed if they are equal or smaller
    than the verbose setting in 'config.py'. For example, if the verbose
    setting is level 2, both messages of type 'v1' and 'v2' will be printed,
    but if verbose level is 1, only messages of type 'v1' are printed.

    If logging is activated (config.logging > 0) the message will be appended
    (line by line) to the log file.

    If the 'nowarn' flag is set in the project module (or given as optional
    command line argument) no warning messages are printed, but they are still
    written to the log file if logging is activated.

    The name of the log file is defined in 'config.py' and can be overwritten
    if the option '-l [file_name]' is given in 'pyStandard.py' or
    'pyLinkage.py'.

    The routine appends a line separator from the 'os' module (os.linesep) at
    the end of each line when a message is logged.
  """

  if (config.verbose == 0) and (config.logging == 0) and \
     (msg_type not in ['err','warn']):
    return None # A normal information message and no verbose output or logging

  if (msg_type == 'v1'):
    msg_level = 1
    msg_type = 'info'
  elif (msg_type == 'v2'):
    msg_level = 2
    msg_type = 'info'
  elif (msg_type not in ['err','warn']): 
    print
    print '***** Error **********'
    print '***** Illegal message type:', msg_type
    raise Exception()

  # Create message header for warning and error messages  - - - - - - - - - - -
  #
  if (msg_type == 'err'):
    if (config.curr_line_no >= 0):
      m = ['', '***** Error! Line number: '+str(config.curr_line_no)]
    else:
      m = ['', '***** Error!']
    if (config.curr_line != ''):
      m.append('***** Line: "'+config.curr_line+'"')
  elif (msg_type == 'warn'):
    config.num_warning += 1  # Increment number of warnings
    if (config.curr_line_no >= 0):
      m = ['', '***** Warning! Line number: '+str(config.curr_line_no)]
    else:
      m = ['', '***** Warning!']
    if (config.curr_line != ''):
      m.append('***** Line: "'+config.curr_line+'"')
  else:
    m = []

  # Make sure message is of type 'list' - - - - - - - - - - - - - - - - - - - -
  #
  if (type(message) == types.StringType):
    message = [message]
  elif (type(message) != types.ListType):  # Illegal type of the message
    print
    print '***** Error **********'
    print '***** Illegal type of message data:', type(message)
    print '***** Must be either of type string or list'
    raise TypeError()

  for msg in message:
    if (msg_type == 'err') or (msg_type == 'warn'):
      m.append('***** '+msg)
    else:
      m.append(msg)

  if (msg_type in ['warn','err']) or (config.verbose >= msg_level):

    if (config.logging > 0):  # Write message to log file
      try:
        f_log = open(config.log_file, 'a')  # Open file for appending
      except:
        print
        print '***** Error **********'
        print '***** Cannot open log file:', log_file
        raise IOError()

      for line in m:
        f_log.write(line+os.linesep)  # Append message to log file

      f_log.close()

    if (msg_type != 'warn') or ((msg_type == 'warn') and (config.nowarn == 0)):
      for line in m:  # Print message
        print line

  return None

# -----------------------------------------------------------------------------
