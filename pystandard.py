# =============================================================================
# pyStandard.py - Main module for data parsing and standardisation
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
# The Original Software is "pyStandard.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University), Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health) and Drs Markus Hegland, Stephen Roberts and Ole Nielsen
# (Mathematical Sciences Insitute, Australian National University). Copyright
# (C) 2002 the Australian National University and others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module PYstandard.py - Main module for data parsing and standardisation

   USAGE:
     pyStandard.py [project_module] [first_rec] [num_rec] [options]

   ARGUMENTS:
     project_module  A Python module containing all project settings.
     first_rec       Number of the first record to be processed
     num_rec         Number of records to be processed
     options         Optional arguments, see below.

   OPTIONS (after arguments);
     -h                     Write header (output field names) to output file
                            (first line)
     -hmm-name [file_name]  Load and use a Hidden Markov Model (HMM) for name
                            standardisation (overwrite default from config.py)
     -hmm-loc [file_name]   Load and use a Hidden Markov Model (HMM) for
                            locality standardisation (overwrite default from
                            config.py)
     -l [file_name]         Log warnings and errors into the given file
     -v1                    Verbose output level 1 (low output)
     -v2                    Verbose output level 2 (high output)
     -nowarn                Don't print warning messages

   DESCRIPTION:
     This is the main module to standardise a text file containing records
     with personal information, like names, addresses and dates.

     All necessary project options (like input and output file names and types,
     field and component definitions, etc) are defined within a project module.

     Records are skipped until 'first_rec' is reached, from then on 'num_rec'
     records are processed. It is assumed that one line in the input file
     corresponds to one record.

     By using the '-hmm-name' ond/or '-hmm-loc' options it is possible to
     overwrite the settings in the project module and set the standardisation
     methods to HMM processing.

     For verbose output, add the option '-v1' or '-v2' with the second option
     producing more detailed output than the first one.

     To suppress the printing of warning messages used the '-nowarn' option.
     If logging is activated, warning messages are always logged.

     If you want to save warning and error messages into a file, use the '-l'
     option followed by the file name for the log file. If the log file already
     exists, new information will be appended, otherwise the file will be
     created. Logging is done according to the verbose level. If no verbose
     level is given, only warning and error messages will be saved to the log
     file.

     Currently, only sequential processing is possible. Parallel functionality
     will be included in a future release of this software. 

   EXAMPLE:
     python pyStandard.py my-project.py 0 100000 -v1 -l
"""

# -----------------------------------------------------------------------------

import sys
import os
import xreadlines
import types
import string
import time

import config
import inout
import simplehmm

import name
import locality
import date

# -----------------------------------------------------------------------------

def standard():
  """Main routine, open file, read lines, standardise them and write into file.

  USAGE:
    standard()

  ARGUMENTS:
    None

  DESCRIPTION:
    Main routine, see description of module above.
  """

  # Process command line arguments and check for correctness  - - - - - - - - -
  #
  if (len(config.options) < 2):
    print '***** Error: %s needs at least three arguments:'% (sys.argv[0])
    print '*****        - Name of the project module'
    print '*****        - Number of the first record to be processed'
    print '*****        - Number of records to be processed'
    print '*****          plus options'
    raise Exception()

  first_rec = int(config.options[0])
  num_rec   = int(config.options[1])
  in_file_name = config.in_file_name
  out_file_name = config.out_file_name

  # Check for optional arguments and process if any - - - - - - - - - - - - - -
  #
  config.verbose = 0  # Default: No verbose output
  config.logging = 0  # Default: No logging into a file
  write_header   = 0  # Write header (output field names) to output file
                      # (default: Don't)
  config.nowarn  = 0  # Deactivate no warning flag (print/log warning messages)

  if (len(config.options) > 2):
    options = config.options[2:]
    while (options != []):  # Do a loop processing all options

      if (options[0] == '-nowarn'):
        config.nowarn = 1  # Activate no warning flag
        options = options[1:]  # Remove processed '-nowarn' option

      elif (options[0] == '-v1'):
        config.verbose = 1  # Set to verbose output level 1
        options = options[1:]  # Remove processed '-v1' option

      elif (options[0] == '-v2'):
        config.verbose = 2  # Set to verbose output level 2
        options = options[1:]  # Remove processed '-v2' option

      elif (options[0] == '-l'):
        config.logging = 1
        if (len(options) > 1):
          if (options[1][0] != '-'):  # Not another option, must be a file name
            config.log_file = options[1]  # Get name of log file
            options = options[1:]  # Remove file_name
        options = options[1:]  # Remove processed -'l' option only

        try:
          f_log = open(config.log_file,'a')  # Test if file is appendable
        except:
          print '***** Error ********************',
          print '***** Cannot write to log file:', config.log_file
          raise IOError()

        # Write (append) header to log file
        #
        f_log.write(os.linesep)
        f_log.write('##################################################')
        f_log.write("############"+os.linesep)
        f_log.write("#"+os.linesep)
        f_log.write("# 'pyStandard.py - Version 0.1' process started at: ")
        f_log.write(time.ctime(time.time())+os.linesep)
        f_log.write("#"+os.linesep)
        f_log.write("# Input file name:  "+in_file_name+os.linesep)
        f_log.write("# Output file name: "+out_file_name+os.linesep)
        f_log.write(os.linesep)
        f_log.close()

      elif (options[0] == '-h'):
        write_header = 1
        options = options[1:]  # Remove processed -'h' option

      elif (options[0] == '-hmm-name'):
        hmm_name_file = options[1]  # Get file name of the name HMM to use
        try:
          f_in = open(hmm_name_file,'r')  # Test if file is available
        except:
          print '***** Error ********************',
          print '***** Cannot open HMM file in "-hmm-name" option:',
          print hmm_name_file
          raise IOError()

        f_in.close()
        options = options[2:]  # Remove processed option and file name
        config.name_standard_method = 'hmm'
        config.name_hmm_file_name = hmm_name_file
        config.name_hmm = simplehmm.hmm([],[])  # Create new empty HMM object
        config.name_hmm.load_hmm(config.name_hmm_file_name)

      elif (options[0] == '-hmm-loc'):
        hmm_loc_file = options[1]  # Get file name of the locality HMM to use
        try:
          f_in = open(hmm_loc_file,'r')  # Test if file is available
        except:
          print '***** Error ********************',
          print '***** Cannot open HMM file in "-hmm-loc" option:',
          print hmm_loc_file
          raise IOError()
        f_in.close()
        options = options[2:]  # Remove processed option and file name
        config.geoloc_standard_method == 'hmm'
        config.geoloc_hmm_file_name = hmm_loc_file
        config.geoloc_hmm = simplehmm.hmm([],[])  # Create new HMM object
        config.geoloc_hmm.load_hmm(config.geoloc_hmm_file_name)

      else:
        print '***** Error: Illegal option:', options[0]
        raise Exception()

  # Open input file and check number of available records - - - - - - - - - - -
  #
  try:
    f_in = open(in_file_name,'r')
  except:
    inout.log_message('Cannot open input file: '+in_file_name,'err')
    raise IOError()

  line_count = 0
  for line in f_in.xreadlines():
    line_count += 1
  f_in.close()

  if ((first_rec+num_rec) > line_count):  # Illegal value for last record
    print '***** Error: Illegal values for number of records to process:',
    print num__rec, ', with start record:', start_rec
    print '*****        File only contains',line_count, 'lines/records'
    raise Exception()

  # Open files  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  try:
    f_in = open(in_file_name,'r')
  except:
    inout.log_message('Cannot open input file: '+in_file_name,'err')
    raise IOError()

  try:
    f_out = open(out_file_name,'w')
  except:
    inout.log_message('Cannot open output file: '+out_file_name,'err')
    raise IOError()

  # Write header (name of output fields) into output file - - - - - - - - - - -
  #
  if (write_header == 1):
    header_dict = {}
    for n in config.output_field_names:
      header_dict.update({n:n})  # Dictionary where values are field names

    header_line = inout.compose_line(header_dict,header=1)
    f_out.write(header_line+os.linesep)

  # Skip over records - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (first_rec > 0):
    for i in range(first_rec):
      f_in.readline()

  # Read lines, process them and write into output files  - - - - - - - - - - -
  #
  line_read = 0  # Number of read lines

  while (line_read < num_rec):  # Loop until 'num_rec' records processed
    line = f_in.readline()

    # Print process indicator message
    #
    if (config.proc_ind >= 0) and (line_read > 0):  # Only print if activated
      if (line_read % config.proc_ind == 0):
        print 'Processed line', line_read, 'of', num_rec

    line = line.strip()  # Remove line separators
    config.curr_line = line  # Make a copy of the unprocessed current line

    line = line.lower()  # Make all characters lower case

    inout.log_message(['Record '+str(line_read+first_rec)],'v1')
    config.curr_line_no = line_read+first_rec  # Store current line number

    # Process line and extract content into components (name, geocode, etc.)
    #
    [name_comp, geocode_comp, locality_comp, date1_comp, date2_comp] = \
           inout.process_line(line)

    # Make a local empty working copy of the output field dictionary  - - - - -
    #
    output_fields = config.output_field.copy()
    output_fields_keys = output_fields.keys()
    for k in output_fields_keys:
      output_fields[k] = ''  # Set all fields to an empty string

    # Standardise name component  - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (type(name_comp) == types.ListType):  # Givenname and surname separate

      givenname_comp = name_comp[0].strip()
      surname_comp   = name_comp[1].strip()

      if (givenname_comp != ''):  # There is a givenname  - - - - - - - - - - -

        inout.log_message('  Givenname component: |'+givenname_comp+'|','v1')

        givenname_comp = name.clean_name_component(givenname_comp)
        [name_list, tag_list] = name.tag_name_component(givenname_comp)
        output_fields['gender_guess'] = name.get_gender_guess(name_list, \
                                        tag_list)
        [name_list, tag_list, output_fields['title']] = \
                                         name.get_title(name_list, tag_list)

        [output_fields['givenname'], output_fields['alt_givenname']] = \
                       name.get_name_component(name_list, tag_list, 'gname')

      if (surname_comp != ''):  # There is a surname  - - - - - - - - - - - - -

        inout.log_message('  Surname component: |'+surname_comp+'|','v1')

        surname_comp = name.clean_name_component(surname_comp)
        [name_list, tag_list] = name.tag_name_component(surname_comp)
        [output_fields['surname'], output_fields['alt_surname']] = \
                        name.get_name_component(name_list, tag_list, 'sname')

    elif (name_comp.strip() != ''):  # Given- and surname both in one field - -

      inout.log_message('  Name component: |'+name_comp+'|','v1')

      name_comp = name.clean_name_component(name_comp)
      [name_list, tag_list] = name.tag_name_component(name_comp)

      output_fields['gender_guess'] = name.get_gender_guess(name_list,tag_list)

      [name_list, tag_list, output_fields['title']] = \
                                        name.get_title(name_list, tag_list)

      if (config.name_standard_method == 'rules'):
        name_dict = name.get_names_rules(name_list, tag_list, 'gname')

      elif (config.name_standard_method == 'hmm'):
        name_dict = name.get_names_hmm(name_list, tag_list)

      else:
        inout.log_message('Illegal name standardisation method:'+ \
                          config.name_standard_method,'err')
        raise Exception()

      for (field,value) in name_dict.items(): # Assign to output dictionary
          output_fields[field] = value 

    # Standardise geocode and locality components using HMM - - - - - - - - - -
    #
    if (config.geoloc_standard_method == 'hmm') and \
       ((geocode_comp.strip() != '') or (locality_comp.strip() != '')):

      geoloc_comp = geocode_comp.strip()+' '+locality_comp.strip()
      inout.log_message('  Geocode and locality component: |'+geoloc_comp+'|',\
                        'v1')

      geoloc_comp = locality.clean_geoloc_component(geoloc_comp)
      [geoloc_words, geoloc_tags] = locality.tag_geoloc_component(geoloc_comp)

      if (geoloc_words != []):  # Component not empty, do HMM standardisation

        geoloc_dict = locality.get_geoloc_hmm(geoloc_words,geoloc_tags)

        for (field,value) in geoloc_dict.items(): # Assign to output dictionary
          output_fields[field] = value

    # Standardise geocode component using rules - - - - - - - - - - - - - - - -
    #
    elif (config.geoloc_standard_method == 'rules') and \
         (geocode_comp.strip() != ''):
      inout.log_message('  Geocode component: |'+geocode_comp+'|','v1')

      ### TO BE DONE
      inout.log_message('Rules based standardisation for geocode is' + \
                        'not implemented yet','err')
      raise Exception()

    # Standardise locality component using rules  - - - - - - - - - - - - - - -
    #
    elif (config.geoloc_standard_method == 'rules') and \
         (locality_comp.strip() != ''):
      inout.log_message('  Locality component: |'+locality_comp+'|','v1')

      ### TO BE FINALISED
      inout.log_message('Rules based standardisation for locality is' + \
                        'not implemented yet','err')
      raise Exception()

#      locality_comp = locality.clean_geoloc_component(locality_comp)
#      [loc_words, loc_tags] = locality.tag_geoloc_component(locality_comp)
#
#      [terr,loc_words2,loc_tags2] = locality.get_territory(loc_words,loc_tags)
#      if (terr != ''):
#        output_fields['territory'] = terr
#
#      [pc,loc_words3,loc_tags3] = locality.get_postcode(loc_words2,loc_tags2)
#      if (pc != ''):
#        output_fields['postcode'] = pc
#
#      [loc_name, loc_quali, loc_words4, loc_tags4] = \
#         locality.get_localityname_qualifier(loc_words3, loc_tags3)
#      if (loc_name != ''):
#        output_fields['locality_name'] = loc_name
#      if (loc_quali != ''):
#        output_fields['locality_quali'] = loc_quali
#
#      if (loc_words4 != []):  # Not all words are standardised yet
#        print '  # Remaining word list:', loc_words4  ###### TEST
#        print '  # Remaining tag list: ', loc_tags4   ###### TEST

    # Standardise date strings  - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (date1_comp != ''):
      inout.log_message('  Date1 component: |'+date1_comp+'|','v1')

      [day1,month1,year1,status1] = date.parse_datestr(date1_comp)
      if (day1 != -1):
        output_fields['day1'] = str(day1)
      if (month1 != -1):
        output_fields['month1'] = str(month1)
      if (year1 != -1):
        output_fields['year1'] = str(year1)

    if (date2_comp != ''):
      inout.log_message('  Date2 component: |'+date2_comp+'|','v1')

      [day2,month2,year2,status2] = date.parse_datestr(date2_comp)
      if (day2 != -1):
        output_fields['day2'] = str(day2)
      if (month2 != -1):
        output_fields['month2'] = str(month2)
      if (year2 != -1):
        output_fields['year2'] = str(year2)

    # Create log message of output fields - - - - - - - - - - - - - - - - - - -
    #
    msg = ['  Standardised record output fields:']
    for (field,value) in output_fields.items():
      if (value != '') and (value != []):
        msg.append('    '+field+':'+str(value))
    inout.log_message(msg,'v1')

    # Save standardised record into output field
    #
    out_line = inout.compose_line(output_fields)
    f_out.write(out_line+os.linesep)

    # Increment line counter and go to beginning of loop  - - - - - - - - - - -
    #
    line_read += 1

    inout.log_message('','v1')  # Print empty lines between records

  # Close files - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  f_in.close()
  f_out.close()

  msg = ['','Number of warnings: '+str(config.num_warning), \
         'Number of corrected word spillings: '+str(config.num_word_spills)]
  inout.log_message(msg,'v1')

  print msg[1]
  print msg[2]

  inout.log_message('End.','v1')


# ----------------------------------------------------------------------------

# Start main standardisation routine
#
standard()

# ----------------------------------------------------------------------------
