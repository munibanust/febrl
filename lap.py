# =============================================================================
# lap.py - Routines for linear assignment procedures
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
# The Original Software is "lap.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module lap.py - Routines for linear assignment procedures

   Currently, the following lap algorithm is implemented:

   - auction

     Implementation of the asymmtric auction algortihm by Bertsekas as
     described in:

     "Auction Algorithms for Network Flow Problems: A Tutorial Introduction"
     Dimitri P. Bertsekas, Computational Optimization and Applications, Vol.
     1, pp. 7-66, 1992.
"""

# =============================================================================
# The following flag can be set to True or False
# It is used to test if various internal data structures and intermediate
# results are correct. It should be set to False for normal use.

DO_TESTS = False  # Perform various tests for immediate results

SAVE_PARALLEL_TEST_FILES = False  # Write intermediate results to file for
                                  # inspection

MAX_ROW_ELEMENTS = 10  # Following an idea by William Winkler (US Census
                       # Bureau) only the elements with the largest weights
                       # in a row (plus the diagonal element) are used for the
                       # linear assignment procedure.
                       # Winkler suggests a value of 5, but 10 seems to be good
                       # as well.
                       # If you want all elements to be kept, set this to a
                       # very large number.
                       # Note that setting MAX_ROW_ELEMENTS to a small number
                       # (less 20 or so) result in a changed assignment

# =============================================================================
# Imports go here

import time
import os  # For testing only

import output
import parallel

# =============================================================================

def do_lap(lap_method, results_dict, process_type, threshold):
  """Linear sum assignment procedure.

     This routine calculates a linear assignments for one-to-one matching.

     The routine do_lap does all kinds of preprocessing, including the
     extraction of unique record pairs (which can be removed before the lap is
     applied) and the extraction of sub-set which can be solved independently.
     These sub-sets are then given to the chosen lap routine.

     The routine takes as input a results dictionary, as produced by a
     classifier (see classification.py), and returns a dictionary with the
     assigned record pair numbers as keys and the corresponding weight as
     values.

     Possible methods are 'auction'

     The process_type attribute must either be set to 'deduplication' or to
     'linkage' in order to be able to preprocess the classifier data prior to
     the lap procedure.
  """

  if (lap_method not in ['auction']):
    print 'error:Illegal method for lap method: %s' % (str(lap_method))
    raise Exception

  if (process_type not in ['deduplication', 'linkage']):
    print 'error:Illegal value for attribute "process_type": %s' %\
          (str(process_type))
    raise Exception

  if (results_dict == {}):
    print 'error:Empty results dictionary'
    raise Exception

  lap_start_time = time.time()  # Start timer

  lap_results = {}  # Result dictionary with final record pair matches

  # Make one (or two) disctionary of all assigned rercord numbers
  #
  if (process_type == 'deduplication'):
    used_rec_nums = {}
  else:
    used_rec_nums_a = {}
    used_rec_nums_b = {}

  # Make sure the threshold is a number if it is defined
  #
  if (threshold != None):
    if (not (isinstance(threshold, int) or isinstance(threshold, float))):
      print 'error:Threshold is not a number: %s' % (str(threshold))
      raise Exception

  print '1:  Start linear assignment procedure using method: %s' % (lap_method)
  print '1:    Original length of results dictionary: %i' % (len(results_dict))

  # Step 1: Filter out record pairs with weight lower than the threshold  - - -
  #
  if (threshold != None):
    print '1:    Remove record pairs with weight less than: %f' % (threshold)
  else:
    threshold = -999999999999.999  # Make it a very very small number

  work_dict = {}  # Make an empty working dictionary

  for row_num in results_dict:  # Loop over all record numbers (keys)
    row_dict = results_dict[row_num]  # Get corresponding record dictionary

    new_row_dict = {}  # Start a new record dictionary

    for col_num in row_dict:  # Loop over all records in this dictionary
      weight = row_dict[col_num]

      if (weight >= threshold):
        new_row_dict[col_num] = weight  # Copy to new dictionary

    if (new_row_dict != {}):  # Only insert non empty dictionaries
      work_dict[row_num] = new_row_dict

  results_len = len(work_dict)  # Save results length (after filtering)

  if (threshold > -999999999999.999):
    print '1:    Length of working dictionary after filtering: %i' % \
          (results_len)

  # Step 2: Remove all matches (record pairs) which are unique  - - - - - - - -
  #         (i.e. which don't have matches with other records)
  #
  row_num_dict = {}  # Count occurences of record numbers in rows
  col_num_dict = {}  # Count occurences of record numbers in columns

  for row_num in work_dict:  # First count occurences of rows and columns

    # Insert a count for the row number
    #
    row_num_dict[row_num] = row_num_dict.get(row_num, 0) + 1

    row_dict = work_dict[row_num]

    for col_num in row_dict:

      # Increase a count for a column number
      #
      col_num_dict[col_num] = col_num_dict.get(col_num, 0) + 1

      if (process_type == 'deduplication'):

        # For deduplication, insert symmetric record numbers as well
        #
        row_num_dict[col_num] = row_num_dict.get(col_num, 0) + 1
        col_num_dict[row_num] = col_num_dict.get(row_num, 0) + 1

  for row_num in work_dict.keys():  # Secondly remove unique rows and column

    row_dict = work_dict[row_num]  # Get corresponding record dictionary

    if (len(row_dict) == 1):  # Only one record pair for this record

      col_num, weight = row_dict.items()[0]  # Get the only element in row

      if (row_num_dict[row_num] == 1) and (col_num_dict[col_num] == 1):

        #################### START TEST CODE ##################################

        if (DO_TESTS == True):
          if (process_type == 'deduplication'):
            if (row_num in used_rec_nums):
              print 'warning:Record number %i already used for deduplication' \
                    % (row_num)
            if (col_num in used_rec_nums):
              print 'warning:record number %i already used for deduplication' \
                    % (col_num)
          else:
            if (row_num in used_rec_nums_a):
              print 'warning:Record number A %i already used for linkage' % \
                    (row_num)
            if (col_num in used_rec_nums_b):
              print 'warning:record number B %i already used for linkage' % \
                    (col_num)

        #################### END TEST CODE ##################################

        lap_results[(row_num,col_num)] = True  # Insert into final results
        del work_dict[row_num]  # And delete the record in the results

        if (process_type == 'deduplication'):
          used_rec_nums[row_num] = True
          used_rec_nums[col_num] = True
        else:
          used_rec_nums_a[row_num] = True
          used_rec_nums_b[col_num] = True

  print '1:    Found and extracted %i unique record ' % (len(lap_results)) + \
        'pairs in results dictionary'

  for rec_pair in lap_results:
    print '3:      %s' % (str(rec_pair))
  print '3:'

  lap_pair_extract_time = time.time() - lap_start_time

  #################### START PARALLEL TEST CODE ###########################

  if (SAVE_PARALLEL_TEST_FILES == True):
    tmp_list = lap_results.items()
    tmp_list.sort()
    f = open('one2one-unique-dedup-'+str(parallel.rank())+'-'+ \
             str(parallel.size()),'w')
    for c in tmp_list:
      f.write(str(c)+os.linesep)
    f.close()

    tmp_list = work_dict.keys()
    tmp_list.sort()
    f = open('work-dict-'+str(parallel.rank())+'-'+str(parallel.size()),'w')
    for c in tmp_list:
      cc = work_dict[c].items()
      cc.sort()
      f.write(str(c)+':: '+str(cc)+os.linesep)
    f.close()

  #################### END PARALLEL TEST CODE #############################

  #################### START TEST CODE ########################################
  # Test if a record only appears once in the lap results dictionary
  #
  if (DO_TESTS == True):
    if (process_type == 'deduplication'):
      test_dict = {}
      for (rec_a, rec_b) in lap_results:
        if (test_dict.has_key(rec_a)):
          print 'warning:Record %s is already in the test dictionary' % \
                (str(rec_a))
        else:
          test_dict[rec_a] = 1
        if (test_dict.has_key(rec_b)):
          print 'warning:Record %s is already in the test dictionary' % \
                (str(rec_b))
        else:
          test_dict[rec_b] = 1

    else:  # Linkage process
      test_dict_a = {}
      test_dict_b = {}
      for (rec_a, rec_b) in lap_results:
        if (test_dict_a.has_key(rec_a)):
          print 'warning:Record %s is already in test dictionary A' % \
                (str(rec_a))
        else:
          test_dict_a[rec_a] = 1
        if (test_dict_b.has_key(rec_b)):
          print 'warning:Record %s is already in test dictionary B' % \
                (str(rec_b))
        else:
          test_dict_b[rec_b] = 1

  #################### END TEST CODE ##########################################

  if (len(work_dict) == 0):  # All record pairs are processed - - - - - - - - -
    return lap_results

  print '1:    Remaining number of records in working dictionary: %i' % \
        (len(work_dict)) + ' (down from: %i)' % (results_len)

  # Step 3: Find connected sub-sets in the results dictionary - - - - - - - - -
  #         (using depth-first search)
  #
  visited =  {}  # Dictionary which will contain all so far visited rows
  sub_sets = {}  # Dictionary which will contain the sub-sets extracted 

  print '1:    Find connected sub-graphs in results dictionary'

  lap_subset_start_time = time.time()

  max_sub_set_length = -1
  num_visited = 0  # Number of rows visited so far
  row_num_done = 0
  work_dict_len = len(work_dict)

  work_dict_rows = work_dict.keys()
  work_dict_rows.sort()

  # Create a column oriented work dictionary  - - - - - - - - - - - - - - - - -
  #
  col_work_dict = {}

  for row_num in work_dict_rows:  # Loop over all rows
    row_dict = work_dict[row_num]

    for col_num in row_dict:
      col_dict = col_work_dict.get(col_num,{})
      col_dict[row_num] = True  # Only position is needed, but not the weight
      col_work_dict[col_num] = col_dict

  for row_num in work_dict_rows:  # Loop over all rows

    if (not visited.has_key(row_num)):  # This row has not been visited yet

      visited[row_num] = row_num  # Mark visited as 'seeding' row
      num_visited += 1
      print '2:      Create sub-set with seeding record %i' % (row_num)

      process_queue = [row_num]  # Start a new queue of rows to process
      row_sub_set = {row_num:1}  # Row numbers connected to this row

      while (process_queue != []): # Process rows until all connected rows done
        print '3:        Process queue: %s' % (str(process_queue))

        next_row = process_queue.pop(0)  # Get and remove first row to process
        row_col_numbers = work_dict[next_row].keys()  # Get columns in this row

        # For deduplication, also insert row number into this column numbers
        #
        if (process_type == 'deduplication'):
          row_col_numbers.append(next_row)

        print '3:          Row %i with column numbers: %s' % \
              (next_row, str(row_col_numbers))

        # Get the row numbers from all column numbers
        #
        for col_num in row_col_numbers:

          # Get list of all row numbers in this column
          #
          row_num_dict = col_work_dict.get(col_num, {})
          row_num_list = row_num_dict.keys()

          if (process_type == 'deduplication') and (col_num in work_dict) and \
            (col_num not in row_num_list):
            row_num_list.append(col_num)

          print '3:          Column: %i with row numbers: %s' % \
                (col_num, str(row_num_list))

          for row_num2 in row_num_list:
            row_sub_set[row_num2] = 1
            if (not visited.has_key(row_num2)):  # Check if it's a new row
              process_queue.append(row_num2)
              print '3:          Appended row number %i to process queue' % \
                    (row_num2)

              visited[row_num2] = row_num  # Mark row as visited by seeding row
              num_visited += 1
              print '3:          Row %i connected to row %i' % \
                    (row_num2, row_num)

      sub_sets[row_num] = row_sub_set.keys()  # Only store keys

      if (len(row_sub_set) > max_sub_set_length):
        max_sub_set_length = len(row_sub_set)

      print '3:        Sub-set contains records: %s' % \
            (str(row_sub_set.keys()))

    row_num_done += 1

    # Now determine timing and print progress report (every 10%)  - - - - - - -
    # (only if more than 100 records in the work dictionary)
    #
    if (work_dict_len >= 100) and (row_num_done % int(work_dict_len/10) == 0):
      used_time =     time.time() - lap_subset_start_time
      perc_done =     100.0 * row_num_done / work_dict_len
#      todo_time =     (work_dict_len - row_num_done) * \  #################
#                      (used_time / row_num_done)
      todo_time =     (work_dict_len - num_visited) * \
                      (used_time / row_num_done)

      used_time_string =       output.time_string(used_time)
      todo_time_string =       output.time_string(todo_time)

      print '1:      Processed %.1f%% of records in %s (%i/%i records ' % \
            (perc_done, used_time_string, num_visited, work_dict_len) + \
            'visited)'
      print '1:        Estimated %s until finished' % (todo_time_string)

  del col_work_dict  # Delete the column oriented work dictionary

  num_sub_sets = len(sub_sets)  # Get the total number of sub-sets

  lap_subset_total_time = time.time() - lap_subset_start_time
  lap_subset_total_time_string = output.time_string(lap_subset_total_time)

  print '1:    Extracted %i sub-sets in %s' % \
        (num_sub_sets, lap_subset_total_time_string)
  print '1:      Longest sub-set contains %i rows' % (max_sub_set_length)

  #################### START TEST CODE ########################################
  # Test if all the sub-sets are mutually exclusive, and if the seed rows are
  # in the sub-set record lists
  #
  if (DO_TESTS == True):
    for seed_row in sub_sets:
      row_list = sub_sets[seed_row]
      if (seed_row not in row_list):
        print 'warning:Seed row %s not in sub-set row list: %s' % \
              (str(seed_rec), str(row_list))
      for rec_num in row_list:
        for seed_row2 in sub_sets:
          row_list2 = sub_sets[seed_row2]
          if (seed_row != seed_row2):  # Don't test itself
            if (rec_num in row_list2):
              print 'warning:Record %s in more than one sub-set: %s, %s' % \
                    (str(rec_num), str(row_list), str(row_list2))

  #################### END TEST CODE ##########################################

  #################### START PARALLEL TEST CODE ###########################

  if (SAVE_PARALLEL_TEST_FILES == True):
    tmp_list = sub_sets.keys()
    tmp_list.sort()
    f = open('sub-sets-'+str(parallel.rank())+'-'+ \
             str(parallel.size()),'w')
    for s in tmp_list:
      tmp_sub_set = sub_sets[s]
      tmp_sub_set.sort()

      f.write(str(s)+'::'+str(tmp_sub_set)+os.linesep)
    f.close()

  #################### END PARALLEL TEST CODE #############################

  # Now loop over all sub-sets  - - - - - - - - - - - - - - - - - - - - - - - -
  # (pre-process them first before giving them to the actual linear assignment
  # method)
  #
  lap_lap_start_time = time.time()
  lap_comm_time = 0.0

  sub_set_cnt = 0  # A round robin counter, used for parallelism

  sub_set_rows = sub_sets.keys()
  sub_set_rows.sort()  # Needed to make the same on all processes

  for seed_row in sub_set_rows:

    # Distribute sub-sets equally to all processors
    #
    if ((sub_set_cnt % parallel.size()) == parallel.rank()):

      row_list = sub_sets[seed_row]
      row_list.sort()

      print '1:'
      print '1:    Sub-set %i of %i with seed row %i contains %i rows' % \
            (sub_set_cnt, num_sub_sets, seed_row, len(row_list))
      print '3:      Sub-set rows:  %s' % (str(row_list))

      if (len(row_list) == 1):  # Special case: One row only  - - - - - - - - -
        max_weight = -99999.9
        max_col = -1
        row_dict = work_dict[row_list[0]]  # Get the dictionary for this row

        # Find element with largest weight
        #
        for col_num in row_dict:
          weight = row_dict[col_num]
          if (weight > max_weight):
            max_weight = weight
            max_col = col_num

        # Assignment dictionary is of form col_num:row_num
        #
        tmp_assign_dict = {max_col:row_list[0]}  # Make record pair dictionary

        print '2:      Special case sub-set with one row only, ' + \
              'assignment pair: (%i,%i)' % (row_list[0], max_col)

      else:  # General case with more than one row  - - - - - - - - - - - - - -

        # Get minimal and maximal weights, and lists with row and column
        # numbers
        #
        min_weight =  999999.9
        max_weight = -999999.9
        col_numbers = {}
        row_col_numbers = {}

        for row_num in row_list:  # Loop over rows in this sub-set
          row_dict = work_dict[row_num]  # Get the dictionary for this row
          row_col_numbers[row_num] = 1

          for col_num in row_dict:
            weight = row_dict[col_num]
            col_numbers[col_num] = 1
            row_col_numbers[col_num] = 1

            if (weight < min_weight):
              min_weight = weight
            if (weight > max_weight):
              max_weight = weight

        print '3:      Minimal and maximal weight: %.3f / %.3f' % \
              (min_weight, max_weight)

        row_numbers = work_dict.keys()
        col_numbers = col_numbers.keys()
        row_numbers.sort()
        col_numbers.sort()
        num_rows = len(row_numbers)
        num_cols = len(col_numbers)

        row_col_numbers = row_col_numbers.keys()
        row_col_numbers.sort()

        #print '1:      Row numbers:    %s' % (str(row_numbers))
        #print '1:      Column numbers: %s' % (str(col_numbers))
        #print '1:      Row/colum numbers: %s' % (str(row_col_numbers))
        #print '1:      Number of unique weights: %i' % (len(weight_dict))

        # Deal with the special case that there is only one column number - - -
        #
        if (num_cols == 1):
          max_weight = -99999.9
          max_row = -1

          col_num = col_numbers[0]  # Get the column number

          # Find element with largest weight
          #
          for row_num in row_list:  # Loop over rows

            # Get only weight in row
            #
            row_weight = work_dict[row_num].values()[0]

            if (row_weight > max_weight):
              max_weight = row_weight
              max_row = row_num

          # Assignment dictionary is of form col_num:row_num
          #
          tmp_assign_dict = {col_num:max_row}  # Make record pair dictionary

          print '2:      Special case sub-set with one column only, ' + \
                'assignment pair: (%i,%i)' % (max_row, col_num)

        else:  # General case with more than one row and column - - - - - - - -

          # Construct the cost dictionary - - - - - - - - - - - - - - - - - - -
          #
          cost_dict = {}
          dim = len(row_col_numbers)  # Final dimension of the LAP

          min_cost = -max_weight * (dim + 1)  # Use original wieghts

          for row_num in row_list:  # Loop over rows

            row_dict = work_dict[row_num]

            # Get the column numbers in this row
            #
            col_list = row_dict.keys()
            col_list.sort()

            row_cost_dict = cost_dict.get(row_num, {})

            for col_num in col_list:
              weight = row_dict[col_num]

              cost = weight * (dim + 1)  # Use original weights
              row_cost_dict[col_num] = cost  # And store into row dictionary

              # Insert symmetric element as well (if not on diagonal)
              #
              if (row_num != col_num):
                row_cost_dict2 = cost_dict.get(col_num,{})

                if (row_num not in row_cost_dict2):  # Only insert if not there

                  if (process_type == 'deduplication'):
                    row_cost_dict2[row_num] = cost  # Insert symmetric cost
                  else:  # Linkage process
                    row_cost_dict2[row_num] = min_cost  # Insert minimal cost

                  # And insert diagonal element if there is none
                  #
                  if (not row_cost_dict2.has_key(col_num)):
                    row_cost_dict2[col_num] = min_cost

                  cost_dict[col_num] = row_cost_dict2

            # Make sure there is a diagonal element (for feasibility)
            #
            if (not row_cost_dict.has_key(row_num)):
              row_cost_dict[row_num] = min_cost

            # If more than MAX_ROW_ELEMENTS elements in row only take the
            # largest (following an idea by William Winkler)
            #
            if (len(row_cost_dict) > MAX_ROW_ELEMENTS):
              row_col_numbers = row_cost_dict.keys()
              row_weights =     row_cost_dict.values()
              row_elem_list = map(None, row_weights, row_col_numbers)
              row_elem_list.sort()

              diag_weight = row_cost_dict[row_num]  # Keep diagonal element
              row_cost_dict = {row_num:diag_weight}

              ##print '   '
              ##print '****** row_elem_list: %s' % (str(row_elem_list))
              ##print '   '

              for (weight, col_num) in row_elem_list[-MAX_ROW_ELEMENTS:]:
                row_cost_dict[col_num] = weight

            # Insert row into cost dictionary
            #
            cost_dict[row_num] = row_cost_dict

          # Get the final row and column numbers  - - - - - - - - - - - - - - -
          #
          row_numbers = cost_dict.keys()
          col_numbers = {}
          for row_dict in cost_dict.values():
            col_numbers.update(row_dict)

          col_numbers = col_numbers.keys()
          row_numbers.sort()
          col_numbers.sort()

          # Check if number of rows and columns are equal - - - - - - - - - - -
          #
          if (len(row_numbers) != len(col_numbers)):
            print 'error:Different number of rows (%i) and columns (%i)' \
                  % (len(row_numbers), len(col_numbers))
            raise Exception

          print '1:      Cost dictionary with %i rows/columns given to ' % \
                (len(row_numbers)) + 'assignment method %s:' % (lap_method)
          print '2:        Row numbers:    %s' % (str(row_numbers))
          print '2:        Column numbers: %s' % (str(row_numbers))
          print '2:        Minimal weight: %3f' % (min_weight)
          print '2:        Maximal weight: %3f' % (max_weight)
          print '3:        Cost dictionary: %s' % (str(cost_dict))
          print '3:        Process type:    %s' % (process_type)

          #################### START PARALLEL TEST CODE #######################

          if (SAVE_PARALLEL_TEST_FILES == True):

            tmp_list = cost_dict.keys()
            tmp_list.sort()
            tmp_str = str(sub_set_cnt)+':: '+ str(min_weight) + ' / ' + \
                      str(max_weight) + ', ' + str(row_numbers) + ', ' + \
                      process_type + '::'
            for k in tmp_list:
              tmp_list2 = cost_dict[k].items()
              tmp_list2.sort()
              tmp_str = tmp_str + ' ' + str(tmp_list2) + ' / '

            f = open('lap-calling-'+str(parallel.rank())+'-'+ \
                str(parallel.size()),'a')
            f.write(tmp_str+os.linesep)
            f.close()

          #################### END PARALLEL TEST CODE #########################

          # Call the lap method which returns an assignment dictionary  - - - -
          #
          if (lap_method == 'auction'):
            tmp_assign_dict = auction(cost_dict, min_weight, max_weight, \
                                  row_numbers, col_numbers)
          else:
             print 'error:LAP method %s not implemented' % (lap_method)
             raise Exception

      # If run in parallel, send temporary assignment dictionary process 0  - -
      #
      if (parallel.rank() > 0):
        tmp_time = time.time()
        parallel.send(tmp_assign_dict, 0)
        lap_comm_time += (time.time() - tmp_time)
        print '1:      Sent assignment dictionary with %i entries to process' \
              % (len(tmp_assign_dict)) + ' 0'

    # Only process 0 inserts temporary assignment dictionary into results - - -
    #
    if (parallel.rank() == 0):

      # Receive assignment dictionary from other process if necessary
      #
      p = (sub_set_cnt % parallel.size())  # Process number to receive from

      if (p != 0):
        tmp_time = time.time()
        tmp_assign_dict = parallel.receive(p)
        lap_comm_time += (time.time() - tmp_time)
        print '1:    Received subset %i of %i assignment dictionary with ' % \
              (sub_set_cnt, num_sub_sets) + '%i entries from process %i' % \
              (len(tmp_assign_dict), p)

      # Post-process the assignment dictionary  - - - - - - - - - - - - - - - -
      #
      assign_pairs = {}

      for rec_num_b in tmp_assign_dict:
        rec_num_a = tmp_assign_dict[rec_num_b]

        # Now check if this record pair is in the original results dictionary
        #
        if (rec_num_a in results_dict):
          row_dict = results_dict[rec_num_a]
          if (rec_num_b in row_dict):
            weight = row_dict[rec_num_b]

            # Insert into dictionary of potential record pairs
            #
            assign_pairs[(rec_num_a, rec_num_b)] = weight

      #################### START PARALLEL TEST CODE ###########################

      if (SAVE_PARALLEL_TEST_FILES == True):

        tmp_list = tmp_assign_dict.items()
        tmp_list.sort()
        tmp_list2 = assign_pairs.items()
        tmp_list2.sort()
        tmp_list3 = sub_sets[seed_row]
        tmp_list3.sort()

        f = open('assignments-'+str(parallel.rank())+'-'+ \
             str(parallel.size()),'a')
        f.write(str(sub_set_cnt)+', '+str(seed_row)+os.linesep)
        f.write(str(tmp_list)+os.linesep)
        f.write(str(tmp_list2)+os.linesep)
        f.write(str(tmp_list3)+os.linesep)
        f.write(os.linesep)
        f.close()

      #################### END PARALLEL TEST CODE #############################

      # Sort the assigned pairs according to their weight
      #
      assign_weights = assign_pairs.values()  # Get the weights in a list
      assign_rec_pairs = assign_pairs.keys()  # And the record pairs

      assign_pair_list = map(None, assign_weights, assign_rec_pairs)
      assign_pair_list.sort()

      num_assigned_pairs = 0  # Number of assigned pairs for this sub-set
      dedup_check_rec_nums = {}  # Already assigned record numbers

      while (assign_pair_list != []):  # Now check all record pairs
        check_pair = assign_pair_list.pop()  # Get largest weight record pair

        weight = check_pair[0]
        rec_num_a = check_pair[1][0]
        rec_num_b = check_pair[1][1]
        rec_pair = (rec_num_a, rec_num_b)

        # Now check if a record pair has already been used in an assignment
        # and for a deduplication process also check if any of the two
        # records has been used in an assignment
        #
        if ((process_type == 'linkage') and \
            (rec_num_a not in used_rec_nums_a) and \
            (rec_num_b not in used_rec_nums_b)) or \
           ((process_type == 'deduplication') and \
            (rec_num_a not in used_rec_nums) and \
            (rec_num_b not in used_rec_nums)):

          # For deduplication insert record numbers into used record numbers
          #
          if (process_type == 'deduplication'):
            used_rec_nums[rec_num_a] = True
            used_rec_nums[rec_num_b] = True
          else:
            used_rec_nums_a[rec_num_a] = True
            used_rec_nums_b[rec_num_b] = True

          if (rec_pair not in lap_results):
            lap_results[rec_pair] = True
            num_assigned_pairs += 1
          else:
            print 'warning:Record pair (%i,%i) already in LAP results' \
                  % (rec_num_a, rec_num_b)

      print '2:      Inserted %i (out of %i) record pairs into LAP ' % \
            (num_assigned_pairs, len(tmp_assign_dict)) + 'results'

    sub_set_cnt += 1

    # Report progress every 10% (only if more than 100 sub-sets)  - - - - - - -
    #
    if (num_sub_sets >= 100) and (sub_set_cnt % int(num_sub_sets / 10) == 0):
        used_time =    time.time() - lap_lap_start_time
        perc_done =    100.0 * sub_set_cnt / num_sub_sets
        sub_set_time = used_time / sub_set_cnt
        todo_time =    (num_sub_sets - sub_set_cnt) * sub_set_time

        used_time_string =    output.time_string(used_time)
        todo_time_string =    output.time_string(todo_time)
        sub_set_time_string = output.time_string(sub_set_time)

        print '1:      Processed %.1f%% (%i/%i) of sub-sets in %s' % \
                (perc_done, sub_set_cnt, num_sub_sets, used_time_string) + \
                ' (%s per sub-set)' % (sub_set_time_string)
        print '1:        Estimated %s until finished' % (todo_time_string)

  print '1:  Total number of assignments: %i' % (len(lap_results))
  print '1:    Number of rows in original results dictionary: %i' % \
        (len(results_dict))

  #################### START TEST CODE ########################################
  # Test if a record only appears once in the lap results dictionary
  #
  if (DO_TESTS == True) and (parallel.rank() == 0):
    if (process_type == 'deduplication'):
      test_dict = {}
      for (rec_a, rec_b) in lap_results:
        if (test_dict.has_key(rec_a)):
          print 'warning:Record %i is already in the test dictionary' % \
                (rec_a)+' rec_pair: (%i,%i)' % (rec_a, rec_b)
        else:
          test_dict[rec_a] = True
        if (test_dict.has_key(rec_b)):
          print 'warning:Record %i is already in the test dictionary' % \
                (rec_b)+' rec_pair: (%i,%i)' % (rec_a, rec_b)
        else:
          test_dict[rec_b] = 1

    else:  # Linkage process
      test_dict_a = {}
      test_dict_b = {}
      for (rec_a, rec_b) in lap_results:
        if (test_dict_a.has_key(rec_a)):
          print 'warning:Record %s is already in test dictionary A' % \
                (str(rec_a))
        else:
          test_dict_a[rec_a] = 1
        if (test_dict_b.has_key(rec_b)):
          print 'warning:Record %s is already in test dictionary B' % \
                (str(rec_b))
        else:
          test_dict_b[rec_b] = 1

  #################### END TEST CODE ##########################################

  lap_stop_time = time.time()
  lap_lap_time = lap_stop_time - lap_lap_start_time
  lap_total_time = lap_stop_time - lap_start_time

  lap_pair_extract_time_string = output.time_string(lap_pair_extract_time)
  lap_subset_total_time_string = output.time_string(lap_subset_total_time)
  lap_lap_time_string =          output.time_string(lap_lap_time)
  if (parallel.size() > 1):
    lap_comm_time_string =         output.time_string(lap_comm_time)
  lap_total_time_string =        output.time_string(lap_total_time)

  print '1:'
  print '1:  Finished linear record pair assignment procedure'
  print '1:    Time for extracting unique record pairs: %s' % \
        (lap_pair_extract_time_string)
  print '1:    Time for creating record sub-sets:       %s' % \
        (lap_subset_total_time_string)
  print '1:    Time for linear assignment algorithm:    %s' % \
        (lap_lap_time_string)
  if (parallel.size() > 1):
    print '1:    Time for communication:                  %s' % \
          (lap_comm_time_string)
  print '1:    Total time for linear assignment:        %s' % \
        (lap_total_time_string)
  print '1:'

  return lap_results

# =============================================================================

def auction(cost_dict, min_cost, max_cost, row_numbers, col_numbers):
  """Linear sum assignment procedure based on symetric auction algorithm.

     Re-implementation of a FORTRAN code, taken from: 

     http://web.mit.edu/afs/athena.mit.edu/user/d/i/dimitrib/www/auction.txt

     and as described in:

     "Auction Algorithms for Network Flow Problems: A Tutorial Introduction"
     Dimitri P. Bertsekas, Computational Optimization and Applications, Vol.
     1, pp. 7-66, 1992.

     and other papers by Dimitri P. Bertsekas, see:

     http://www.mit.edu:8001//people/dimitrib/publ.html

     Takes as input a cost dictionary, the minimum and maximum weights in this
     dictionary, and two sorted lists with the row and column numbers.
  """

  #################### START PARALLEL TEST CODE ###############################

  if (SAVE_PARALLEL_TEST_FILES == True):

    tmp_list = cost_dict.keys()
    tmp_list.sort()
    tmp_str = str(min_cost) + ' /'  + str(max_cost)+', '+ str(row_numbers) + \
              ' / ' + str(col_numbers) + ':: '
    for k in tmp_list:
      tmp_list2 = cost_dict[k].items()
      tmp_list2.sort()
      tmp_str = tmp_str + ' ' + str(tmp_list2) + ' / '

    f = open('lap-auction-'+str(parallel.rank())+'-'+str(parallel.size()),'a')
    f.write(tmp_str+os.linesep)
    f.close()

    #################### END PARALLEL TEST CODE ###############################

  i_large = 100000000  # A value larger than max_cost

  auction_results = {}  # Result dictionary with final record pair matches

  num_rows = len(row_numbers)
  num_cols = len(col_numbers)

  if (num_rows != num_cols):
    print 'error:Asymmetric problem given to symmetric "auction" algorithm'
    raise Exception

  # Set parameters for auction algorithm  - - - - - - - - - - - - - - - - - - -
  #
  if (max_cost > int(i_large / (num_rows + 1))):
    print 'error:Cost range too large to work with integer epsilon'
    raise Exception

  max_cost *= (num_rows + 1)

  beg_eps =    max_cost / 5   # Maybe smaller, can even be 1, but not smaller
  end_eps =    num_rows / 10  # Must be smaller than beg_eps, can be 1
  if (end_eps < 1):
    end_eps = 1
  elif (end_eps > beg_eps):
    end_eps = beg_eps
  factor =     5              # Must be greater than 1
  start_incr = beg_eps / 10   # Maybe even be 1, but not smaller
  if (start_incr < 1):
    start_incr = 1

  # Initialisation
  #
  eps =         beg_eps
  i_small =     -i_large
  large_incr =  int(i_large / 10)
  thresh =      min(int(num_rows / 5), 100)  # Maximal value 100
  incr_factor = 2
  cycles =      1
  average =     num_cols
  num_phases =  1

  # Initialise dictionaries for prices and row assignments
  #
  pcol =     {}
  assigned = {}

  for col_num in col_numbers:
    pcol[col_num] =     i_small
    # assigned[col_num] = -1#######

  # Initialise list of un-assigned rows (all rows at the beginning)
  #
  list =     row_numbers[:]
  no_list =  num_rows

  do_phase = True  # Set flag so a first phase is performed

  # Start sub problem (scaling phase with new epsilon)  - - - - - - - - - - - -
  #
  while (do_phase == True):

    if (eps == 1):
      thresh = 0
    incr = max(start_incr, eps)  # Increment must not be larger than epsilon

    print '2:      Start of a scaling phase with epsilon: %f' % (eps)

    do_cycle = True  # Set flag so a first cycle is performed
    cycle_count = 0

    while (do_cycle == True):

      # Start forward auction cycle - - - - - - - - - - - - - - - - - - - - - -
      #
      no_new_list = 0  # Initialise count of next list of un-assigned rows

      # Cycle through the current list of un-assigned rows
      #
      for i in range(no_list):
        row_num =  list[i]
        row_dict = cost_dict[row_num]
        row_list = row_dict.items()
        row_list.sort()  ################ Maybe not needed ??? #######
        row_len = len(row_dict)

        # Get first and second column number and cost
        #
        col_num, cost =   row_list[0]
        col_num2, cost2 = row_list[1]

        max1 = cost -  pcol[col_num]
        max2 = cost2 - pcol[col_num2]
        if (max1 > max2):
          best_col_num = col_num
        elif (max1 < max2):  # Swap maximum values
          max1, max2 = max2, max1
          best_col_num = col_num2
        else:  # Both are the same
          if (col_num < col_num2):  # Make sure best column number is smallest
            best_col_num = col_num
          else:
            best_col_num = col_num2

        if (row_len > 2):  # Row has more than two elements
          for c in range(2, row_len):
            col_num3, cost3 = row_list[c] # Loop through cols

            max_tmp = cost3 - pcol[col_num3]
            if (max_tmp > max2):
              if (max_tmp > max1):
                best_col_num = col_num3  # New best column
                max2 = max1
                max1 = max_tmp
              elif (max_tmp == max1) and (col_num3 < best_col_num):
                  best_col_num = col_num3  # Best column number is smallest
              else:
                max2 = max_tmp

        # Row bids for best column increasing its price, and gets - - - - - -
        # assigned to best column, while any row assigned to best
        # column gets un-asssigned
        #
        pcol[best_col_num] = pcol[best_col_num] + max1 - max2 + incr

        old_row = assigned.get(best_col_num, -1)  #########
        assigned[best_col_num] = row_num
        if (old_row >= 0):  # Row has been assigned
          list[no_new_list] = old_row  # Save un-assigned row into list
          no_new_list += 1

      cycle_count += 1
      if ((cycle_count % 10000) == 0):
        print '1:        Finished %i cycles, %i out of %i rows un-assigned' % \
              (cycle_count, no_new_list, num_rows)
        if (no_new_list > 0):
          print '1:          Un-assigned rows: %s' % (str(list[:no_new_list]))

      # Collect statistics
      #
      average = (cycles * average + no_list) / (cycles+1)
      cycles += 1

      # Check if there are still 'many' unassignedrows, i.e. if the - - - - - -
      # number of unassignd rows is greater than the parameter 'thresh'.
      # If not, replace current list with the new list and go for another
      # cycle. Otherwise, if epsilon > 1, reduce epsilon, reset the
      # asignment to empty and restart auction.
      # If epsilon == 1 terminate.
      # Also increase the minimal bidding increment up to a maximun value
      # of epsilon (this is the adaptive feature)
      #
      incr *= incr_factor
      if (incr > eps):
        incr = eps
      if (no_new_list > thresh):
        no_list = no_new_list
      else:
        do_cycle = False  # Set flag so cycle is left

    # End of sub-problem (scaling phase) - - - - - - - - - - - - - - - - - - -
    #
    if (eps == 1):
      do_phase = False  # Set flag so phase is left
    else:
      num_phases += 1
      eps = int(eps / factor)
      if (eps > incr):
        eps = int(eps / factor)
      if (eps < 1) or (eps < end_eps):
        eps = 1
      thresh = int(thresh / factor)

      print '1:        End of a scaling phase, new epsilon: %3f' % (eps)

      t_min = min(pcol.values()+[i_large])

      col_num_list = assigned.keys()   ##########
      col_num_list.sort()              ##########

      for col_num in assigned:    ###########
#      for col_num in col_num_list:
        row_num = assigned[col_num]
        if (row_num >= 0):
          list[no_new_list] = row_num
          no_new_list += 1
          assigned[col_num] = -1

      incr = t_min - i_small  # Reset minimum price to i_small
      for col_num in col_numbers:  # Update all prices
        pcol[col_num] = pcol[col_num] - incr

      # Final parameter updates before starting another scaling phase
      #
      no_list = no_new_list

      if (start_incr < eps):
        start_incr *= factor

  return assigned

# =============================================================================
