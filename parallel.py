# =============================================================================
# parallel.py - Module for parallel imports and definitions.
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
# The Original Software is "parallel.py".
# The Initial Developers of the Original Software are Dr Peter Christen
# (Department of Computer Science, Australian National University) and Dr Tim
# Churches (Centre for Epidemiology and Research, New South Wales Department
# of Health). Copyright (C) 2002, 2003 the Australian National University and
# others. All Rights Reserved.
# Contributors:
#
# =============================================================================

"""Module parallel.py - Module for parallel imports and definitions.

   Written by Ole M. Nielsen, January 2003, ANU/MSI
   Extended by Peter Christen, February 2003

   Use PyPar for parallelism if it is installed, otherwise define a
   rudimentary interface for sequential execution.
"""

# =============================================================================

# Set the mode for printing ('host' or 'all')
#
printmode = 'host'  # 'all'

# Set mode for saving (writing) files (data sets) ('host' or 'all')
#
writemode = 'host'  # 'all'

# -----------------------------------------------------------------------------
# Conditional import of PyPar
#
try:
  import pypar

# -----------------------------------------------------------------------------
# PyPar is not installed, so define sequential interface for parallel functions
#
except:
  print 'warning:Could not import module "PyPar", defining sequential '+ \
        'interface'
  def size(): return 1
  def rank(): return 0

  def Get_processor_name():  # - - - - - - - - - - - - - - - - - - - - - - - -
    import os
    try:
      hostname = os.environ['HOST']
    except:
      try:
        hostname = os.environ['HOSTNAME']
      except:
        hostname = 'Unknown'

    return hostname

  def Abort():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import sys
    sys.exit()

  def Finalize():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pass

  def Barrier():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pass

  def Wtime():  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import time
    return time.time()

# -----------------------------------------------------------------------------
# PyPar is installed, so import it
#
else:
  from pypar import *

# -----------------------------------------------------------------------------
# Now define processor information for print statements (prompt)

if (size() > 1):
  prompt = 'P%d/%d: ' % (rank(),size())
else:
  prompt = ''  # Empty prompt for sequential runs

# =============================================================================
