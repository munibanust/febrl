README for generate.py      June 2003
----------------------

The example data sets in this directory were created using the
'generate.py' database generator and the frequency files given
in this directory. For more information on 'generate.py' please
see Chapter 'Auxiliary Programs' in the Febrl manual.

1) dataset1.csv

   python generate.py dataset1.csv 500 500 1 uniform

   This data set contains 1000 records (500 original and 500
   duplicates, with exactly one duplicate per original record.


2) dataset2.csv

   python generate.py dataset2.csv 4000 1000 5 poisson

   This data set contains 5000 records (4000 originals and 1000
   duplicates), with a maximum of 5 duplicates based on one original
   record (and a poisson distribution of duplicate records).
   Distribution of duplicates:
     19 originals records have 5 duplicate records
     47 originals records have 4 duplicate records
    107 originals records have 3 duplicate records
    141 originals records have 2 duplicate records
    114 originals records have 1 duplicate record
    572 originals records have no duplicate record


3) dataset3.csv

   python generate.py dataset3.csv 2000 3000 5 zipf

   This data set contains 5000 records (2000 originals and 3000
   duplicates), with a maximum of 5 duplicates based on one original
   record (and a Zipf distribution of duplicate records).
   Distribution of duplicates:
    168 originals records have 5 duplicate records
    161 originals records have 4 duplicate records
    212 originals records have 3 duplicate records
    256 originals records have 2 duplicate records
    368 originals records have 1 duplicate record
   1835 originals records have no duplicate record


4) dataset4a.csv and dataset4b.csv

   python generate.py dataset4.csv 5000 5000 1 uniform

   Generated as one data set with 10000 records (5000 originals and
   5000 duplicates, with one duplicate per original), the originals
   have been split from the duplicates, into

   dataset4a.csv (containing the 5000 original records)
   dataset4b.csv (containing the 5000 duplicate records)

   These two data sets can be used for testing linkage procedures.
