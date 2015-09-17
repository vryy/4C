# Script for comparing result of the parellel and serial post processing by the post_drt_ensight filter
#
# Maintainer: A. Nagler

from paraview import servermanager as sm
import sys
import csv
import operator

if sm.vtkSMProxyManager.GetVersionMajor() <= 3:
  sm.Connect()

  # read in of parallel calculatet case file
  ens_reader1 = sm.sources.ensight()
  #p1 = ['displacement', '1']
  #p2 = ['nodal_EA_strains_xyz', '1']
  #p3 = ['nodal_cauchy_stresses_xyz', '1']

  c1 = ['element_EA_strains_xyz', '0']
  c2 = ['element_cauchy_stresses_xyz', '0']

  ens_reader1.CaseFileName = sys.argv[1]
  #ens_reader1.PointArrayStatus = p1+p2+p3
  ens_reader1.CellArrayStatus = c1+c2

  ens_reader1.UpdatePipelineInformation()

  converter1 = sm.filters.AttributeDataToTableFilter(Input=ens_reader1, FieldAssociation="Points", AddMetaData=1)

  csvWriter1 = sm.writers.CSVWriter(Input=converter1)
  csvWriter1.FileName = './xxx_par.csv'
  csvWriter1.UpdatePipeline()


  csv_reader1 = csv.reader(open('./xxx_par0.csv','r'), delimiter = ',')
  head1   = csv_reader1.next()
  index1  = head1.index('Points:0')


  # read in of serial calculatet case file
  ens_reader2 = sm.sources.ensight()

  ens_reader2.CaseFileName = sys.argv[2]
  #ens_reader2.PointArrayStatus = p1+p2+p3
  ens_reader2.CellArrayStatus = c1+c2

  ens_reader2.UpdatePipelineInformation()

  converter2 = sm.filters.AttributeDataToTableFilter(Input=ens_reader2, FieldAssociation="Points", AddMetaData=1)

  csvWriter2 = sm.writers.CSVWriter(Input=converter2)
  csvWriter2.FileName = './xxx_ser.csv'
  csvWriter2.UpdatePipeline()

  csv_reader2 = csv.reader(open('./xxx_ser0.csv', 'r'), delimiter = ',')
  head2   = csv_reader2.next()
  index2  = head1.index('Points:0')

  # read in of reference file
  ref_sortlist_name = sys.argv[3]
  csv_reader_ref = csv.reader(open(ref_sortlist_name,'r'), delimiter = ',')
  head_ref   = csv_reader_ref.next()
  index_ref  = head_ref.index('Points:0')

  # sorting of data according to xyz coodinates of the Points
  sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2))
  sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2))
  sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2))

  # comparison
  for ele1, ele2 in zip(sortlist1, sortlist2):
      if ele1 != ele2:
	  print 'csv lists are not equal for: '
	  print ele1
	  print ele2
	  sys.exit(1)
  print 'files are identical'

  # verification of results
  for ref_line in sortlist_ref:
    if not ref_line in sortlist1:
      print 'results in csv lists are NOT correct'
      sys.exit(1)
  print'results in csv lists are correct'

else:
  sm.Connect()

  # read in of parallel calculatet case file
  ens_reader1 = sm.sources.EnSightReader()
  #p1 = ['displacement', '1']
  #p2 = ['nodal_EA_strains_xyz', '1']
  #p3 = ['nodal_cauchy_stresses_xyz', '1']

  c1 = ['element_EA_strains_xyz', '0']
  c2 = ['element_cauchy_stresses_xyz', '0']

  ens_reader1.CaseFileName = sys.argv[1]
  #ens_reader1.PointArrayStatus = p1+p2+p3
  ens_reader1.CellArrays = c1+c2

  ens_reader1.UpdatePipelineInformation()

  converter1 = sm.filters.ExtractAttributes(Input=ens_reader1, FieldAssociation="Points", AddMetaData=1)

  csvWriter1 = sm.writers.CSVWriter(Input=converter1)
  csvWriter1.FileName = './xxx_par.csv'
  csvWriter1.UpdatePipeline()


  csv_reader1 = csv.reader(open('./xxx_par0.csv','r'), delimiter = ',')
  head1   = csv_reader1.next()
  index1  = head1.index('Points:0')


  # read in of serial calculatet case file
  ens_reader2 = sm.sources.EnSightReader()

  ens_reader2.CaseFileName = sys.argv[2]
  #ens_reader2.PointArrayStatus = p1+p2+p3
  ens_reader2.CellArrays = c1+c2

  ens_reader2.UpdatePipelineInformation()

  converter2 = sm.filters.ExtractAttributes(Input=ens_reader2, FieldAssociation="Points", AddMetaData=1)

  csvWriter2 = sm.writers.CSVWriter(Input=converter2)
  csvWriter2.FileName = './xxx_ser.csv'
  csvWriter2.UpdatePipeline()

  csv_reader2 = csv.reader(open('./xxx_ser0.csv', 'r'), delimiter = ',')
  head2   = csv_reader2.next()
  index2  = head1.index('Points:0')

  # read in of reference file
  ref_sortlist_name = sys.argv[3]
  csv_reader_ref = csv.reader(open(ref_sortlist_name,'r'), delimiter = ',')
  head_ref   = csv_reader_ref.next()
  index_ref  = head_ref.index('Points:0')

  # sorting of data according to xyz coodinates of the Points
  sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2))
  sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2))
  sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2))

  # comparison
  for ele1, ele2 in zip(sortlist1, sortlist2):
      if ele1 != ele2:
	  print 'csv lists are not equal for: '
	  print ele1
	  print ele2
	  sys.exit(1)
  print 'files are identical'

  # verification of results
  for ref_line in sortlist_ref:
    if not ref_line in sortlist1:
      print 'results in csv lists are NOT correct'
      sys.exit(1)
  print'results in csv lists are correct'
