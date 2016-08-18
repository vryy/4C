# Script for comparing result of the parallel and serial post processing by the post_drt_ensight filter
#
# Maintainer: A. Nagler

from paraview import servermanager as sm
import sys
import csv
import operator

def isEqual(ref,line,tol):
  if len(ref) != len(line):
    print 'check your reference csv-file'
    sys.exit(1)
  for j in range(0,len(ref)):
    if not (float(line[j]) >= (float(ref[j])-float(tol))) and (float(line[j]) <= (float(ref[j])+float(tol))):
      return False
  return True


if sm.vtkSMProxyManager.GetVersionMajor() <= 3:
  sm.Connect()

  # read in of parallel calculated case file
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


  # read in of serial calculated case file
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

  # sorting of data according to xyz coordinates of the points
  # sortlist* are lists with values of the corresponding csv-files. Each list contains further lists for each line in the csv-file.
  # In scalar transport problems involving scatra-scatra interface coupling, several nodes may be located at
  # exactly the same position. For the sorting to still be unique, we include the concentration values 'phi_1'
  # or 'c_1', and possibly also the flux components 'flux_phi_1:{0;1;2}', as additional sorting criteria.
  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1'), head1.index('flux_domain_phi_1:0'), head1.index('flux_domain_phi_1:1'), head1.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1'), head1.index('flux_boundary_phi_1:0'), head1.index('flux_boundary_phi_1:1'), head1.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1')))
  elif 'c_1' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('c_1')))
  else:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2))

  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1'), head2.index('flux_domain_phi_1:0'), head2.index('flux_domain_phi_1:1'), head2.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1'), head2.index('flux_boundary_phi_1:0'), head2.index('flux_boundary_phi_1:1'), head2.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1')))
  elif 'c_1' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('c_1')))
  else:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2))

  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1'), head_ref.index('flux_domain_phi_1:0'), head_ref.index('flux_domain_phi_1:1'), head_ref.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1'), head_ref.index('flux_boundary_phi_1:0'), head_ref.index('flux_boundary_phi_1:1'), head_ref.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1')))
  elif 'c_1' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('c_1')))
  else:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2))

  # remove last element (tolerance) in each list in sortlist_ref
  tol = []
  for i in range(0,len(sortlist_ref)):
    index_tol = len(sortlist_ref[i]) - 1
    tol.append(sortlist_ref[i][index_tol])
    sortlist_ref[i] = sortlist_ref[i][0:index_tol]

  # comparison
  for ele1, ele2 in zip(sortlist1, sortlist2):
    if ele1 != ele2:
      print 'csv lists are not equal for: '
      print ele1
      print ele2
      sys.exit(1)
  print 'files are identical'

  # verification of results
  for i in range(0,len(sortlist_ref)):
    for j in range(0,len(sortlist1)):
      if isEqual(sortlist_ref[i],sortlist1[j],tol[i]):
        break
      if j == (len(sortlist1)-1):
        print 'results in csv-files are NOT correct'
        print sortlist_ref[i],'in',sys.argv[3],'is not being found in xxx_par.csv'
        sys.exit(1)
  print'results in csv-files are correct'

else:
  sm.Connect()

  # read in of parallel calculated case file
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


  # read in of serial calculated case file
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
  index2  = head2.index('Points:0')
  
  # read in of reference file
  ref_sortlist_name = sys.argv[3]
  csv_reader_ref = csv.reader(open(ref_sortlist_name,'r'), delimiter = ',')
  head_ref   = csv_reader_ref.next()
  index_ref  = head_ref.index('Points:0')

  # sorting of data according to xyz coordinates of the points
  # sortlist* are lists with values of the corresponding csv-files. Each list contains further lists for each line in the csv-file.
  # In scalar transport problems involving scatra-scatra interface coupling, several nodes may be located at
  # exactly the same position. For the sorting to still be unique, we include the concentration values 'phi_1'
  # or 'c_1', and possibly also the flux components 'flux_phi_1:{0;1;2}', as additional sorting criteria.
  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1'), head1.index('flux_domain_phi_1:0'), head1.index('flux_domain_phi_1:1'), head1.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1'), head1.index('flux_boundary_phi_1:0'), head1.index('flux_boundary_phi_1:1'), head1.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('phi_1')))
  elif 'c_1' in head1:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2, head1.index('c_1')))
  else:
    sortlist1 = sorted(csv_reader1, key = operator.itemgetter(index1, index1 + 1, index1 + 2))

  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1'), head2.index('flux_domain_phi_1:0'), head2.index('flux_domain_phi_1:1'), head2.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1'), head2.index('flux_boundary_phi_1:0'), head2.index('flux_boundary_phi_1:1'), head2.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('phi_1')))
  elif 'c_1' in head2:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2, head2.index('c_1')))
  else:
    sortlist2 = sorted(csv_reader2, key = operator.itemgetter(index2, index2 + 1, index2 + 2))

  if 'phi_1' and 'flux_domain_phi_1:0' and 'flux_domain_phi_1:1' and 'flux_domain_phi_1:2' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1'), head_ref.index('flux_domain_phi_1:0'), head_ref.index('flux_domain_phi_1:1'), head_ref.index('flux_domain_phi_1:2')))
  elif 'phi_1' and 'flux_boundary_phi_1:0' and 'flux_boundary_phi_1:1' and 'flux_boundary_phi_1:2' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1'), head_ref.index('flux_boundary_phi_1:0'), head_ref.index('flux_boundary_phi_1:1'), head_ref.index('flux_boundary_phi_1:2')))
  elif 'phi_1' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('phi_1')))
  elif 'c_1' in head_ref:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2, head_ref.index('c_1')))
  else:
    sortlist_ref = sorted(csv_reader_ref, key = operator.itemgetter(index_ref, index_ref + 1, index_ref + 2))

  # remove last element (tolerance) in each list in sortlist_ref
  tol = []
  for i in range(0,len(sortlist_ref)):
    index_tol = len(sortlist_ref[i]) - 1
    tol.append(sortlist_ref[i][index_tol])
    sortlist_ref[i] = sortlist_ref[i][0:index_tol]

  # comparison
  for ele1, ele2 in zip(sortlist1, sortlist2):
    if ele1 != ele2:
      print 'csv lists are not equal for: '
      print ele1
      print ele2
      sys.exit(1)
  print 'files are identical'

  # verification of results
  for i in range(0,len(sortlist_ref)):
    for j in range(0,len(sortlist1)):
      if isEqual(sortlist_ref[i],sortlist1[j],tol[i]):
        break
      if j == (len(sortlist1)-1):
        print 'results in csv-files are NOT correct'
        print sortlist_ref[i],'in',sys.argv[3],'is not being found in xxx_par.csv'
        sys.exit(1)
  print'results in csv-files are correct'