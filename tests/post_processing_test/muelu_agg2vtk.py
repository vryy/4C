#
# Maintainer: Rui Fang
# fang@lnm.mw.tum.de
# http://www.lnm.mw.tum.de/
# 089 - 289-15251
#
import argparse
import os
import sys
import vtk
sys.path.append(os.path.join(os.path.dirname(__file__),'../code_test'))
from read_ccarat_NIGHTLYTESTCASES import read_ccarat
from vtk import *

def main(argv=None):

  # create parser for input arguments
  parser = argparse.ArgumentParser(description='Python script for converting [prefix]aggs_level%LEVELID_block%BLOCKID_proc%PROCID.out files with aggregation information from MueLu into Paraview-readable *.vtp files')

  # add input arguments to parser
  parser.add_argument('DAT',help='*.dat input file')
  parser.add_argument('PATH',help='*.out and *.vtp file path')
  parser.add_argument('NUMLEVEL',help='number of multigrid levels',type=int)
  parser.add_argument('NUMPROC',help='number of processors which generated *.out files',type=int)
  parser.add_argument("-p","--prefix",help='*.out and *.vtp file name prefix',default='')
  parser.add_argument("-s","--scaling",help='geometric scaling factor',type=float,default=1.)

  # read input arguments
  arguments = parser.parse_args()

  # read *.dat input file
  sectiontitles,sections = read_ccarat(arguments.DAT)

  # extract *.out and *.vtp file path
  filepath = os.path.abspath(arguments.PATH)

  # extract number of multigrid levels
  numlevel = arguments.NUMLEVEL
  if(numlevel < 1):
    raise ValueError("Invalid number of multigrid levels!")

  # extract number of processors which generated *.out files
  numproc = arguments.NUMPROC
  if(numproc < 1):
    raise ValueError("Invalid number of processors which generated *.out files!")

  # extract *.out and *.vtp file name prefix
  prefix = arguments.prefix

  # extract geometric scaling factor
  scaling = arguments.scaling

  # extract problem dimension
  dimension = -1
  for iline in sections['PROBLEM SIZE']:
    if iline[0] == 'DIM':
      dimension = int(iline[1])
      break
  if dimension == -1:
    raise ValueError("Couldn't extract problem dimension from *.dat input file!")
  if dimension != 3:
    raise ValueError("Script only works for 3D geometries so far!")

  # extract spatial coordinates of all nodes from *.dat input file
  nodescoords = sections['NODE COORDS']

  # initialize dictionary for spatial coordinates of nodes associated with single matrix blocks
  blockids2nodecoords = {}
  for iblock in sections['DESIGN S2I COUPLING VOL CONDITIONS / PARTITIONING']:
    if iblock[0] == 'E':
      blockids2nodecoords[int(iblock[1])] = {}

  # fill dictionary by looping over all nodes from *.dat input file
  for inode in sections['DVOL-NODE TOPOLOGY']:

    # extract spatial coordinates of current node
    inodecoords = nodescoords[int(inode[1]) - 1]

    # add scaled spatial coordinates of current node to dictionary
    blockids2nodecoords[int(inode[3])][int(inode[1])] = (float(inodecoords[3]) * scaling,float(inodecoords[4]) * scaling,float(inodecoords[5]) * scaling)

  # remove offset of node IDs associated with each matrix block to obtain minimum node ID of zero
  for iblock in blockids2nodecoords:

    # compute offset for current matrix block
    offset = min(blockids2nodecoords[iblock])

    # loop over all node IDs associated with current matrix block
    for inode in sorted(blockids2nodecoords[iblock]):

      # subtract offset from current node ID
      blockids2nodecoords[iblock][inode - offset] = blockids2nodecoords[iblock].pop(inode)

  # process all multigrid levels
  for ilevel in range(0,numlevel):

    # screen output
    print ""
    print "  Processing multigrid level " + str(ilevel) + " ..."

    # initialize offset for aggregate IDs
    aggid_offset = 0

    # process all matrix blocks associated with current multigrid level
    for iblock in blockids2nodecoords:

      # screen output
      print "  |"
      print "  | Processing matrix block " + str(iblock - 1) + " ..."

      # initialize dictionaries for nodes and processors associated with single node aggregates in current matrix block of current multigrid level
      aggids2nodeids = {}
      aggids2procids = {}

      # loop over all processors
      for iproc in range(0,numproc):

        # assemble path of *.out file associated with current processor, multigrid level, and matrix block
        outfilepath = filepath + "/" + prefix + "aggs_level" + str(ilevel) + "_block" + str(iblock - 1) + "_proc" + str(iproc) + ".out"

        # check existence of *.out file
        if os.path.exists(outfilepath):

          # screen output
          print "  | | Reading output file " + outfilepath + " ...",

          # loop over all aggregates, i.e., over all lines of current *.out file
          for iline in file(outfilepath):

            # remove whitespace characters from beginning and end of current line
            line = iline.strip()

            # check syntax of current line
            if line.find("Agg ") == 0:

              # split current line into two parts before and after colon symbol
              line = line.split(": ")

              # extract aggregate ID from first part of current line
              agginfo = line[0].split(" ")
              aggid = int(agginfo[1])

              # add IDs of nodes associated with current aggregate to dictionary
              aggids2nodeids[aggid] = map(int,line[1].split(" "))

              # add ID of processor associated with current aggregate to dictionary
              aggids2procids[aggid] = int(agginfo[3])

            else:
              raise ValueError("Output file " + outfilepath + " contains line with invalid syntax!")

          print "Done!"

        else:
          raise ValueError("Output file " + outfilepath + " doesn't exist!")

      # initialize data set for *.vtp file
      vtpdata = vtk.vtkAppendPolyData()

      # loop over all aggregates in current matrix block of current multigrid level
      for aggid,nodeids in aggids2nodeids.iteritems():

        # initialize data set for points associated with current aggregate
        points = vtk.vtkPoints()

        # loop over all nodes associated with current aggregate
        for inode in range(0,len(nodeids)):

          # extract spatial coordinates of current node
          nodecoords = blockids2nodecoords[iblock][nodeids[inode]]

          # add current node to data set
          points.InsertNextPoint(nodecoords[0],nodecoords[1],nodecoords[2])

        # perform Delaunay triangulation of points associated with current aggregate
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        delaunay = vtk.vtkDelaunay3D()
        delaunay.SetInputData(polydata)
        delaunay.Update()

        # create surface filter for Delaunay triangulation
        surfacefilter = vtk.vtkDataSetSurfaceFilter()
        surfacefilter.SetInputConnection(delaunay.GetOutputPort())
        surfacefilter.Update()

        # create integer array for aggregate IDs associated with triangulated points
        aggids = vtk.vtkUnsignedIntArray()
        aggids.SetNumberOfComponents(1)
        aggids.SetName("Aggregate ID")
        agggid = aggid + aggid_offset
        for ipoint in range(0,points.GetNumberOfPoints()):
          aggids.InsertNextTuple1(agggid)

        # create integer array for processor IDs associated with triangulated points
        procids = vtk.vtkUnsignedCharArray()
        procids.SetNumberOfComponents(1)
        procids.SetName("Processor ID")
        procid = aggids2procids[aggid]
        for ipoint in range(0,points.GetNumberOfPoints()):
          procids.InsertNextTuple1(procid)

        # add surface filter for Delaunay triangulation to data set for *.vtp file
        polydata = surfacefilter.GetOutput()
        polydata.GetPointData().SetScalars(aggids)
        polydata.GetPointData().AddArray(procids)
        vtpdata.AddInputData(polydata)

        # add points associated with current aggregate to data set for *.vtp file
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.GetPointData().SetScalars(aggids)
        polydata.GetPointData().AddArray(procids)
        vtpdata.AddInputData(polydata)

      # finalize data set for *.vtp file
      vtpdata.Update()

      # assemble path of *.vtp file associated with current multigrid level and matrix block
      vtpfilepath = filepath + "/" + prefix + "aggs_level" + str(ilevel) + "_block" + str(iblock - 1) + ".vtp"

      # screen output
      print "  | | Creating result file " + vtpfilepath + " ...",

      # write *.vtp file associated with current multigrid level and matrix block
      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName(vtpfilepath)
      writer.SetInputData(vtpdata.GetOutput())
      writer.Write()

      # screen output
      print "Done!"

      # determine root nodes, i.e., centroids of aggregates in current matrix block of current multigrid level, for aggregates in current matrix block of next coarser multigrid level
      if numlevel - ilevel > 1:

        # screen output
        print "  | | Calculating root nodes for aggregation on next coarser multigrid level ... ",

        # initialize dictionary for spatial coordinates of root nodes
        rootnodescoords = {}

        # loop over all aggregates in current matrix block of current multigrid level
        for aggid in range(len(aggids2nodeids)):

          # extract IDs of nodes associated with current aggregate
          nodeids = aggids2nodeids[aggid]

          # extract number of nodes associated with current aggregate
          numnode = len(nodeids)

          # initialize list for spatial coordinates of root node associated with current aggregate
          rootnodecoords = [0.] * dimension

          # loop over all nodes associated with current aggregate
          for inode in nodeids:

            # extract spatial coordinates of current node
            nodecoords = blockids2nodecoords[iblock][inode]

            # add spatial coordinates of current node to list
            for idim in range(dimension):
              rootnodecoords[idim] += nodecoords[idim]

          # compute and store spatial coordinates of centroid associated with current aggregate
          rootnodescoords[aggid] = tuple([rootnodecoords[idim] / len(nodeids) for idim in range (dimension)])

        # replace spatial coordinates of nodes associated with current matrix block of current multigrid level by spatial coordinates of root nodes
        blockids2nodecoords[iblock] = rootnodescoords

        # screen output
        print "Done!"

      # update offset for aggregate IDs
      aggid_offset += len(aggids2nodeids)

      # screen output
      print "  | ... Done!"

  # final linebreak
  print ""

if __name__ == "__main__":
  sys.exit(main())
