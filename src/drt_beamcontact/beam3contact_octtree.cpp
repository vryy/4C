/*!----------------------------------------------------------------------
\file beam3contact.cpp
\brief Octtree for beam contact search

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "beam3contact_octtree.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_Time.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include <map>
#include <math.h>

#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif
#ifdef D_BEAM3II
#include "../drt_beam3ii/beam3ii.H"
#endif

// measure time for octree build and intermediate steps
//#define MEASURETIME

using namespace std;

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/11|
 *----------------------------------------------------------------------*/
Beam3ContactOctTree::Beam3ContactOctTree(ParameterList& params, DRT::Discretization& discret,DRT::Discretization& searchdis, const int& dofoffset):
discret_(discret),
searchdis_(searchdis),
basisnodes_(discret.NumGlobalNodes()),
dofoffset_(dofoffset)
{
  // define max tree depth (maybe, set this as input file parameter)
  maxtreedepth_ = 5;
  // set flag signaling the existence of periodic boundary conditions
	Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
	if(statmechparams.get<double>("PeriodLength",0.0)>0.0)
		periodicBC_ = true;
	else
		periodicBC_ = false;
  // determine bounding box type
  string boundingbox = params.get<string>("BEAMS_OCTREE","None");
  if(boundingbox == "octree_axisaligned")
  {
  	if(!discret_.Comm().MyPID())
  		cout<<"Search routine:\nOctree + Axis Aligned BBs"<<endl;
  	boundingbox_ = Beam3ContactOctTree::axisaligned;
  }
  else if(boundingbox == "octree_cyloriented")
  {
  	if(!discret_.Comm().MyPID())
  		cout<<"Search routine:\nOctree + Cylindrical Oriented BBs"<<endl;
  	boundingbox_ = Beam3ContactOctTree::cyloriented;
  }
  else
  	dserror("No Octree declared in your Input file!");

  // initialize beam diameter
  diameter_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap())));

  // initialize vector mapping bounding boxes to octants with -1.0 for empty
	// initialize with 4 columns (the maximum number of octants a single bounding box can belong to
	bbox2octant_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),4));
	bbox2octant_->PutScalar(-1.0);
	// initialize vector counting the number of shifts across volume boundaries in case of periodic boundary conditions
	// used to optimize bounding box intersection
	if(periodicBC_)
		numshifts_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap()),true));

  return;
}


/*----------------------------------------------------------------------*
 |  calls the almighty Octtree (public)                      meier 01/11|
 *----------------------------------------------------------------------*/
vector<RCP<Beam3contact> > Beam3ContactOctTree::OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions)
{
  // initialize beam diameter
  diameter_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap())));

  // initialize vector mapping bounding boxes to octants with -1.0 for empty
	// initialize with 4 columns (the maximum number of octants a single bounding box can belong to
	bbox2octant_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),4));
	bbox2octant_->PutScalar(-1.0);
	// initialize vector counting the number of shifts across volume boundaries in case of periodic boundary conditions
	// used to optimize bounding box intersection
	if(periodicBC_)
		numshifts_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap()),true));

  //beam diameter
  for(int i=0; i<searchdis_.ElementColMap()->NumMyElements(); i++)
  {
		DRT::Element* beamelement = searchdis_.lColElement(i);
		const DRT::ElementType & eot = beamelement->ElementType();

#ifdef D_BEAM3
		if (eot == DRT::ELEMENTS::Beam3Type::Instance())
			(*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(beamelement))->Izz()) / M_PI));
#endif // #ifdef BEAM3II

#ifdef D_BEAM3II
		if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
			(*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3ii*>(beamelement))->Izz()) / M_PI));
#endif

		// feasibility check
		if ((*diameter_)[i] <= 0.0) dserror("ERROR: Did not receive feasible element radius.");
  }
  // initialize multivector storage of Bounding Boxes
  // (components 0,...,5 contain bounding box limits)
  // (components 6,...,23 contain bounding box limits in case of periodic boundary conditions (2nd part of the box)):
  // a box may b subject to a boundary shift up to 3 times -> 4 segments -> 24 values + 1 bounding box ID
  // (component 24 containts element ID)
  // in case of periodic boundary conditions, we need 12+1 vectors within allbboxes_, without periodic BCs 6+1
  allbboxes_ = Teuchos::null;
  if(periodicBC_)
  	allbboxes_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),4*6+1, true));
  else
  	allbboxes_ = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),7, true));

  // build axis aligned bounding boxes
  CreateBoundingBoxes(currentpositions);

  // call recursive octtree build
  // get the root octant
  LINALG::Matrix<1,6> rootoctantlim = GetRootOctant();
  // clear vector for assigning bounding boxes to octants to be on the safe side before (re)assigning bounding boxes
  std::vector<std::vector<int> > bboxesinoctants;
  bboxesinoctants.clear();
  locateAll(rootoctantlim, bboxesinoctants);

  // intersection checks
  vector<RCP<Beam3contact> > contactpairs;
  BoundingBoxIntersection(currentpositions, bboxesinoctants, &contactpairs);

  return contactpairs;
}// OctTreeSearch()
/*----------------------------------------------------------------------*
 |  Return the octants to which this bounding box belongs               |
 |  (public)                                               mueller 01/11|
 *----------------------------------------------------------------------*/
std::vector<int> Beam3ContactOctTree::InWhichOctantLies(const int& thisBBoxID)
{
	std::vector<int> octants(bbox2octant_->NumVectors(),-1);
	int bboxcolid = searchdis_.ElementColMap()->LID(thisBBoxID);
	for(int i=0; i<bbox2octant_->NumVectors(); i++)
		octants[i] = (int)(*bbox2octant_)[i][bboxcolid];
	return octants;
}
/*----------------------------------------------------------------------*
 |  Intersect the bounding boxes of a certain octant with a given       |
 |  bounding box (public)                                  mueller 01/11|
 *----------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectBBoxesWith(Epetra_SerialDenseMatrix& nodecoords, LINALG::Matrix<2,1>& nodeLID)
{
	/* note:
	 * 1) do not apply this before having constructed the octree. This is merely a query tool
	 * 2)"boxid" does not necessarily coincide with the bounding box we are going to intersect with the other boxes
	 * in the octant. The reason: The bounding box may actually not exist.
	 * Of course, if it does exist, "boxid" will be the id of the bounding box we actually want to check the
	 * other boxes against. However, if the bounding box is merely a hypothetical construct (i.e. there is no actual beam element), then
	 * we have to give a box id that does exist in order to find the correct octant. Ideally, that means that "boxid" should be the ID of
	 * a bounding box which is a direct neighbor of our (hypothetical) bounding box.
	 * 3) nodecoords are the coordinates of the nodes of the (non-)existing element.
	 */
	bool intersection = false;

	// determine bounding box limits
	RCP<Epetra_SerialDenseMatrix> bboxlimits;

	// build bounding box according to given type
	switch(boundingbox_)
	{
		case Beam3ContactOctTree::axisaligned:
			bboxlimits = CreateAABBNoElement(nodecoords);
		break;
		case Beam3ContactOctTree::cyloriented:
			bboxlimits = CreateCOBBNoElement(nodecoords);
		break;
		default: dserror("No or an invalid Octree type was chosen. Check your input file!");
	}

	// retrieve octants in which the bounding box with ID thisBBoxID is located
	std::vector<std::vector<int> > octants;
	octants.clear();
	// get the octants for two bounding boxe (element) GIDs adjacent to each given node LID
	for(int i=0; i<(int)nodeLID.M(); i++)
		octants.push_back(InWhichOctantLies(searchdis_.lColNode((int)nodeLID(i))->Elements()[0]->Id()));

	// intersection of given bounding box with all other bounding boxes in given octant
	for(int ibox=0; ibox<(int)octants.size(); ibox++)
	{
		for(int oct=0; oct<(int)octants[ibox].size(); oct++)
		{
			if(octants[ibox][oct]!=-1)
			{
				for(int i=0; i<bboxesinoctants_->NumVectors(); i++)
				{
					// take only values of existing bounding boxes and not the filler values (-9)
					if((int)(*bboxesinoctants_)[i][octants[ibox][oct]]>-0.9)
					{
						// get the second bounding box ID
						int bboxinoct = (int)(*bboxesinoctants_)[i][octants[ibox][oct]];
						/*check for adjacent nodes: if there are adjacent nodes, then, of course, there
						 * the intersection test will turn out positive. We skip those cases.*/
						// note: bounding box IDs are equal to element GIDs
						bool sharednode = false;
						for(int j=0; j<searchdis_.gElement(bboxinoct)->NumNode(); j++)
						{
							for(int k=0; k<(int)nodeLID.M(); k++)
							{
								if(searchdis_.NodeColMap()->LID(searchdis_.gElement(bboxinoct)->NodeIds()[j])==(int)nodeLID(k))
								{
									sharednode = true;
									break;
								}
							}
							if(sharednode)
								break;
						}
						// apply different bounding box intersection schemes
						if(!sharednode)
						{
							switch(boundingbox_)
							{
								case Beam3ContactOctTree::axisaligned:
									intersection = IsecAABBNoElement(bboxinoct, bboxlimits);
								break;
								case Beam3ContactOctTree::cyloriented:
									intersection = IsecCOBBNoElement(bboxinoct, bboxlimits);
								break;
								default: dserror("No or an invalid Octree type was chosen. Check your input file!");
							}
						}

						if(intersection)
							break;
					}
					else // loop reached the first bogus value (-9)
						break;
				}
			}
			else
				break;
			if(intersection)
				break;
		}
		if(intersection)
			break;
	}

	return intersection;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)                 meier 01/11|
 |  generates bounding boxes extended with factor 1.05                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions)
{
#ifdef MEASURETIME
  double t_AABB = Teuchos::Time::wallTime();
#endif
  // Initialize Variables....................
  //edge length of root box in case of periodic boundary conditions
  double PeriodLength = 0.0;
  if(periodicBC_)
  {
  	Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
		PeriodLength = statmechparams.get<double>("PeriodLength", 0.0);
  }

  // Get Nodes from discretization....................
	// build bounding boxes according to input parameter
	for (int elecolid=0; elecolid<searchdis_.ElementColMap()->NumMyElements(); elecolid++)
	{
		int elegid = searchdis_.ElementColMap()->GID(elecolid);
		// only do stuff for Row Elements
		if(searchdis_.ElementRowMap()->LID(elegid)>-1)
		{
			// get the element with local ID (LID) elecolid
			DRT::Element* element = searchdis_.lColElement(elecolid);

			// vector for the global IDs (GID) of element
			std::vector<int> nodelids(2,0);
			for(int i=0; i<(int)nodelids.size(); i++)
			{
				int gid = element->Nodes()[i]->Id();
				nodelids[i] = searchdis_.NodeColMap()->LID(gid);
			}

			//store nodal positions into matrix coords
			LINALG::SerialDenseMatrix coord(3,2,true);
			for(int i=0; i<(int)nodelids.size(); i++)
			{
				const map<int, LINALG::Matrix<3,1> >::const_iterator nodepos = currentpositions.find(nodelids.at(i));
				for(int j=0; j<coord.M(); j++)
					coord(j,i) = (nodepos->second)(j);
			}

			// build bounding box according to given type
    	switch(boundingbox_)
    	{
    		case Beam3ContactOctTree::axisaligned:
    			CreateAABB(coord, elecolid);
				break;
    		case Beam3ContactOctTree::cyloriented:
    			CreateCOBB(coord, elecolid);
				break;
    		default: dserror("No or an invalid Octree type was chosen. Check your input file!");
    	}
		}
	} //end for-loop which goes through all elements

  // communication of findings
  Epetra_MultiVector allbboxesrow(*(searchdis_.ElementRowMap()),allbboxes_->NumVectors(),true);
  Epetra_Export exporter(*(searchdis_.ElementColMap()),*(searchdis_.ElementRowMap()));
  Epetra_Import importer(*(searchdis_.ElementColMap()),*(searchdis_.ElementRowMap()));
  allbboxesrow.Export(*allbboxes_, exporter, Add);
  allbboxes_->Import(allbboxesrow,importer,Insert);
  if(periodicBC_)
  {
    Epetra_Vector numshiftsrow(*(searchdis_.ElementRowMap()),true);
    numshiftsrow.Export(*numshifts_, exporter, Add);
    numshifts_->Import(numshiftsrow, importer, Insert);
  }

  //Test: print allbboxes_->...................
  //cout << "\n\tTest extendedAABB" << endl;

  /*/Visualization print allbboxes_ to .dat-file and plot with Matlab....................
  std::ostringstream filename;
  filename << "extendedAABB.dat";
  FILE* fp = NULL;
  //open file to write output data into
  fp = fopen(filename.str().c_str(), "w");
  // write output to temporary stringstream;
  std::stringstream myfile;
  for (int u = 0; u < allbboxes_->MyLength(); u++)
  {
    for (int v = 0; v < allbboxes_->NumVectors(); v++)
    {
      myfile <<scientific<<setprecision(10)<<(*allbboxes_)[v][u] <<" ";
    }
    myfile <<endl;
  }
  //write content into file and close it
  fprintf(fp, myfile.str().c_str());
  fclose(fp);*/

#ifdef MEASURETIME
  double bbgentimelocal = Teuchos::Time::wallTime() - t_AABB;
  double bbgentimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&bbgentimelocal, &bbgentimeglobal, 1);

  if(!searchdis_.Comm().MyPID())
    cout << "\n\nBounding Box generation time:\t" << bbgentimeglobal<< " seconds";
#endif

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create an Axis Aligned Bounding Box   (private)                           mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateAABB(const Epetra_SerialDenseMatrix& coord, const int& elecolid)
{
  // factor by which the box is extruded in each dimension (->input file parameter??)
  const double extrusionfactor = 1.05;

  //Decision if periodic boundary condition is applied....................
  //number of spatial dimensions
  const int ndim = 3;

	int elegid = searchdis_.ElementColMap()->GID(elecolid);

  /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
   * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
   * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
   * at PeriodLength);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
   * 0 close to boundary at PeriodLength);*/
  LINALG::Matrix<3,1> cut;
  cut.Clear();

  /* "coord" currently holds the shifted set of coordinates.
   * In order to determine the correct vector "dir" of the visualization at the boundaries,
   * a copy of "coord" with adjustments in the proper places is introduced
   * in unshift, always the second node lies outside of the volume*/
  LINALG::SerialDenseMatrix unshift = coord;

  // Compute "cut"-matrix (only in case of periodic BCs)
  double PeriodLength = 0.0;
  // number of overall shifts
  int numshifts = 0;
  // dof at which the bounding box segment is shifted (used in case of a single shift)
  int shiftdof = -1;
  if(periodicBC_)
  {
		Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
		PeriodLength = statmechparams.get<double>("PeriodLength",0.0);

    for (int dof=0; dof<ndim; dof++)
    {
    	// initialize unshift with coord values
    	unshift(dof,0) = coord(dof,0);
    	unshift(dof,1) = coord(dof,1);
      if (fabs(coord(dof,1)-PeriodLength-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
      {
        cut(dof) = 1.0;
        shiftdof = dof;
        unshift(dof,1) -= PeriodLength;
        numshifts++;
      }
      if (fabs(coord(dof,1)+PeriodLength - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
      {
        cut(dof) = 2.0;
        shiftdof = dof;
        unshift(dof,1) += PeriodLength;
        numshifts++;
      }
    }
    // store number of shifts
  	(*numshifts_)[elecolid] = numshifts;
  }

	/* take action according to number of shifts
	 * this may seem not too elegant, but consider that among the cut elements the majority is only shifted once.
	 * A single shift can be performed with much less computational effort than multiple shifts since in the case of
	 * multiple shifts, one has to reiterate the dof-wise shifts, determine the position of the nodes after each shift
	 * and calculate the coordinates of the found segments.
	 */
	// standard unshifted bounding box
	switch(numshifts)
	{
		case 0:
		{
			//do normal process with nodecoords0 and nodecoords1....................
			//Calculate Center Point of AABB
			LINALG::Matrix<3,1> midpoint;
			for(int i=0; i<(int)midpoint.M(); i++)
					midpoint(i) = 0.5*(coord(i,0) + coord(i,1));

			//Calculate edgelength of AABB
			LINALG::Matrix<3,1> edgelength;
			for(int i=0; i<(int)edgelength.M(); i++)
					edgelength(i) = fabs(coord(i,1) - coord(i,0));

			//Check for edgelength of AABB
			for(int i=0; i<(int)edgelength.M(); i++)
				if (edgelength(i)<(*diameter_)[elecolid])
					edgelength(i) = (*diameter_)[elecolid];

			// Calculate limits of AABB with extrusion around midpoint
			for(int i=0; i<6; i++)
				if(i%2==0)
					(*allbboxes_)[i][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
				else if(i%2==1)
					(*allbboxes_)[i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
			// fill all latter entries  except for the last one (->ID) with bogus values (in case of periodic BCs)
			if(periodicBC_)
				for(int i=6; i<allbboxes_->NumVectors()-1; i++)
					(*allbboxes_)[i][elecolid] = -1e9;

			(*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;
		}
		break;
		default: // broken bounding boxes due to periodic BCs
		{
			Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
			PeriodLength = statmechparams.get<double>("PeriodLength",0.0);

			// directional vector
			LINALG::Matrix<3, 1> dir;
			for (int dof = 0; dof < ndim; dof++)
				dir(dof) = unshift(dof,1) - unshift(dof,0);
			dir.Scale(1.0/dir.Norm2());

			/* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
			 * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
			 * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
			LINALG::Matrix<3,2> lambdaorder;
			lambdaorder.PutScalar(1e6);
			// collect lambdas
			for(int dof=0; dof<(int)lambdaorder.M(); dof++)
			{
				switch((int)cut(dof))
				{
					case 1:
					{
						lambdaorder(dof,0) = -coord(dof, 0) / dir(dof);
						lambdaorder(dof,1) = dof;
					}
					break;
					case 2:
					{
						lambdaorder(dof,0) = (PeriodLength - coord(dof,0)) / dir(dof);
						lambdaorder(dof,1) = dof;
					}
					break;
					default:
					{
						lambdaorder(dof,1) = dof;
					}
				}
			}

			// sort the lambdas (ascending values) and indices accordingly
			// in case of multiple shifts
			if(numshifts>1)
			{
				for(int j=0; j<(int)lambdaorder.M()-1; j++)
					for(int k=j+1; k<(int)lambdaorder.M(); k++)
						if(lambdaorder(k,0)<lambdaorder(j,0))
						{
							double temp = lambdaorder(j,0);
							int tempindex = (int)lambdaorder(j,1);
							lambdaorder(j,0) = lambdaorder(k,0);
							lambdaorder(j,1) = lambdaorder(k,1);
							lambdaorder(k,0) = temp;
							lambdaorder(k,1) = tempindex;
						}
			}
			else	// for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
				for(int n=0; n<(int)lambdaorder.N(); n++)
					lambdaorder(0,n) = lambdaorder(shiftdof,n);

			// calculate segment lambdas
			for(int dof=numshifts-1; dof>0; dof--)
				lambdaorder(dof,0) -= lambdaorder(dof-1,0);

			// the idea is to gradually shift the matrix "unshift" back into the volume and, while doing so,
			// calculate the segments except for the last one
			// determine closest boundary component-wise
			for(int shift=0; shift<numshifts; shift++)
			{
				//second point
				for(int j=0 ;j<unshift.M(); j++)
					unshift(j,1) = unshift(j,0) + lambdaorder(shift,0)*dir(j);

				//Calculate Center Point and edgelength of AABB
				LINALG::Matrix<3,1> midpoint;
				LINALG::Matrix<3,1> edgelength;
				for(int i=0; i<(int)midpoint.M(); i++)
				{
					midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
					edgelength(i) = unshift(i,1) - unshift(i,0);
					//Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
					if (edgelength(i)<(*diameter_)[elecolid])
						edgelength(i) = (*diameter_)[elecolid];
				}
				// Calculate limits of AABB of the current segment (which definitely lies in the volume) with extrusion around midpoint
				for(int i=0; i<6; i++)
					if(i%2==0)
					 (*allbboxes_)[shift*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
					else if(i%2==1)
					 (*allbboxes_)[shift*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;

				int currshift = (int)lambdaorder(shift,1);
				// shift the coordinates of the second point
				if(cut(currshift)==1.0)
					unshift(currshift,1) += PeriodLength;
				else if(cut(currshift)==2.0)
					unshift(currshift,1) -= PeriodLength;
				// make second point the first and calculate new second point in the next iteration!
				for(int j=0; j<unshift.M(); j++)
					unshift(j,0) = unshift(j,1);
			}

			// the last segment
			for(int dof=0; dof<unshift.M(); dof++)
				unshift(dof,1) = coord(dof,1);

			//Calculate Center Point and edgelength of AABB
			LINALG::Matrix<3,1> midpoint;
			LINALG::Matrix<3,1> edgelength;
			for(int i=0; i<(int)midpoint.M(); i++)
			{
				midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
				edgelength(i) = unshift(i,1) - unshift(i,0);
				//Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
				if (edgelength(i)<(*diameter_)[elecolid])
					edgelength(i) = (*diameter_)[elecolid];
			}

			// limits of the last bounding box
			for(int i=0; i<6; i++)
				if(i%2==0)
				 (*allbboxes_)[numshifts*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
				else if(i%2==1)
				 (*allbboxes_)[numshifts*6+i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;

			//fill up the rest of the 24 values with bogus values
			if(numshifts<3)
				for(int i=allbboxes_->NumVectors()-1; i>(numshifts+1)*6-1; i--)
					(*allbboxes_)[i][elecolid] = -1e9;

			 // store GID (=box number)
			 (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;
		}
	}

  //Bring coordinates in case of periodic boundary condition in right order ( "-1e9 signals the bogus value from above)
  //[xmin xmax ymin ymax zmin zmax ...]
  if (periodicBC_)
  {
  	if((*allbboxes_)[6][elecolid]!=-1e9)
  	{
			double minimum =0.0;    double maximum =0.0;

			for(int l=0; l<(allbboxes_->NumVectors()-1)/2;l++)
			{
				// leave loop at first bogus entry
				if((2*l)%6==0 && (*allbboxes_)[2*l][elecolid]==-1e9)
					break;
				minimum = min((*allbboxes_)[2*l][elecolid],(*allbboxes_)[2*l+1][elecolid]);
				maximum = max((*allbboxes_)[2*l][elecolid],(*allbboxes_)[2*l+1][elecolid]);
				//cout << minimum << endl;
				(*allbboxes_)[2*l][elecolid] = minimum;    (*allbboxes_)[2*l+1][elecolid] = maximum;
			}// end of correct
  	}
  } // end of normal AAABB

	return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box   (private)                   mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateCOBB(const Epetra_SerialDenseMatrix& coord, const int& elecolid)
{
	return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create an Axis Aligned Bounding Box without a finite element (private)    mueller 11/11|
 *----------------------------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> Beam3ContactOctTree::CreateAABBNoElement(const Epetra_SerialDenseMatrix& coord)
{
	RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null;

  // factor by which the box is extruded in each dimension (->input file parameter??)
  const double extrusionfactor = 1.05;

  //Decision if periodic boundary condition is applied....................
  //number of spatial dimensions
  const int ndim = 3;

  /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
   * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
   * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
   * at PeriodLength);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
   * 0 close to boundary at PeriodLength);*/
  LINALG::Matrix<ndim,1> cut;
  cut.Clear();

  /* "coord" currently holds the shifted set of coordinates.
   * In order to determine the correct vector "dir" of the visualization at the boundaries,
   * a copy of "coord" with adjustments in the proper places is introduced*/
  LINALG::SerialDenseMatrix unshift = coord;

  // number of shifts for this bounding box
  int numshifts = 0;
  // Compute "cut"-matrix (only in case of periodic BCs)
  double PeriodLength = 0.0;
  // dof at which the bounding box segment is shifted (used in case of a single shift)
	int shiftdof = -1;
	if(periodicBC_)
	{
		Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
		PeriodLength = statmechparams.get<double>("PeriodLength",0.0);

		for (int dof=0; dof<ndim; dof++)
		{
			// initialize unshift with coord values
			unshift(dof,0) = coord(dof,0);
			unshift(dof,1) = coord(dof,1);
			if (fabs(coord(dof,1)-PeriodLength-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
			{
				cut(dof) = 1.0;
				shiftdof = dof;
				unshift(dof,1) -= PeriodLength;
				numshifts++;
			}
			if (fabs(coord(dof,1)+PeriodLength - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
			{
				cut(dof) = 2.0;
				shiftdof = dof;
				unshift(dof,1) += PeriodLength;
				numshifts++;
			}
		}
		bboxlimits = rcp(new Epetra_SerialDenseMatrix((numshifts+1)*6,1));
	}
	else
		bboxlimits = rcp(new Epetra_SerialDenseMatrix(6,1));

	/* take action according to number of shifts
	 * this may seem not too elegant, but consider that among the cut elements the majority is only shifted once.
	 * A single shift can be performed with much less computational effort than multiple shifts since in the case of
	 * multiple shifts, one has to reiterate the dof-wise shifts, determine the position of the nodes after each shift
	 * and calculate the coordinates of the found segments.
	 */
	// standard unshifted bounding box
	switch(numshifts)
	{
		case 0:
		{
			//do normal process with nodecoords0 and nodecoords1....................
			//Calculate Center Point of AABB
			LINALG::Matrix<3,1> midpoint;
			for(int i=0; i<(int)midpoint.M(); i++)
					midpoint(i) = 0.5*(coord(i,0) + coord(i,1));

			//Calculate edgelength of AABB
			LINALG::Matrix<3,1> edgelength;
			for(int i=0; i<(int)edgelength.M(); i++)
					edgelength(i) = fabs(coord(i,1) - coord(i,0));

			//Check for edgelength of AABB
			for(int i=0; i<(int)edgelength.M(); i++)
				if (edgelength(i)<(*diameter_)[0])
					edgelength(i) = (*diameter_)[0];

			// Calculate limits of AABB with extrusion around midpoint
			for(int i=0; i<6; i++)
				if(i%2==0)
					(*bboxlimits)(i,0) = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
				else if(i%2==1)
					(*bboxlimits)(i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
			// fill all latter entries  except for the last one (->ID) with bogus values (in case of periodic BCs)
			if(periodicBC_)
				for(int i=6; i<bboxlimits->M(); i++)
					(*bboxlimits)(i,0) = -1e9;
		}
		break;
		default: // broken bounding boxes due to periodic BCs
		{
			Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
			PeriodLength = statmechparams.get<double>("PeriodLength",0.0);

			// directional vector
			LINALG::Matrix<3, 1> dir;
			for (int dof = 0; dof < ndim; dof++)
				dir(dof) = unshift(dof,1) - unshift(dof,0);
			dir.Scale(1.0/dir.Norm2());

			/* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
			 * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
			 * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
			LINALG::Matrix<3,2> lambdaorder;
			lambdaorder.PutScalar(1e6);
			// collect lambdas
			for(int dof=0; dof<(int)lambdaorder.M(); dof++)
			{
				switch((int)cut(dof))
				{
					case 1:
					{
						lambdaorder(dof,0) = -coord(dof, 0) / dir(dof);
						lambdaorder(dof,1) = dof;
					}
					break;
					case 2:
					{
						lambdaorder(dof,0) = (PeriodLength - coord(dof,0)) / dir(dof);
						lambdaorder(dof,1) = dof;
					}
					break;
					default:
					{
						lambdaorder(dof,1) = dof;
					}
				}
			}
			// sort the lambdas (ascending values) and indices accordingly
			// in case of multiple shifts
			if(numshifts>1)
			{
				for(int j=0; j<(int)lambdaorder.M()-1; j++)
					for(int k=j+1; k<(int)lambdaorder.M(); k++)
						if(lambdaorder(k,0)<lambdaorder(j,0))
						{
							double temp = lambdaorder(j,0);
							int tempindex = (int)lambdaorder(j,1);
							lambdaorder(j,0) = lambdaorder(k,0);
							lambdaorder(j,1) = lambdaorder(k,1);
							lambdaorder(k,0) = temp;
							lambdaorder(k,1) = tempindex;
						}
			}
			else	// for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
				for(int n=0; n<(int)lambdaorder.N(); n++)
					lambdaorder(0,n) = lambdaorder(shiftdof,n);

			// calculate segment lambdas
			for(int dof=numshifts-1; dof>0; dof--)
				lambdaorder(dof,0) -= lambdaorder(dof-1,0);

			// the idea is to gradually shift the matrix "unshift" back into the volume and, while doing so,
			// calculate the segments except for the last one
			// determine closest boundary component-wise
			for(int shift=0; shift<numshifts; shift++)
			{
				//second point
				for(int j=0 ;j<unshift.M(); j++)
					unshift(j,1) = unshift(j,0) + lambdaorder(shift,0)*dir(j);

				//Calculate Center Point and edgelength of AABB
				LINALG::Matrix<3,1> midpoint;
				LINALG::Matrix<3,1> edgelength;
				for(int i=0; i<(int)midpoint.M(); i++)
				{
					midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
					edgelength(i) = unshift(i,1) - unshift(i,0);
					//Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
					if (edgelength(i)<(*diameter_)[0])
						edgelength(i) = (*diameter_)[0];
				}
				// Calculate limits of AABB of the current segment (which definitely lies in the volume) with extrusion around midpoint
				for(int i=0; i<6; i++)
					if(i%2==0)
						(*bboxlimits)(shift*6+i,0) = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
					else if(i%2==1)
						(*bboxlimits)(shift*6+i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;

				int currshift = (int)lambdaorder(shift,1);
				// shift the coordinates of the second point
				if(cut(currshift)==1.0)
					unshift(currshift,1) += PeriodLength;
				else if(cut(currshift)==2.0)
					unshift(currshift,1) -= PeriodLength;
				// make second point the first and calculate new second point in the next iteration!
				for(int j=0; j<unshift.M(); j++)
					unshift(j,0) = unshift(j,1);
			}

			// the last segment
			for(int dof=0; dof<unshift.M(); dof++)
				unshift(dof,1) = coord(dof,1);

			//Calculate Center Point and edgelength of AABB
			LINALG::Matrix<3,1> midpoint;
			LINALG::Matrix<3,1> edgelength;
			for(int i=0; i<(int)midpoint.M(); i++)
			{
				midpoint(i) = 0.5*(unshift(i,0) + unshift(i,1));
				edgelength(i) = unshift(i,1) - unshift(i,0);
				//Check for edgelength of AABB if too small (bbox parallel to one of the spatial axes)
				if (edgelength(i)<(*diameter_)[0])
					edgelength(i) = (*diameter_)[0];
			}
			// limits of the last bounding box
			for(int i=0; i<6; i++)
				if(i%2==0)
					(*bboxlimits)(numshifts*6+i,0) = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
				else if(i%2==1)
					(*bboxlimits)(numshifts*6+i,0) = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
		}
	}
	//Bring coordinates in case of periodic boundary condition in right order ( "-1e9 signals the bogus value from above)
	//[xmin xmax ymin ymax zmin zmax ...]
	if (periodicBC_)
	{
		double minimum =0.0;    double maximum =0.0;

		for(int l=0; l<(bboxlimits->M())/2;l++)
		{
			// leave loop at first bogus entry
			if((2*l)%6==0 && (*bboxlimits)(2*l,0)==-1e9)
				break;
			minimum = min((*bboxlimits)(2*l,0),(*bboxlimits)(2*l+1,0));
			maximum = max((*bboxlimits)(2*l,0),(*bboxlimits)(2*l+1,0));
			//cout << minimum << endl;
			(*bboxlimits)(2*l,0) = minimum;    (*bboxlimits)(2*l+1,0) = maximum;
		}
	}

	return bboxlimits;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box without a beam (private)      mueller 11/11|
 *----------------------------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> Beam3ContactOctTree::CreateCOBBNoElement(const Epetra_SerialDenseMatrix& coord)
{
	RCP<Epetra_SerialDenseMatrix> boxlimits=Teuchos::null;

	return boxlimits;
}

/*-----------------------------------------------------------------------------------------*
 |  locateAll function (private); Recursive division of a 3-dimensional set.    meier 02/11|
 |  [Oct_ID, element-ID] = locateAll( &allbboxes_)                                          |
 |  Performs recursive tree-like division of a set of Bounding Boxes.                      |
 |  N0 is maximum permissible number of "counted" boxes in the leaf octant.                |
 |  Returns vector IND of the same size as rows of allbboxes_ showing which region each     |
 |  box of a set belongs to; binary matrices BX, BY, BZ where each row shows               |
 |  "binary address" of each region are written to .dat- files each.                       |
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateAll(LINALG::Matrix<1,6>& rootoctantlim, std::vector<std::vector<int> >& bboxesinoctants)
{
#ifdef MEASURETIME
  double t_octree = Teuchos::Time::wallTime();
#endif
  //cout << "\n\tTest locateAll" << endl;
  
  // Parameters and Initialization....................
  // Convert Epetra_MultiVector allbboxes_ to vector(vector(vector)....................
  std::vector<std::vector<double> > allbboxesstdvec(allbboxes_->MyLength(), std::vector<double>(allbboxes_->NumVectors(),0.0));
  for(int i=0; i < allbboxes_->MyLength(); i++)
  {
  	// swapped first and last position of the vectors (the GID)
    allbboxesstdvec[i][0] = (*allbboxes_)[allbboxes_->NumVectors()-1][i];
    for(int j=0; j<allbboxes_->NumVectors()-1; j++)
      allbboxesstdvec[i][j+1] = (*allbboxes_)[j][i];
  }
  //initial tree depth value (will be incremented with each recursive call of locateBox()
  int treedepth = 0;
  //Initialize Vector of Vectors containing limits of Octants for Visualization
  std::vector< std::vector <double> > OctreeLimits;

  // Recursively construct octree; Proc 0 only (parallel computing impossible)
  if(searchdis_.Comm().MyPID()==0)
    locateBox(allbboxesstdvec, rootoctantlim, &OctreeLimits, bboxesinoctants, treedepth);
  else // communicate bbox2octant_
  	bbox2octant_->PutScalar(0.0);

  Epetra_Export exporter(*(searchdis_.ElementColMap()), *(searchdis_.ElementRowMap()));
  Epetra_Import importer(*(searchdis_.ElementColMap()), *(searchdis_.ElementRowMap()));

  Epetra_MultiVector bbox2octantrow(*(searchdis_.ElementRowMap()),4,true);
  bbox2octantrow.Export(*bbox2octant_, exporter, Add);
  bbox2octant_->Import(bbox2octantrow, importer, Insert);

  /*// Write allbboxesstdvec to .dat-file allBBoxesstdvec.dat
  std::ostringstream filename;
  filename << "allbboxesstdvec.dat";
  FILE* fp = NULL;
  //open file to write output data into
  fp = fopen(filename.str().c_str(), "w");
  // write output to temporary stringstream;
  std::stringstream myfile;
  for (int u = 0; u < (int)allbboxesstdvec.size(); u++)
      {
        for (int v = 0; v < (int)allboxesstdvec[0].size(); v++)
        {
          myfile <<scientific<<allbboxesstdvec[u][v] <<" ";
        }
        myfile <<endl;
      }
  //write content into file and close it
  fprintf(fp, myfile.str().c_str());
  fclose(fp);*/


  /*/ For Octree Visualization: Write OctreeLimits to.dat-file OctreeLimits...................
  std::ostringstream filename3;
  filename3 << "OctreeLimits.dat";
  FILE* fp3 = NULL;
  //open file to write output data into
  fp3 = fopen(filename3.str().c_str(), "w");
  // write output to temporary stringstream;
  std::stringstream myfile3;
  for (int u = 0; u < (int)OctreeLimits.size(); u++)
      {
        for (int v = 0; v < (int)OctreeLimits[u].size(); v++)
        {
          myfile3 <<scientific<<OctreeLimits[u][v] <<" ";
        }
        myfile3 <<endl;
      }
   //write content into file and close it
   fprintf(fp3, myfile3.str().c_str());
   fclose(fp3);*/
#ifdef MEASURETIME
   if(!searchdis_.Comm().MyPID())
     cout << "\nOctree building time:\t\t" << Teuchos::Time::wallTime() - t_octree<< " seconds" << endl;
#endif
  return;
}// end of method locateAll


/*-----------------------------------------------------------------------------------------*
 |  locateBox function (private);                                               meier 02/11|
 |  [ind, bx, by, bz] = locateBox(&allbboxes, lim, n0)                                     |
 |  Primitive for locateAll                                                                |
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateBox(std::vector<std::vector<double> > allbboxesstdvec,
                                    LINALG::Matrix<1,6> lim,
                                    std::vector<std::vector<double> >* OctreeLimits,
                                    std::vector<std::vector<int> >& bboxesinoctants,
                                    int& treedepth)
{
  //cout << "\n\n\nTest locateBox...................." << endl;

  // Initialization
  // Number of Bounding Boxes in leaf octant
  int n0 = 10;

  // Zeilenanzahl der Eingangsmatrix
  //cout << "INPUT Number of Boxes:   " << allbboxesstdvec.size() << endl;

  // Check for further recursion by checking current treedepth (second criterion)....................
  //cout << "INPUT Treedepth:  " << treedepth <<"\n";

  // Divide further....................
  // Center of octant....................
  double xcenter = 0.0; double ycenter = 0.0; double zcenter = 0.0;
  xcenter = (lim(0)+lim(1))/2.0;
  ycenter = (lim(2)+lim(3))/2.0;
  zcenter = (lim(4)+lim(5))/2.0;

  double edgelength = fabs(lim(1)-(lim)(0));
  std::vector<LINALG::Matrix<1,6> > limits;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        LINALG::Matrix<1,6> sublim;
        sublim(0) = xcenter + (i-1)*edgelength/2.0;
        sublim(1) = xcenter +  i   *edgelength/2.0;
        sublim(2) = ycenter + (j-1)*edgelength/2.0;
        sublim(3) = ycenter +  j   *edgelength/2.0;
        sublim(4) = zcenter + (k-1)*edgelength/2.0;
        sublim(5) = zcenter +  k   *edgelength/2.0;

        limits.push_back(sublim);
      }
  
  /* Decision to which child box belongs....................
  *
  *           5 ======================== 7
  *           //|                       /||
  *          // |                      //||
  *         //  |                     // ||
  *        //   |                    //  ||
  *       //    |                   //   ||
  *      //     |                  //    ||
  *     //      |                 //     ||
  *    1 ========================= 3     ||
  *    ||       |                ||      ||
  *    ||       |                ||      ||
  *    ||       |      o (center)||      ||
  *    ||      4 ----------------||------ 6
  *    ||      /                 ||     //
  *    ||     /                  ||    //
  *    ||    /                   ||   //
  *    ||   /                    ||  //
  *    ||  /                     || //      y  z
  *    || /                      ||//       | /
  *    ||/                       ||/        |/
  *    0 ========================= 2        ---> x
  *
  */

  //Goes through all suboctants
  for( int oct=0; oct<8; oct++)
  {
    // Define temporary vector of same size as current allbboxesstdvec
    std::vector<std::vector<double> > bboxsubset;
    bboxsubset.clear();

    // Goes through all lines of input allbboxesstdvec....................
    for( int i=0; i<(int)allbboxesstdvec.size(); i++)
    {
      if(periodicBC_)
      {
      	// a bounding box is at maximum divided into 4 subsegments due to periodic boundary conditions
      	for(int isub=0; isub<4; isub++)
      	{
      		// 1)check the first value of the next subsegment. If not used, break ("1+..." because of the bbox gid being at position 0)
    			// 2)loop over the limits of the current octant and check if the current bounding box lies within this octant.
  				// 3)Then, check componentwise and leave after first "hit"
      		if(allbboxesstdvec[i][1+6*isub]!=-1e9)
      		{
						// check for intersection. if yes, there's no need to check further segments ->break
            if(!((limits[oct](0) >= allbboxesstdvec[i][6*isub+2]) || (allbboxesstdvec[i][6*isub+1] >= limits[oct](1)) || (limits[oct](2) >= allbboxesstdvec[i][6*isub+4]) || (allbboxesstdvec[i][6*isub+3] >= limits[oct](3)) || (limits[oct](4) >= allbboxesstdvec[i][6*isub+6]) || (allbboxesstdvec[i][6*isub+5] >= limits[oct](5))))
            {
            	bboxsubset.push_back(allbboxesstdvec[i]);
            	break;
            }
      		}
      		else
      			break;
      	}
      }
      else // standard procedure without periodic boundary conditions
      {
        // Processes colums indices 1 to 6
  			// 2)loop over the limits of the current octant and check if the current bounding box lies within this octant.
  			// 3)Then, check componentwise and leave after first "hit"
        if(!((limits[oct](0) >= allbboxesstdvec[i][2]) || (allbboxesstdvec[i][1] >= limits[oct](1)) || (limits[oct](2) >= allbboxesstdvec[i][4]) || (allbboxesstdvec[i][3] >= limits[oct](3)) || (limits[oct](4) >= allbboxesstdvec[i][6]) || (allbboxesstdvec[i][5] >= limits[oct](5))))
        	bboxsubset.push_back(allbboxesstdvec[i]);
      }

    }// end of for-loop which goes through all elements of input

    // current tree depth
    int currtreedepth = treedepth+1;
    // Check for further recursion by checking number of boxes in octant (first criterion)....................
    int N = (int)bboxsubset.size();
    //cout << "Number of Bounding Boxes in suboctant "<<oct<<":  "<< bboxsubset.size() <<"\t";

    //If to divide further, let LocateBox call itself with updated inputs
    if ((N > n0) && (currtreedepth < maxtreedepth_-1))
      locateBox(bboxsubset, limits[oct], OctreeLimits, bboxesinoctants, currtreedepth);
    else
    {
    	// no further discretization of the volume because either the maximal tree depth or the minimal number of bounding
    	// boxes per octant has been reached
      // Store IDs of Boxes of bboxsubset in vector
      std::vector<int> Box_IDs;
      if(bboxsubset.size()!=0)
      {
        //Push back Limits of suboctants to OctreeLimits
        std::vector<double> suboct_limVec(6,0);
        for(int k=0; k<6; k++)
          suboct_limVec[k]=limits[oct](k);
        OctreeLimits->push_back(suboct_limVec);

        for (int m = 0; m < (int)bboxsubset.size(); m++)
        {
        	// note: the Bounding Box ID is the first column entry of the vector bboxsubset
          Box_IDs.push_back((int)bboxsubset[m][0]);
          // assign current octant number to the bounding box
          for(int n=0; n<bbox2octant_->NumVectors(); n++)
          	if((*bbox2octant_)[n][searchdis_.ElementColMap()->LID((int)bboxsubset[m][0])]<-0.9)
          	{
          		(*bbox2octant_)[n][searchdis_.ElementColMap()->LID((int)bboxsubset[m][0])] = (double)bboxesinoctants.size();
          		break; //leave after finding first empty slot
          	}
        }
        // add bounding box IDs of this octant to the global vector
        bboxesinoctants.push_back(Box_IDs);
      }
    }
  }// end of loop which goes through all suboctants
} // end of method locateBox

/*-----------------------------------------------------------------------------------*
 |  Calculate limits of the root octant (public)                        mueller 11/11|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<1,6> Beam3ContactOctTree::GetRootOctant()
{
	LINALG::Matrix<1,6> lim;
	// if peLriodic BCs are applied
	if(periodicBC_)
	{
		Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
		for(int i=0; i<(int)lim.N(); i++)
		{
			if(i%2==0)
				lim(i)= 0.0;
			else
				lim(i) =  statmechparams.get<double>("PeriodLength", 0.0);
		}
	}
	else // standard procedure to find root box limits
	{
		// initialize
		for(int i=0; i<(int)lim.N(); i++)
			if(i%2==0)
				lim(i) = 1e13;
			else
				lim(i) = -1e13;

		// loop over allbboxes_ and determine the extremes
		for(int i=0; i<allbboxes_->MyLength(); i++)
			for(int j=0; j<allbboxes_->NumVectors()-1; j++)
				if(j%2==0 && (*allbboxes_)[j][i]<lim(j)) // minimal values
					lim(j) = (*allbboxes_)[j][i];
				else if(j%2!=0 && (*allbboxes_)[j][i]>lim(j))
					lim(j) = (*allbboxes_)[j][i];
	}
	return lim;
}

/*-----------------------------------------------------------------------------------*
 |  Bounding Box Intersection function (private)                          meier 01/11|
 |  Intersects Bounding Boxes in same line of map OctreeMap                          |
 |  Gives back vector of intersection pairs                                          |
 *----------------------------------------------------------------------------------*/
void Beam3ContactOctTree::BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
																					 std::vector<std::vector<int> >& bboxesinoctants,
                                           vector<RCP<Beam3contact> >* contactpairs)
{
#ifdef MEASURETIME
  double t_search = Teuchos::Time::wallTime();
#endif
  //cout << "\n\nIntersectionAABB Test..................." << endl;

  //determine maximum depth of OctreeMap
  int maxdepthlocal = 0;
  int bboxlengthlocal = 0;
  if(discret_.Comm().MyPID()==0)
  {
    bboxlengthlocal = (int)bboxesinoctants.size();
    for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
      if((int)bboxesinoctants[i].size()>maxdepthlocal)
        maxdepthlocal = (int)bboxesinoctants[i].size();
  }

  int maxdepthglobal = 0;
  int bboxlengthglobal = 0;
  discret_.Comm().MaxAll(&maxdepthlocal, &maxdepthglobal, 1);
  discret_.Comm().MaxAll(&bboxlengthlocal, &bboxlengthglobal, 1);

  // build Epetra_MultiVector from OtreeMap
  // build temporary, fully overlapping map and row map
  // create crosslinker maps
  std::vector<int> gids;
  for (int i=0 ; i<bboxlengthglobal; i++ )
    gids.push_back(i);
  // crosslinker column and row map
  Epetra_Map octtreerowmap((int)gids.size(), 0, discret_.Comm());
  Epetra_Map octtreemap(-1, (int)gids.size(), &gids[0], 0, discret_.Comm());
  // build Epetra_MultiVectors which hold the BBs of the OctreeMap; for communication
  bboxesinoctants_ = rcp(new Epetra_MultiVector(octtreemap,maxdepthglobal+1));
  Epetra_MultiVector bboxinoctrow(octtreerowmap,maxdepthglobal+1, true);
  // fill bboxinoct for Proc 0
  if(searchdis_.Comm().MyPID()==0)
  {
  	bboxesinoctants_->PutScalar(-9.0);
    for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
      for(int j=0; j<(int)bboxesinoctants[i].size(); j++)
        (*bboxesinoctants_)[j][i] = bboxesinoctants[i][j];
  }
  else
  	bboxesinoctants_->PutScalar(0.0);

  // Communication
  Epetra_Export exporter(octtreemap, octtreerowmap);
  Epetra_Import importer(octtreemap, octtreerowmap);
  bboxinoctrow.Export(*bboxesinoctants_,exporter,Add);
  bboxesinoctants_->Import(bboxinoctrow,importer,Insert);

  //Algorithm begins

  // Build contact pair Map
  std::map<int, std::vector<int> > contactpairmap;
  // create contact pair vector, redundant on all Procs; including redundant pairs
  //for-loop lines of map
  for (int i=0 ; i<bboxesinoctants_->MyLength(); i++ )
  {
    //for-loop index first box
    for(int j=0; j<bboxesinoctants_->NumVectors(); j++)
    {
      // first box ID
      std::vector<int> bboxIDs(2,0);
      bboxIDs[0] = (int)(*bboxesinoctants_)[j][i];

      //for-loop second box
      for(int k=j+1; k<bboxesinoctants_->NumVectors(); k++)
      {
        bboxIDs[1] = (int)(*bboxesinoctants_)[k][i];

        // exclude element pairs sharing one node
        // contact flag
        bool considerpair = false;
        // only consider existing bounding boxes, i.e. no dummy entries "-9.0"
        if(bboxIDs[0]>-1 && bboxIDs[1]>-1)
        {
        	considerpair = true;
          DRT::Element* element1 = searchdis_.gElement(bboxIDs[0]);
          DRT::Element* element2 = searchdis_.gElement(bboxIDs[1]);
					for(int k=0; k<element1->NumNode(); k++)
					{
						for(int l=0; l<element2->NumNode(); l++)
							if(element1->NodeIds()[k]==element2->NodeIds()[l])
							{
								considerpair = false;
								break;
							}

						// break after first found match
						if(!considerpair)
							break;
					}
        }

        if (considerpair)
        {
        	// apply different bounding box intersection schemes
        	bool intersection = false;
        	switch(boundingbox_)
        	{
        		case Beam3ContactOctTree::axisaligned:
        			intersection = IntersectionAABB(bboxIDs);
						break;
        		case Beam3ContactOctTree::cyloriented:
        			intersection = IntersectionCOBB(bboxIDs);
						break;
        		default: dserror("No or an invalid Octree type was chosen. Check your input file!");
        	}

          if (intersection)
          {
            // note: creation of unique "first" entries in map, attention: IDs identical to crosslinker GIDs!!
            int mapfirst = (bboxIDs[0] + 1)*basisnodes_ + bboxIDs[1];
            contactpairmap.insert ( pair<int, vector<int> > (mapfirst, bboxIDs));
          }
        }
      }
    }
  }

  // build Pair Vector from contactpairmap
  std::map<int, std::vector<int> >::iterator it;
  int counter = 0;
  for ( it=contactpairmap.begin() ; it !=contactpairmap.end(); it++ )
  {
    counter++;
    //if(!discret_.Comm().MyPID())
      //cout << scientific << (*it).first <<"  "<< ((*it).second)[0]<<" "<< ((*it).second)[1]<<endl;
    int collid1 = searchdis_.ElementColMap()->LID(((*it).second)[0]);
    int collid2 = searchdis_.ElementColMap()->LID(((*it).second)[1]);

    DRT::Element* tempele1 = searchdis_.lColElement(collid1);
    DRT::Element* tempele2 = searchdis_.lColElement(collid2);

    // matrices to store nodal coordinates
    Epetra_SerialDenseMatrix ele1pos(3,tempele1->NumNode());
    Epetra_SerialDenseMatrix ele2pos(3,tempele2->NumNode());

    // store nodal coordinates of element 1
    for (int m=0;m<tempele1->NumNode();++m)
    {
     int tempGID = (tempele1->NodeIds())[m];
     LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
     for(int n=0;n<3;n++) ele1pos(n,m) = temppos(n);
    }

    // store nodal coordinates of element 2
    for (int m=0;m<tempele2->NumNode();++m)
    {
     int tempGID = (tempele2->NodeIds())[m];
     LINALG::Matrix<3,1> temppos = currentpositions[tempGID];
     for(int n=0;n<3;n++) ele2pos(n,m) = temppos(n);
    }

    // add to pair vector
    contactpairs->push_back(rcp (new Beam3contact(discret_,searchdis_,dofoffset_,tempele1,tempele2,ele1pos,ele2pos)));
  }
  //if(!discret_.Comm().MyPID())
    //cout<<"number of boxes: "<<counter<<endl;

#ifdef MEASURETIME
  double isectimelocal = Teuchos::Time::wallTime() - t_search;
  double isectimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&isectimelocal, &isectimeglobal, 1);
  discret_.Comm().Barrier();
  if(!searchdis_.Comm().MyPID())
    cout << "Intersection time:\t\t" << isectimeglobal << " seconds\n\n";
#endif

  /*/Print ContactPairs to .dat-file and plot with Matlab....................
  std::ostringstream filename2;
  filename2 << "ContactPairs.dat";
  FILE* fp2 = NULL;
  fp2 = fopen(filename2.str().c_str(), "w");
  fclose(fp2);
  //open file to write output data into
  // write output to temporary stringstream;
  std::stringstream myfile2;
  fp2 = fopen(filename2.str().c_str(), "a");
  for (int i=0;i<(int)contactpairs->size();i++)
    myfile2 << ((*contactpairs)[i]->Element1())->Id() <<"  "<< ((*contactpairs)[i]->Element2())->Id() <<endl;
  fprintf(fp2, myfile2.str().c_str());
  fclose(fp2);*/
  return;
}//end of method BoundingBoxIntersection()

/*-----------------------------------------------------------------------------------*
 |  Axis Aligned Bounding Box Intersection function when both bounding boxes         |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionAABB(std::vector<int>& bboxIDs)
{
	bool intersection = false;
	// translate box / element GIDs to ElementColMap()-LIDs
	// note: GID and ColumnMap LID are usually the same except for crosslinker elements from statistical mechanics
  int entry1 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int entry2 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  //Initialization....................
  double a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax;
  double b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax;

	if(periodicBC_)
	{
		// note: n shifts means n+1 segments
		for(int i=0; i<(int)(*numshifts_)[entry1]+1; i++)
		{
			for(int j=0; j<(int)(*numshifts_)[entry2]+1; j++)
			{
			  //Intersection Test
			  a_xmin=(*allbboxes_)[i*6][entry1];   	 a_xmax=(*allbboxes_)[i*6+1][entry1];
			  a_ymin=(*allbboxes_)[i*6+2][entry1];   a_ymax=(*allbboxes_)[i*6+3][entry1];
			  a_zmin=(*allbboxes_)[i*6+4][entry1];   a_zmax=(*allbboxes_)[i*6+5][entry1];

			  b_xmin=(*allbboxes_)[j*6][entry2];   	 b_xmax=(*allbboxes_)[j*6+1][entry2];
			  b_ymin=(*allbboxes_)[j*6+2][entry2];   b_ymax=(*allbboxes_)[j*6+3][entry2];
			  b_zmin=(*allbboxes_)[j*6+4][entry2];   b_zmax=(*allbboxes_)[j*6+5][entry2];

			  // if intersection exists, return true
			  if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
			  {
			  	intersection = true;
			  	break;
			  }
			}
			if(intersection)
				break;
		}
	}
	else	// standard procedure without periodic boundary conditions
	{
	  //Intersection Test
	  a_xmin=(*allbboxes_)[0][entry1];   a_xmax=(*allbboxes_)[1][entry1];
	  a_ymin=(*allbboxes_)[2][entry1];   a_ymax=(*allbboxes_)[3][entry1];
	  a_zmin=(*allbboxes_)[4][entry1];   a_zmax=(*allbboxes_)[5][entry1];

	  b_xmin=(*allbboxes_)[0][entry2];   b_xmax=(*allbboxes_)[1][entry2];
	  b_ymin=(*allbboxes_)[2][entry2];   b_ymax=(*allbboxes_)[3][entry2];
	  b_zmin=(*allbboxes_)[4][entry2];   b_zmax=(*allbboxes_)[5][entry2];

	  // if intersection exists, return true
	  if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
	  	intersection = true;
	}


  return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Cylindrical Oriented Bounding Box Intersection function when both bounding boxes |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionCOBB(std::vector<int>& bboxIDs)
{
	bool intersection = false;
	return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Axis Aligned Bounding Box Intersection function when only one of the             |
 |  bounding boxes represents an actual finite element (private)        mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IsecAABBNoElement(const int& bboxinoct, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  bool intersection = false;
	// translate box / element GIDs to ElementColMap()-LIDs
	// note: GID and ColumnMap LID are usually the same except for crosslinker elements from statistical mechanics
	int bboxcolid = searchdis_.ElementColMap()->LID(bboxinoct);

	// number of shifts in bboxlimits
	int numsegments = bboxlimits->M()/6;

	//Initialization....................
	double a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax;
	double b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax;

	if(periodicBC_)
	{
		// note: n shifts means n+1 segments
		for(int i=0; i<(int)(*numshifts_)[bboxcolid]+1; i++)
		{
			for(int j=0; j<numsegments; j++)
			{
				//Intersection Test
				a_xmin=(*allbboxes_)[i*6][bboxcolid];   	a_xmax=(*allbboxes_)[i*6+1][bboxcolid];
				a_ymin=(*allbboxes_)[i*6+2][bboxcolid];   a_ymax=(*allbboxes_)[i*6+3][bboxcolid];
				a_zmin=(*allbboxes_)[i*6+4][bboxcolid];   a_zmax=(*allbboxes_)[i*6+5][bboxcolid];

				b_xmin=(*bboxlimits)(j*6,0);   	 b_xmax=(*bboxlimits)(j*6+1,0);
				b_ymin=(*bboxlimits)(j*6+2,0);   b_ymax=(*bboxlimits)(j*6+3,0);
				b_zmin=(*bboxlimits)(j*6+4,0);   b_zmax=(*bboxlimits)(j*6+5,0);

				// if intersection exists, return true
				if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
				{
					intersection = true;
					break;
				}
			}
			if(intersection)
				break;
		}
	}
	else	// standard procedure without periodic boundary conditions
	{
		//Intersection Test
		a_xmin=(*allbboxes_)[0][bboxcolid];   a_xmax=(*allbboxes_)[1][bboxcolid];
		a_ymin=(*allbboxes_)[2][bboxcolid];   a_ymax=(*allbboxes_)[3][bboxcolid];
		a_zmin=(*allbboxes_)[4][bboxcolid];   a_zmax=(*allbboxes_)[5][bboxcolid];

		b_xmin=(*bboxlimits)(0,0);   b_xmax=(*bboxlimits)(1,0);
		b_ymin=(*bboxlimits)(2,0);   b_ymax=(*bboxlimits)(3,0);
		b_zmin=(*bboxlimits)(4,0);   b_zmax=(*bboxlimits)(5,0);

		// if intersection exists, return true
		if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
			intersection = true;
	}

	return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Cylindrical Oriented Bounding Box Intersection function when only one of the     |
 |  bounding boxes represents an actual finite element (private)        mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IsecCOBBNoElement(const int& bboxinoct, RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
	bool intersection = false;
	return intersection;
}
#endif /*CCADISCRET*/
