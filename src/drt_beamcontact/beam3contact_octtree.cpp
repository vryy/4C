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
  {
  	Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
		if(statmechparams.get("PeriodLength",0.0)>0.0)
			periodicBC_ = true;
		else
			periodicBC_ = false;
  }
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
}


/*----------------------------------------------------------------------*
 |  calls the almighty Octtree (public)                      meier 01/11|
 *----------------------------------------------------------------------*/
vector<RCP<Beam3contact> > Beam3ContactOctTree::OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions)
{
  // initialize beam diameter
  diameter_ = rcp(new Epetra_Vector(*(searchdis_.ElementColMap())));

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
  // (components 6,...,11 contain bounding box limits in case of periodic boundary conditions (2nd part of the box))
  // (component 12 containts element ID)
  // in case of periodic boundary conditions, we need 12+1 vectors within allBBoxes, without periodic BCs 6+1
  RCP<Epetra_MultiVector> allBBoxes = Teuchos::null;
  if(periodicBC_)
  	allBBoxes = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),13, true));
  else
  	allBBoxes = rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),7, true));

  // build axis aligned bounding boxes
  CreateBoundingBoxes(currentpositions, allBBoxes);

  // call recursive octtree build
  // vector for assigning bounding boxes to octants
  std::vector<std::vector<int> > bboxesinoctants;
  locateAll(allBBoxes, bboxesinoctants);

  // intersection checks
  vector<RCP<Beam3contact> > contactpairs;
  BoundingBoxIntersection(currentpositions, &contactpairs, &bboxesinoctants, allBBoxes);

  return contactpairs;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)                 meier 01/11|
 |  generates bounding boxes extended with factor 1.05                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, RCP<Epetra_MultiVector> allBBoxes)
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
    			CreateAABB(coord, allBBoxes, elecolid, PeriodLength);
				break;
    		case Beam3ContactOctTree::cyloriented:
    			CreateCOBB(coord, allBBoxes, elecolid, PeriodLength);
				break;
    		default: dserror("No or an invalid Octree type was chosen. Check your input file!");
    	}
		}
	} //end for-loop which goes through all elements

  // communication of findings
  Epetra_MultiVector allBBoxesrow(*(searchdis_.ElementRowMap()),allBBoxes->NumVectors(),true);
  Epetra_Export exporter(*(searchdis_.ElementColMap()),*(searchdis_.ElementRowMap()));
  Epetra_Import importer(*(searchdis_.ElementColMap()),*(searchdis_.ElementRowMap()));
  allBBoxesrow.Export(*allBBoxes, exporter, Add);
  allBBoxes->Import(allBBoxesrow,importer,Insert);

  //Test: print allBBoxes->...................
  //cout << "\n\tTest extendedAABB" << endl;

  /*/Visualization print allBBoxes to .dat-file and plot with Matlab....................
  std::ostringstream filename;
  filename << "extendedAABB.dat";
  FILE* fp = NULL;
  //open file to write output data into
  fp = fopen(filename.str().c_str(), "w");
  // write output to temporary stringstream;
  std::stringstream myfile;
  for (int u = 0; u < allBBoxes->MyLength(); u++)
  {
    for (int v = 0; v < allBBoxes->NumVectors(); v++)
    {
      myfile <<scientific<<setprecision(10)<<(*allBBoxes)[v][u] <<" ";
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
void Beam3ContactOctTree::CreateAABB(const Epetra_SerialDenseMatrix& coord, RCP<Epetra_MultiVector> allBBoxes, const int& elecolid, const double& PeriodLength)
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
  LINALG::SerialDenseMatrix cut(3, 1, true);

  /* "coord" currently holds the shifted set of coordinates.
   * In order to determine the correct vector "dir" of the visualization at the boundaries,
   * a copy of "coord" with adjustments in the proper places is introduced*/
  LINALG::SerialDenseMatrix unshift = coord;

  // Compute "cut"-matrix (only in case of periodic BCs)
  if(periodicBC_)
  {
		for (int dof=0; dof<ndim; dof++)
		{
			if (fabs(coord(dof,1)-PeriodLength-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
			{
				cut(dof,0) = 1;
				unshift(dof,1) -= PeriodLength;
			}
			if (fabs(coord(dof,1)+PeriodLength - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
			{
				cut(dof,0) = 2;
				unshift(dof,1) += PeriodLength;
			}
		}
  }

  // handle broken elements ( without periodic BCs, the "if" remains untapped since sum(cut) will always be zero
  if (cut(0,0) + cut(1,0) + cut(2,0) > 0)
  {
    //compute direction vector between first(i-th) and second(i+1-th) node of element (normed):
    LINALG::Matrix<3, 1> dir;
    double ldir = 0.0;
    for (int dof = 0; dof < ndim; dof++)
    {
      dir(dof) = unshift(dof,1) - unshift(dof,0);
      ldir += dir(dof) * dir(dof);
    }
    for (int dof = 0; dof < ndim; dof++)
      dir(dof) /= ldir;

    //from node 0 to nearest boundary where element is broken you get by vector X + lambda0*dir
    double lambda0 = 1e4;
    for (int dof = 0; dof < ndim; dof++)
    {
      if (cut(dof, 0) == 1)
      {
        if (fabs(-coord(dof, 0) / dir(dof)) < fabs(lambda0))
          lambda0 = -coord(dof, 0) / dir(dof);
      }
      else if (cut(dof, 0) == 2)
      {
        if (fabs((PeriodLength - coord(dof, 0)) / dir(dof)) < fabs(lambda0))
          lambda0 = (PeriodLength - coord(dof, 0)) / dir(dof);
      }
    }

    //from node 1 to nearest boundary where element is broken you get by vector X + lambda1*dir
    double lambda1 = 1e4;
    for (int dof = 0; dof < ndim; dof++)
    {
      if (cut(dof, 0) == 2)
      {
        if (fabs(-coord(dof,1) / dir(dof)) < fabs(lambda1))
          lambda1 = -coord(dof,1) / dir(dof);
      }
      else if (cut(dof, 0) == 1)
      {
        if (fabs((PeriodLength - coord(dof,1)) / dir(dof)) < fabs(lambda1))
          lambda1 = (PeriodLength - coord(dof,1)) / dir(dof);
      }
    }

    //define output coordinates for broken elements, first segment
    LINALG::SerialDenseMatrix coordout=coord;
    for(int i=0 ;i<coordout.M(); i++)
      coordout(i,1) = coord(i,0) + lambda0*dir(i);

    //do process with nodecoords0 and nodecoords_sp1
    //Calculate Center Point of AABB
    LINALG::Matrix<3,1> midpoint;
    for(int i=0; i<(int)midpoint.M(); i++)
      midpoint(i) = 0.5*(coordout(i,0) + coordout(i,1));

    //Calculate edgelength of AABB
    LINALG::Matrix<3,1> edgelength;
    for(int i=0; i<(int)edgelength.M(); i++)
      edgelength(i) = coordout(i,1) - coordout(i,0);

    //Check for edgelength of AABB
    for(int i=0; i<(int)edgelength.M(); i++)
      if (edgelength(i)<(*diameter_)[elecolid])
        edgelength(i) = (*diameter_)[elecolid];

    // Calculate limits of AABB with extrusion around midpoint
    for(int i=0; i<6; i++)
      if(i%2==0)
        (*allBBoxes)[i][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
      else if(i%2!=0)
        (*allBBoxes)[i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;


    //define output coordinates for broken elements, second segment
    coordout = coord;
    for(int j=0; j<coordout.M(); j++)
      coordout(j,0) = coord(j,1)+lambda1*dir(j);

    //do process again with nodecoords_sp2 and nodecoords1
    //Calculate Center Point of AABB
    for(int i=0; i<(int)midpoint.M(); i++)
      midpoint(i) = 0.5*(coordout(i,0) + coordout(i,1));

    //Calculate edgelength of AABB
    for(int i=0; i<(int)edgelength.M(); i++)
      edgelength(i) = coordout(i,1) - coordout(i,0);

    //Check for edgelength of AABB
		for(int i=0; i<(int)edgelength.M(); i++)
			if (edgelength(i)<(*diameter_)[elecolid])
				edgelength(i) = (*diameter_)[elecolid];

    // Calculate limits of AABB with extrusion around midpoint
    for(int i=0; i<6; i++)
      if(i%2==0)
        (*allBBoxes)[i+6][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
      else if(i%2==1)
        (*allBBoxes)[i+6][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;

    // store GID (=box number)
    (*allBBoxes)[allBBoxes->NumVectors()-1][elecolid] = elegid;
  } // end of broken AABB
  else
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
        (*allBBoxes)[i][elecolid] = midpoint((int)floor((double)i/2.0)) - 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
      else if(i%2==1)
        (*allBBoxes)[i][elecolid] = midpoint((int)floor((double)i/2.0)) + 0.5*edgelength((int)floor((double)i/2.0))*extrusionfactor;
    // fill the latter 6 entries with bogus values (in case of periodic BCs)
    if(periodicBC_)
			for(int i=0; i<6; i++)
				(*allBBoxes)[i+6][elecolid] = -1e9;


    (*allBBoxes)[allBBoxes->NumVectors()-1][elecolid] = elegid;
  }

  //18.02.2011 Bring coordinates in case of periodic boundary condition in right order ( "-1.9 signals the bogus value from above)
  //[xmin xmax ymin ymax zmin zmax]
  if (periodicBC_ && (*allBBoxes)[6][elecolid]!=-1e9)
  {
    double minimum =0.0;    double maximum =0.0;

    for(int l=0; l<6;l++)
    {
			minimum = min((*allBBoxes)[2*l][elecolid],(*allBBoxes)[2*l+1][elecolid]);
			maximum = max((*allBBoxes)[2*l][elecolid],(*allBBoxes)[2*l+1][elecolid]);
			//cout << minimum << endl;
			(*allBBoxes)[2*l][elecolid] = minimum;    (*allBBoxes)[2*l+1][elecolid] = maximum;
    }// end of correct
  } // end of normal AAABB

	return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box   (private)                   mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateCOBB(const Epetra_SerialDenseMatrix& coord, RCP<Epetra_MultiVector> allBBoxes, const int& elecolid, const double& PeriodLength)
{
	return;
}

/*-----------------------------------------------------------------------------------------*
 |  locateAll function (private); Recursive division of a 3-dimensional set.    meier 02/11|
 |  [Oct_ID, element-ID] = locateAll( &allBBoxes)                                          |
 |  Performs recursive tree-like division of a set of Bounding Boxes.                      |
 |  N0 is maximum permissible number of "counted" boxes in the leaf octant.                |
 |  Returns vector IND of the same size as rows of allBBoxes showing which region each     |
 |  box of a set belongs to; binary matrices BX, BY, BZ where each row shows               |
 |  "binary address" of each region are written to .dat- files each.                       |
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateAll(RCP<Epetra_MultiVector> allBBoxes, std::vector<std::vector<int> >& bboxesinoctants)
{
#ifdef MEASURETIME
  double t_octree = Teuchos::Time::wallTime();
#endif
  //cout << "\n\tTest locateAll" << endl;
  
  // Parameters and Initialization....................
  // Limits of control volumina [0 PeriodLength 0 PeriodLength 0 PeriodLength]....................
  LINALG::Matrix<1,6> lim;
  // if periodic BCs are applied
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

  	// loop over allBBoxes and determine the extremes
  	for(int i=0; i<allBBoxes->MyLength(); i++)
  		for(int j=0; j<allBBoxes->NumVectors()-1; j++)
  			if(j%2==0 && (*allBBoxes)[j][i]<lim(j)) // minimal values
  				lim(j) = (*allBBoxes)[j][i];
  			else if(j%2!=0 && (*allBBoxes)[j][i]>lim(j))
  				lim(j) = (*allBBoxes)[j][i];
  }

  // Convert Epetra_MultiVector allBBoxes to vector(vector(vector)....................
  std::vector<std::vector<double> > allBBoxesstdvec(allBBoxes->MyLength(), std::vector<double>(allBBoxes->NumVectors(),0.0));
  for(int i=0; i < allBBoxes->MyLength(); i++)
  {
  	// swapped first and last position of the vectors (the GID)
    allBBoxesstdvec[i][0] = (*allBBoxes)[allBBoxes->NumVectors()-1][i];
    for(int j=0; j<allBBoxes->NumVectors()-1; j++)
      allBBoxesstdvec[i][j+1] = (*allBBoxes)[j][i];
  }
  //initial tree depth value (will be incremented with each recursive call of locateBox()
  int treedepth = 0;
  //Initialize Vector of Vectors containing limits of Octants for Visualization
  std::vector< std::vector <double> > OctreeLimits;

  // Recursively construct octree; Proc 0 only (parallel computing impossible)
  if(searchdis_.Comm().MyPID()==0)
    locateBox(allBBoxesstdvec, lim, &OctreeLimits, bboxesinoctants, treedepth);

  /*// Write allBBoxesstdvec to .dat-file allBBoxesstdvec.dat
  std::ostringstream filename;
  filename << "allBBoxesstdvec.dat";
  FILE* fp = NULL;
  //open file to write output data into
  fp = fopen(filename.str().c_str(), "w");
  // write output to temporary stringstream;
  std::stringstream myfile;
  for (int u = 0; u < (int)allBBoxesstdvec.size(); u++)
      {
        for (int v = 0; v < (int)allBBoxesstdvec[0].size(); v++)
        {
          myfile <<scientific<<allBBoxesstdvec[u][v] <<" ";
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
 |  [ind, bx, by, bz] = locateBox(&allBBoxes, lim, n0)                                     |
 |  Primitive for locateAll                                                                |
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateBox(std::vector<std::vector<double> > allBBoxesstdvec,
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
  //cout << "INPUT Number of Boxes:   " << allBBoxesstdvec.size() << endl;

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
    // Define temporary vector of same size as current allBBoxesstdvec
    std::vector<std::vector<double> > tmp_vec2;
    tmp_vec2.clear();

    // Goes through all lines of input allBBoxesstdvec....................
    for( int i=0; i<(int)allBBoxesstdvec.size(); i++)
    {
      // Processes colums indices 1 to 6, e.g. first boxes
      if(!((limits[oct](0) >= allBBoxesstdvec[i][2]) || (allBBoxesstdvec[i][1] >= limits[oct](1)) || (limits[oct](2) >= allBBoxesstdvec[i][4]) || (allBBoxesstdvec[i][3] >= limits[oct](3)) || (limits[oct](4) >= allBBoxesstdvec[i][6]) || (allBBoxesstdvec[i][5] >= limits[oct](5))))
        tmp_vec2.push_back(allBBoxesstdvec[i]);// end first boxes

      // Processes colums indices 7 to 12, e.g. second boxes due to periodic boundary condition
      if(periodicBC_ && allBBoxesstdvec[i][7] !=-1e9)
					if(!((limits[oct](0) >= allBBoxesstdvec[i][8]) || (allBBoxesstdvec[i][7] >= limits[oct](1)) || (limits[oct](2) >= allBBoxesstdvec[i][10]) || (allBBoxesstdvec[i][9] >= limits[oct](3)) || (limits[oct](4) >= allBBoxesstdvec[i][12]) || (allBBoxesstdvec[i][11] >= limits[oct](5))))
						tmp_vec2.push_back(allBBoxesstdvec[i]);//end second boxes

    }// end of for-loop which goes through all elements of input

    // current tree depth
    int currtreedepth = treedepth+1;

    // Check for further recursion by checking number of boxes in octant (first criterion)....................
    int N = (int)tmp_vec2.size();
    //cout << "Number of Bounding Boxes in suboctant "<<oct<<":  "<< tmp_vec2.size() <<"\t";

    //If to divide further, let LocateBox call itself with updated inputs
    if ((N > n0) && (currtreedepth < maxtreedepth_-1))
      locateBox(tmp_vec2, limits[oct], OctreeLimits, bboxesinoctants, currtreedepth);
    else
    {
    	// no further discretization of the volume because either the maximal tree depth or the minimal number of bounding
    	// boxes per octant has been reached
      // Store IDs of Boxes of tmp_vec2 in vector
      std::vector<int> Box_IDs;
      if(tmp_vec2.size()!=0)
      {
        //Push back Limits of suboctants to OctreeLimits
        std::vector<double> suboct_limVec(6,0);
        for(int k=0; k<6; k++)
          suboct_limVec[k]=limits[oct](k);
        OctreeLimits->push_back(suboct_limVec);

        for (int m = 0; m < (int)tmp_vec2.size(); m++)
          Box_IDs.push_back((int)tmp_vec2[m][0]);

        bboxesinoctants.push_back(Box_IDs);
      }
    }
  }// end of loop which goes through all suboctants
} // end of method locateBox


/*-----------------------------------------------------------------------------------*
 |  Bounding Box Intersection function (private)                          meier 01/11|
 |  Intersects Bounding Boxes in same line of map OctreeMap                          |
 |  Gives back vector of intersection pairs                                          |
 *----------------------------------------------------------------------------------*/
void Beam3ContactOctTree::BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
                                           vector<RCP<Beam3contact> >* contactpairs,
                                           std::vector<std::vector<int> >* bboxesinoctants,
                                           RCP<Epetra_MultiVector> allBBoxes)
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
    bboxlengthlocal = (int)bboxesinoctants->size();
    for (int i=0 ; i<(int)bboxesinoctants->size(); i++ )
      if((int)(*bboxesinoctants)[i].size()>maxdepthlocal)
        maxdepthlocal = (int)(*bboxesinoctants)[i].size();
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
  // build Epetra_MultiVectors which hold the AABBs of the OctreeMap; for communication
  Epetra_MultiVector BBoxinOct(octtreemap,maxdepthglobal+1);
  Epetra_MultiVector BBoxinOctRow(octtreerowmap,maxdepthglobal+1, true);
  // fill BBoxinOct for Proc 0
  if(searchdis_.Comm().MyPID()==0)
  {
    BBoxinOct.PutScalar(-9.0);
    for (int i=0 ; i<(int)bboxesinoctants->size(); i++ )
      for(int j=0; j<(int)(*bboxesinoctants)[i].size(); j++)
        BBoxinOct[j][i] = (*bboxesinoctants)[i][j];
  }
  else
    BBoxinOct.PutScalar(0.0);

  // Communication
  Epetra_Export exporter(octtreemap, octtreerowmap);
  Epetra_Import importer(octtreemap, octtreerowmap);
  BBoxinOctRow.Export(BBoxinOct,exporter,Add);
  BBoxinOct.Import(BBoxinOctRow,importer,Insert);

  //Algorithm begins

  // Build contact pair Map
  std::map<int, std::vector<int> > contactpairmap;
  // create contact pair vector, redundant on all Procs; including redundant pairs
  //for-loop lines of map
  for (int i=0 ; i<BBoxinOct.MyLength(); i++ )
  {
    //for-loop index first box
    for(int j=0; j<BBoxinOct.NumVectors(); j++)
    {
      // first box ID
      std::vector<int> bboxIDs(2,0);
      bboxIDs[0] = (int)BBoxinOct[j][i];

      //for-loop second box
      for(int k=j+1; k<BBoxinOct.NumVectors(); k++)
      {
        bboxIDs[1] = (int)BBoxinOct[k][i];

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
        			intersection = IntersectionAABB(bboxIDs, allBBoxes);
						break;
        		case Beam3ContactOctTree::cyloriented:
        			intersection = IntersectionCOBB(bboxIDs, allBBoxes);
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
}//end of method IntersectAABB

/*-----------------------------------------------------------------------------------*
 |  Axis Aligned Bounding Box Intersection function (private)           mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionAABB(std::vector<int>& bboxIDs, RCP<Epetra_MultiVector> allBBoxes)
{
  //Initialization....................
  double a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax;
  double b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax;

	bool intersection = false;
	// translate box / element GIDs to ElementColMap()-LIDs
	// note: GID and ColumnMap LID are usually the same except for crosslinker elements from statistical mechanics
  int entry1 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int entry2 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  //Intersection Test

  a_xmin=(*allBBoxes)[0][entry1];   a_xmax=(*allBBoxes)[1][entry1];
  a_ymin=(*allBBoxes)[2][entry1];   a_ymax=(*allBBoxes)[3][entry1];
  a_zmin=(*allBBoxes)[4][entry1];   a_zmax=(*allBBoxes)[5][entry1];

  b_xmin=(*allBBoxes)[0][entry2];   b_xmax=(*allBBoxes)[1][entry2];
  b_ymin=(*allBBoxes)[2][entry2];   b_ymax=(*allBBoxes)[3][entry2];
  b_zmin=(*allBBoxes)[4][entry2];   b_zmax=(*allBBoxes)[5][entry2];

  // if intersection exists, return true
  if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
  	intersection = true;

  return intersection;
}
/*-----------------------------------------------------------------------------------*
 |  Cylindrical Oriented Bounding Box Intersection function (private)   mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionCOBB(std::vector<int>& bboxIDs, RCP<Epetra_MultiVector> allBBoxes)
{
	bool intersection = false;
	return intersection;
}

#endif /*CCADISCRET*/
