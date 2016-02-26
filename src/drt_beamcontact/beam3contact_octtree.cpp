/*!----------------------------------------------------------------------
\file beam3contact_octtree.cpp
\brief Octtree for beam contact search

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*----------------------------------------------------------------------*/

#include "beam3contact_octtree.H"
#include "beam3contact_defines.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
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

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "beam3contact_manager.H"

//#define OCTREEDEBUG

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/11|
 *----------------------------------------------------------------------*/
Beam3ContactOctTree::Beam3ContactOctTree(Teuchos::ParameterList& params, DRT::Discretization& discret,DRT::Discretization& searchdis):
btsph_(false),
btsol_(false),
discret_(discret),
searchdis_(searchdis),
basisnodes_(discret.NumGlobalNodes())
{
  extrusionvalue_ = Teuchos::rcp(new std::vector<double>);
  extrusionvalue_->clear();

  INPAR::BEAMCONTACT::OctreeType bboxtype_input = INPAR::BEAMCONTACT::boct_none;

  // This OctTree may be used for contact search as well as search for potential-based interaction
  // it will thus be called with slightly different parameter sets
  // first find out which params we got and extract octtree specifications

  if (params.name() == "DAT FILE->BEAM CONTACT")
  {
    // octree specs
    bboxtype_input = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(params, "BEAMS_OCTREE");

    // additive or multiplicative extrusion of bounding boxes
    if(DRT::INPUT::IntegralValue<int>(params, "BEAMS_ADDITEXT"))
      additiveextrusion_ = true;
    else
      additiveextrusion_ = false;

    // extrusion factor(s)
    // depending on the kind of extrusion, more than one extrusion value will be necessary
    // AABB: one value for all dimensions
    // COBB: 1. value for axial extrusion, 2. value for radial extrusion
    // SPBB: one value for radial extrusion

    std::istringstream PL(Teuchos::getNumericStringParameter(params,"BEAMS_EXTVAL"));
    std::string word;
    char* input;
    while (PL >> word)
      extrusionvalue_->push_back(std::strtod(word.c_str(), &input));

    // max tree depth
    maxtreedepth_ = params.get<int>("BEAMS_TREEDEPTH", 6);
    // max number of bounding boxes per leaf octant
    minbboxesinoctant_ = params.get<int>("BEAMS_BOXESINOCT",8);

    btsph_ = DRT::INPUT::IntegralValue<int>(params, "BEAMS_BTSPH");
    btsol_ = DRT::INPUT::IntegralValue<int>(params, "BEAMS_BTSOL");
  }
  else if (params.name() == "DAT FILE->BEAM POTENTIAL")
  {
    bboxtype_input = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OctreeType>(params, "BEAMPOT_OCTREE");

    additiveextrusion_ = true;
    extrusionvalue_->push_back(Teuchos::getDoubleParameter(params,"CUTOFFRADIUS"));
    maxtreedepth_ = params.get<int>("BEAMPOT_TREEDEPTH", 6);
    // max number of bounding boxes per leaf octant
    minbboxesinoctant_ = params.get<int>("BEAMPOT_BOXESINOCT",8);

    btsph_ = DRT::INPUT::IntegralValue<int>(params, "BEAMPOT_BTSPH");
    btsol_ = DRT::INPUT::IntegralValue<int>(params, "BEAMPOT_BTSOL");
  }
  else
  {
    dserror("OctTree called with unknown type of parameter list!");
  }

  // sanity checks for extrusion value(s)
  if((int)extrusionvalue_->size()>2)
    dserror("You gave %i values for BEAMS_EXTVAL! Check your input file.", (int)extrusionvalue_->size());
  if((int)extrusionvalue_->size()==1)
    extrusionvalue_->push_back(extrusionvalue_->at(0));
  for(int i=0; i<(int)extrusionvalue_->size(); i++)
    if(extrusionvalue_->at(i)<1.0 && !additiveextrusion_)
      dserror("BEAMS_EXTVAL( %i ) = %4.2f < 1.0 does not make any sense! Check your input file.", i, extrusionvalue_->at(i));
  if(boundingbox_==Beam3ContactOctTree::cyloriented && (int)extrusionvalue_->size()!= 2)
    dserror("For cylindrical, oriented bounding boxes, TWO extrusion factors are mandatory! Check BEAMS_EXTVAL in your input file.");

  // set flag signaling the existence of periodic boundary conditions
  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // retrieve the dimensions of the periodic boundary box
  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(statmechparams,"PERIODLENGTH"));
    std::string word;
    char* input;
    while (PL >> word)
      periodlength_->push_back(std::strtod(word.c_str(), &input));
  }
  if((int)periodlength_->size()<3 && (int)periodlength_->size()!=1)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());
  else if((int)periodlength_->size()==1)
    for(int i=0; i<2; i++)
      periodlength_->push_back(periodlength_->at(0));

  if(periodlength_->at(0)>1e-12)
    periodicBC_ = true;
  else
    periodicBC_ = false;

  // determine bounding box type
  switch(bboxtype_input)
  {
    case INPAR::BEAMCONTACT::boct_aabb:
    {
      if(!discret_.Comm().MyPID())
        std::cout<<"Search routine:\nOctree + Axis Aligned BBs"<<std::endl;
      boundingbox_ = Beam3ContactOctTree::axisaligned;
    }
      break;
    case INPAR::BEAMCONTACT::boct_cobb:
    {
      if(!discret_.Comm().MyPID())
        std::cout<<"Search routine:\nOctree + Cylindrical Oriented BBs"<<std::endl;
      boundingbox_ = Beam3ContactOctTree::cyloriented;

      if (btsol_)
        dserror("Only axis aligned bounding boxes possible for beam to solid contact!");
      else if (btsph_)
        dserror("Only axis-aligned or spherical bounding boxes possible for beam-to-sphere contact!");
    }
      break;
    case INPAR::BEAMCONTACT::boct_spbb:
    {
      if(!discret_.Comm().MyPID())
        std::cout<<"Search routine:\nOctree + Spherical BBs"<<std::endl;
      boundingbox_ = Beam3ContactOctTree::spherical;

      if (btsol_)
        dserror("Only axis aligned bounding boxes possible for beam to solid contact!");
    }
      break;
    default: dserror("No octree (i.e. none) declared in your input file!");
      break;
  }

  if(!discret_.Comm().MyPID())
  {
    int numextval = (int)extrusionvalue_->size();
    if(additiveextrusion_)
    {
      numextval = 1;
      std::cout<<"additive extrusion     = ";
    }
    else
      std::cout<<"multiplicat. extrusion = ";
    for(int i=0; i<numextval; i++)
      std::cout<<extrusionvalue_->at(i)<<" ";
    std::cout<<std::endl;
    std::cout<<"max. tree depth        = "<<maxtreedepth_<<"\nmax. BB per octant     = "<<minbboxesinoctant_<<std::endl;
  }

  // get line conditions
  bbox2line_ = Teuchos::rcp(new Epetra_Vector(*(searchdis_.NodeColMap())));
  bbox2line_->PutScalar(-1.0);
  std::vector<DRT::Condition*> lines;
  discret_.GetCondition("FilamentNumber", lines);

  for(int i=0; i<(int)lines.size(); i++)
    for(int j=0; j<(int)lines[i]->Nodes()->size(); j++)
    {
      if(searchdis_.NodeColMap()->LID(lines[i]->Nodes()->at(j))<0)
        dserror("ID of node on FILAMENT NUMBERS condition not found in discretization searchdis_. Check your FILAMENT NUMBERS condition!");
      else
      (*bbox2line_)[searchdis_.NodeColMap()->LID(lines[i]->Nodes()->at(j))] = lines[i]->GetInt("Filament Number");
    }

  return;
}

/*----------------------------------------------------------------------*
 |  calls the almighty Octtree (public)                    mueller 01/11|
 *----------------------------------------------------------------------*/
std::vector<std::vector<DRT::Element*> > Beam3ContactOctTree::OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, int step)
{
#ifdef OCTREEDEBUG
  double t_start = Teuchos::Time::wallTime();
#endif
  // initialize class vectors
  InitializeOctreeSearch();

#ifdef OCTREEDEBUG
  double t_01 = Teuchos::Time::wallTime();
#endif
  // build axis aligned bounding boxes
  CreateBoundingBoxes(currentpositions);

#ifdef OCTREEDEBUG
  double t_02 = Teuchos::Time::wallTime();
#endif
  // call recursive octtree build
  // clear vector for assigning bounding boxes to octants to be on the safe side before (re)assigning bounding boxes
  bool bboxesfound = locateAll();

#ifdef OCTREEDEBUG
  double t_03 = Teuchos::Time::wallTime();
  double t_04 = 0.0;
#endif
  // intersection checks
  std::vector<std::vector<DRT::Element*> > contactpairelements;
  if(bboxesfound)
  {
    BoundingBoxIntersection(currentpositions, contactpairelements);

#ifdef OCTREEDEBUG
   t_04 = Teuchos::Time::wallTime();
#endif
    // output
    OctreeOutput(contactpairelements , step);
  }

#ifdef OCTREEDEBUG
  if(!discret_.Comm().MyPID())
  {
    std::cout<<"Initialization     :\t\t"<<t_01-t_start<<" seconds"<<std::endl;
    std::cout<<"Bound. Box Creation:\t\t"<<t_02-t_01<<" seconds"<<std::endl;
    std::cout<<"Octree Search      :\t\t"<<t_03-t_02<<" seconds"<<std::endl;
    std::cout<<"Intersection       :\t\t"<<t_04-t_03<<" seconds"<<std::endl;
    std::cout<<"Octree output      :\t\t"<<Teuchos::Time::wallTime()-t_04<<" seconds"<<std::endl;
    std::cout<<"============================================================="<<std::endl;
    std::cout<<"Octree Search time :\t\t"<<Teuchos::Time::wallTime()-t_start<<" seconds"<<std::endl;
  }
#endif
  //return contactpairs;
  return contactpairelements;
}// OctTreeSearch()
/*----------------------------------------------------------------------*
 |  Return the octants to which this bounding box belongs               |
 |  (public)                                               mueller 01/11|
 *----------------------------------------------------------------------*/
std::vector<int> Beam3ContactOctTree::InWhichOctantLies(const int& thisBBoxID)
{
  std::vector<int> octants;
  octants.clear();
  int bboxcolid = searchdis_.ElementColMap()->LID(thisBBoxID);
  int i=0;
  while((*bbox2octant_)[i][bboxcolid]>=0.0)
  {
    octants.push_back((int)(*bbox2octant_)[i][bboxcolid]);
    i++;
  }
  return octants;
}
/*----------------------------------------------------------------------*
 |  Intersect the bounding boxes of a certain octant with a given       |
 |  bounding box (public)                                  mueller 01/11|
 *----------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectBBoxesWith(Epetra_SerialDenseMatrix& nodecoords, Epetra_SerialDenseMatrix& nodeLID)
{
  dserror("Not in use!");

  /* notes:
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
  Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits = Teuchos::rcp(new Epetra_SerialDenseMatrix(1,1));

  // build bounding box according to given type
  switch(boundingbox_)
  {
    case Beam3ContactOctTree::axisaligned:
      CreateAABB(nodecoords,0,bboxlimits);
    break;
    case Beam3ContactOctTree::cyloriented:
      CreateCOBB(nodecoords,0,bboxlimits);
    break;
    case Beam3ContactOctTree::spherical:
      CreateSPBB(nodecoords,0,bboxlimits);
    default: dserror("No or an invalid Octree type was chosen. Check your input file!");
    break;
  }

  // retrieve octants in which the bounding box with ID thisBBoxID is located
  std::vector<std::vector<int> > octants;
  octants.clear();
  // get the octants for two bounding boxe (element) GIDs adjacent to each given node LID
  for(int i=0; i<nodeLID.M(); i++)
    octants.push_back(InWhichOctantLies(searchdis_.lColNode((int)nodeLID(i,0))->Elements()[0]->Id()));

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
            std::vector<int> bboxinoct(2,-1);
            bboxinoct[0] = (int)(*bboxesinoctants_)[i][octants[ibox][oct]];
            /*check for adjacent nodes: if there are adjacent nodes, then, of course, there
             * the intersection test will turn out positive. We skip those cases.*/
            // note: bounding box IDs are equal to element GIDs
            bool sharednode = false;
            for(int j=0; j<searchdis_.gElement(bboxinoct[0])->NumNode(); j++)
            {
              for(int k=0; k<(int)nodeLID.M(); k++)
              {
                if(searchdis_.NodeColMap()->LID(searchdis_.gElement(bboxinoct[0])->NodeIds()[j])==(int)nodeLID(k,0))
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
                  intersection = IntersectionAABB(bboxinoct, bboxlimits);
                break;
                case Beam3ContactOctTree::cyloriented:
                  intersection = IntersectionCOBB(bboxinoct, bboxlimits);
                break;
                case Beam3ContactOctTree::spherical:
                  intersection = IntersectionSPBB(bboxinoct, bboxlimits);
                break;
                default: dserror("No or an invalid Octree type was chosen. Check your input file!");
                break;
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

/*-----------------------------------------------------------------------------------*
 |  Output of octants, bounding boxes and contact pairs (public)       mueller 01/12 |
 *----------------------------------------------------------------------------------.*/
void Beam3ContactOctTree::OctreeOutput(std::vector<std::vector<DRT::Element*> > cpairelements, int step)
{
  if(!discret_.Comm().MyPID() && step!=-1)
  {
    // active contact pairs
    if((int)cpairelements.size()>0)
    {
      //Print ContactPairs to .dat-file and plot with Matlab....................
      std::ostringstream filename;
      if(step!=-2)
        filename << "ContactPairs"<<std::setw(6) << std::setfill('0') << step <<".dat";
      else
        filename << "ContactPairsInit.dat" <<std::endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int i=0;i<(int)cpairelements.size();i++)
        myfile << ((cpairelements[i])[0])->Id() <<"  "<< ((cpairelements[i][1]))->Id() <<std::endl;
      fprintf(fp, myfile.str().c_str());
      fclose(fp);
    }
    // octant limits output
    if((int)octreelimits_.size()>0)
    {
      std::ostringstream filename;
      if(step!=-2)
        filename << "OctreeLimits"<<std::setw(6) << std::setfill('0') << step <<".dat";
      else
        filename << "OctreeLimitsInit.dat"<<std::endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<(int)octreelimits_.size(); u++)
      {
        for (int v=0; v<(int)octreelimits_[u].M(); v++)
          myfile << std::scientific<<octreelimits_[u](v)<<" ";
        myfile <<std::endl;
      }
      // root box
      for(int u=0; u<(int)rootbox_.M(); u++)
        myfile<<std::scientific<<rootbox_(u)<<" ";
      myfile<<std::endl;
      fprintf(fp, myfile.str().c_str());
      fclose(fp);

#ifdef OCTREEDEBUG
      for (int u=0; u<(int)octreelimits_.size(); u++)
        for (int v=0; v<(int)octreelimits_[u].M(); v++)
          if(v%2==0 && octreelimits_[u](v)<rootbox_(v) && fabs(octreelimits_[u](v)-rootbox_(v))>1e-8)
            dserror("Octant minimum %4.10f below root box minimum %4.10f",octreelimits_[u](v),rootbox_(v));
          else if(v%2==1 && octreelimits_[u](v)>rootbox_(v) && fabs(octreelimits_[u](v)-rootbox_(v))>1e-8)
            dserror("Octant maximum %4.10f above root box minimum %4.10f",octreelimits_[u](v),rootbox_(v));
#endif
    }
    // bounding box coords output
    if(allbboxes_!=Teuchos::null)
    {
      std::ostringstream filename;
      if(step!=-2)
        filename << "BoundingBoxCoords"<<std::setw(6) << std::setfill('0') << step <<".dat";
      else
        filename << "BoundingBoxCoordsInit.dat"<<std::endl;
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<allbboxes_->MyLength(); u++)
      {
        for (int v=0; v<allbboxes_->NumVectors(); v++)
          myfile << std::scientific<< std::setprecision(10)<<(*allbboxes_)[v][u] <<" ";
        myfile <<(*diameter_)[u] <<std::endl;
      }
      fprintf(fp, myfile.str().c_str());
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)               mueller 01/12|
 |  Initialize octree class vectors and specifications                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::InitializeOctreeSearch()
{
#ifdef OCTREEDEBUG
  if(!discret_.Comm().MyPID())
    std::cout<<"Searchdis: "<<searchdis_.ElementColMap()->NumMyElements()<<", Probdis: "<<discret_.NumGlobalElements()<<std::endl;
#endif
  // number of shifts across volume boundaries in case of periodic boundary conditions (for intersection optimization)
  if(periodicBC_)
    numshifts_ = Teuchos::rcp(new Epetra_Vector(*(searchdis_.ElementColMap()),true));

  //determine radius factor by looking at the absolute mean variance of a bounding box (not quite sure...)
  //beam diameter
  diameter_ = Teuchos::rcp(new Epetra_Vector(*(searchdis_.ElementColMap())));
  switch(boundingbox_)
  {
    // for this case, the diameter is calculated in CreateSPBB()
    case Beam3ContactOctTree::spherical:
      diameter_->PutScalar(0.0);
    break;
    default:
    {
      for(int i=0; i<searchdis_.ElementColMap()->NumMyElements(); i++)
      {
        DRT::Element* element = searchdis_.lColElement(i);
        const DRT::ElementType & eot = element->ElementType();

        (*diameter_)[i] = -1.0;

        if (eot == DRT::ELEMENTS::Beam3Type::Instance())
          (*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(element))->Izz()) / M_PI));
        else if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
          (*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element))->Izz()) / M_PI));
        else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
          (*diameter_)[i] = 2.0 * sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3eb*>(element))->Izz()) / M_PI));
        else if (eot == DRT::ELEMENTS::RigidsphereType::Instance())
          (*diameter_)[i] = 2.0 * (dynamic_cast<DRT::ELEMENTS::Rigidsphere*>(element))->Radius();
        //If we have a solid element, we don't need its diameter and can it set to zero:
        else if (!BEAMCONTACT::BeamElement(*element) and !BEAMCONTACT::RigidsphereElement(*element))
          (*diameter_)[i] = 0.0;
        // feasibility check
        if ((*diameter_)[i] < 0.0) dserror("ERROR: Did not receive feasible element radius.");
      }
    }
  }
  // storage of Bounding Boxes
  // (components 0,...,5 contain bounding box limits)
  // (components 6,...,23 contain bounding box limits in case of periodic boundary conditions):
  // a box may b subject to a boundary shift up to 3 times -> 4 bounding boxes -> up to 24 values + 1 bounding box ID
  // note: for spherical bounding boxes, we only need the center point. Hence numbboxcoords = 3.
  int numbboxcoords = 6;
  if(boundingbox_==Beam3ContactOctTree::spherical)
    numbboxcoords = 3;

  // number of shifts across periodic boundaries (std: no periodic boundary conditions -> maxnumshifts = 0)
  int maxnumshifts = 0;
  if(periodicBC_)
    maxnumshifts = 3;
  allbboxes_ = Teuchos::rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),(maxnumshifts+1)*numbboxcoords+1, true));

  return;
}

/*----------------------------------------------------------------------*
 |  Bounding Box creation function (private)               mueller 01/11|
 |  generates bounding boxes extended with factor 1.05                  |
 *----------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions)
{
#ifdef MEASURETIME
  double t_AABB = Teuchos::Time::wallTime();
#endif
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

      //store nodal positions into matrix coords
      LINALG::SerialDenseMatrix coord(3,2,true);

      if(BEAMCONTACT::BeamElement(*element))
      {
        for(int i=0; i<2; i++)
        {
          int gid = element->Nodes()[i]->Id();
          LINALG::Matrix<3,1> coord_aux = currentpositions[gid];
          for(int j=0; j<3; j++)
            coord(j,i) = coord_aux(j);
        }
      }
      else if (BEAMCONTACT::RigidsphereElement(*element))
      {
        int gid = element->Nodes()[0]->Id();
        LINALG::Matrix<3,1> coord_aux = currentpositions[gid];
        for(int j=0; j<3; j++)
        {
          // write sphere center coordinates into both columns, then the creation of bbs will work with the following code for beams
          coord(j,0) = coord_aux(j);
          coord(j,1) = coord_aux(j);
        }
        if (boundingbox_ == Beam3ContactOctTree::cyloriented)
          dserror("Only axis-aligned or spherical bounding boxes possible for beam-to-sphere contact!");
      }
      else
      {
        //for solid elements only axis-aligned bounding boxes are possible
        switch(boundingbox_)
        {
          case Beam3ContactOctTree::axisaligned:
            CalcCornerPos(element, currentpositions, coord);
          break;
          case Beam3ContactOctTree::cyloriented:
            dserror("Only axis aligned bounding boxes possible for beam to solid contact!");
          break;
          case Beam3ContactOctTree::spherical:
            dserror("Only axis aligned bounding boxes possible for beam to solid contact!");
          break;
          default: dserror("No or an invalid Octree type was chosen. Check your input file!");
          break;
        }
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
        case Beam3ContactOctTree::spherical:
          CreateSPBB(coord, elecolid);
        break;
        default: dserror("No or an invalid Octree type was chosen. Check your input file!");
        break;
      }
    }
  } //end for-loop which goes through all elements

  // Communicate diameters to all procs, in case of spherical bounding boxes, this has to be done separately, i.e., here
  if(boundingbox_==Beam3ContactOctTree::spherical)
  {
    Epetra_Vector diameterrow(*(searchdis_.ElementRowMap()),true);
    CommunicateVector(diameterrow,*diameter_);
  }
  // communication of findings
  Epetra_MultiVector allbboxesrow(*(searchdis_.ElementRowMap()),allbboxes_->NumVectors(),true);
  CommunicateMultiVector(allbboxesrow, *allbboxes_);

  if(periodicBC_)
  {
    Epetra_Vector numshiftsrow(*(searchdis_.ElementRowMap()),true);
    CommunicateVector(numshiftsrow, *numshifts_);
  }

#ifdef MEASURETIME
  double bbgentimelocal = Teuchos::Time::wallTime() - t_AABB;
  double bbgentimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&bbgentimelocal, &bbgentimeglobal, 1);

  if(!searchdis_.Comm().MyPID())
    std::cout << "\n\nBBox creation time:\t\t" << bbgentimeglobal<< " seconds"<<std::endl;
#endif

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create an Axis Aligned Bounding Box   (private)                           mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateAABB(Epetra_SerialDenseMatrix&              coord,
                                     const int&                             elecolid,
                                     Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    dserror("Not in use!");

  // factor by which the box is extruded in each dimension
  double extrusionvalue = GetBoundingBoxExtrusionValue();

  if(elecolid<0 || elecolid>=searchdis_.ElementColMap()->NumMyElements())
    dserror("Given Element Column Map ID is %d !", elecolid);

  int elegid = searchdis_.ElementColMap()->GID(elecolid);

  // make the bounding box continuous
  std::vector<int> cut(3,0);
  int numshifts = 0;
  if(periodicBC_)
    UndoEffectOfPeriodicBoundaryCondition(coord,cut, numshifts);

  //directional vector (of the non-interrupted bounding box/ beam element
  LINALG::Matrix<3,1> dir;
  for(int dof=0; dof<(int)dir.M(); dof++)
    dir(dof) = coord(dof,1) - coord(dof,0);

  // extrude bounding box in axial direction
  if(additiveextrusion_)
    dir.Scale(1.0+extrusionvalue);
  else // multiplicative extrusion
    dir.Scale(extrusionvalue);

  // extrude coords
  for(int dof=0; dof<coord.M(); dof++)
  {
    double tmpcoord = coord(dof,1);
    coord(dof,1) = coord(dof,0) + dir(dof);
    coord(dof,0) = tmpcoord     - dir(dof);
  }

  // first bounding box (independent of periodicBC_)
  //Calculate Center Point of AABB
  LINALG::Matrix<3,1> midpoint;
  for(int i=0; i<(int)midpoint.M(); i++)
    midpoint(i) = 0.5*(coord(i,0) + coord(i,1));

  //Calculate edgelength of AABB, select max value
  double maxedgelength = -1.0;
  for(int i=0; i<(int)midpoint.M(); i++)
    if(fabs(coord(i,1) - coord(i,0))>maxedgelength)
      maxedgelength = fabs(coord(i,1) - coord(i,0));
  //Check for edgelength of AABB
  if (maxedgelength<(*diameter_)[elecolid])
    maxedgelength = (*diameter_)[elecolid];

  // Calculate limits of AABB
  // first bounding box
  //mins
  (*allbboxes_)[0][elecolid] = midpoint(0) - 0.5*maxedgelength; // x
  (*allbboxes_)[2][elecolid] = midpoint(1) - 0.5*maxedgelength; // y
  (*allbboxes_)[4][elecolid] = midpoint(2) - 0.5*maxedgelength; // z
  // maxes
  (*allbboxes_)[1][elecolid] = midpoint(0) + 0.5*maxedgelength; // x
  (*allbboxes_)[3][elecolid] = midpoint(1) + 0.5*maxedgelength; // y
  (*allbboxes_)[5][elecolid] = midpoint(2) + 0.5*maxedgelength; // z

  // periodic shifts of bounding box
  if(numshifts>0)
  {
    int shift = 0;

    for(int dof=0; dof<(int)cut.size(); dof++)
    {
      if(cut[dof]>0)
      {
        double periodicshift = 0.0;
        switch(cut[dof])
        {
          case 1:
            periodicshift = (*periodlength_)[dof];
          break;
          case 2:
            periodicshift = -(*periodlength_)[dof];
          break;
        }
        for(int i=0; i<6; i++)
        {
          (*allbboxes_)[6*shift+i][elecolid] = (*allbboxes_)[i][elecolid];
          if(i%3==dof)
            (*allbboxes_)[6*shift+i][elecolid] += periodicshift;
        }
      }
    }
    (*numshifts_)[elecolid] = shift; // store number of shifts
  }

  if(periodicBC_ && numshifts<3)
    for(int i=allbboxes_->NumVectors()-1; i>(numshifts+1)*6-1; i--)
      (*allbboxes_)[i][elecolid] = -1e9;

  // store GID (=box number)
 (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Cylindrical an Oriented Bounding Box   (private)                   mueller 11/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateCOBB(Epetra_SerialDenseMatrix&              coord,
                                     const int&                             elecolid,
                                     Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
 if(bboxlimits!=Teuchos::null)
   dserror("Not in use!");
  // make the bounding box continuous
  std::vector<int> cut(3,0);
  int numshifts = 0;
  if(periodicBC_)
    UndoEffectOfPeriodicBoundaryCondition(coord,cut, numshifts);

  double extrusionvalue = GetBoundingBoxExtrusionValue();


  int elegid = searchdis_.ElementColMap()->GID(elecolid);

  //directional vector (of the non-interrupted bounding box/ beam element
  LINALG::Matrix<3,1> dir;
  for(int dof=0; dof<(int)dir.M(); dof++)
    dir(dof) = coord(dof,1) - coord(dof,0);

  // extrude bounding box in axial direction
  if(additiveextrusion_)
  {
    LINALG::Matrix<3,1> additivedir(dir);
    additivedir.Scale(extrusionvalue/dir.Norm2());
    dir += additivedir;
  }
  else // multiplicative extrusion
    dir.Scale(extrusionvalue);

  LINALG::Matrix<3,1> midpoint;
  for(int i=0; i<(int)midpoint.M(); i++)
    midpoint(i) = 0.5*(coord(i,0) + coord(i,1));

  // First bounding box is the untreated set of coordinates!
  for(int dof=0; dof<coord.M(); dof++)
  {
    (*allbboxes_)[dof][elecolid]   = midpoint(dof) - 0.5*dir(dof);
    (*allbboxes_)[dof+3][elecolid] = midpoint(dof) + 0.5*dir(dof);
  }

  int shift = 0;
  if(numshifts>0)
  {
    // create shifted copies of the bounding boxes if there is an interruption by periodic BCs
    for(int dof=0; dof<(int)cut.size(); dof++)
    {
      if(cut[dof]>0)
      {
        shift++;
        double periodicshift = 0.0;
        switch(cut[dof])
        {
          case 1:
            periodicshift = (*periodlength_)[dof];
          break;
          case 2:
            periodicshift = -(*periodlength_)[dof];
          break;
        }
        // (node) position 1
        (*allbboxes_)[6*shift+dof][elecolid] = (*allbboxes_)[dof][elecolid] + periodicshift;
        // (node) position 2
        (*allbboxes_)[6*shift+dof+3][elecolid] = (*allbboxes_)[dof+3][elecolid] + periodicshift;

        // unshifted DOFs
        // (node) position 1
        (*allbboxes_)[6*shift+((dof+1)%3)][elecolid] = (*allbboxes_)[(dof+1)%3][elecolid];
        (*allbboxes_)[6*shift+((dof+2)%3)][elecolid] = (*allbboxes_)[(dof+2)%3][elecolid];
        // (node) position 2
        (*allbboxes_)[6*shift+((dof+1)%3) + 3][elecolid] = (*allbboxes_)[(dof+1)%3 + 3][elecolid];
        (*allbboxes_)[6*shift+((dof+2)%3) + 3][elecolid] = (*allbboxes_)[(dof+2)%3 + 3][elecolid];
      }
    }
    // fill all latter entries  except for the last one (->ID) with bogus values (in case of periodic BCs)
    (*numshifts_)[elecolid] = numshifts;

//    std::cout<<"Shifted BBOX "<<elecolid<<": ";
//    for(int i=0; i<allbboxes_->NumVectors(); i++)
//      std::cout<<(*allbboxes_)[i][elecolid]<<" ";
//    std::cout<<std::endl;
  }
  if(periodicBC_ && numshifts<3)
    for(int i=allbboxes_->NumVectors()-1; i>(numshifts+1)*6-1; i--)
      (*allbboxes_)[i][elecolid] = -1e9;

  // last entry: element GID
  (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  Create Spherical Bounding Box   (private)                                  mueller 1/12|
 *-----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateSPBB(Epetra_SerialDenseMatrix&              coord,
                                     const int&                             elecolid,
                                     Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    dserror("Not in use!");
  // make the bounding box continuous
  std::vector<int> cut(3,0);
  int numshifts = 0;
  if(periodicBC_)
    UndoEffectOfPeriodicBoundaryCondition(coord,cut, numshifts);

  //Note: diameter_ is used as the sphere diameter in the context of SPBBs
  double extrusionvalue = GetBoundingBoxExtrusionValue();

  int elegid = searchdis_.ElementColMap()->GID(elecolid);

  // get the element with local ID (LID) elecolid
  DRT::Element* element = searchdis_.lColElement(elecolid);
  double diameter =0.0;

  if (BEAMCONTACT::BeamElement(*element))
  {
    if (coord.M()==3 and coord.N()==2)
    {
      diameter = sqrt((coord(0,0)-coord(0,1))*(coord(0,0)-coord(0,1))+
                      (coord(1,0)-coord(1,1))*(coord(1,0)-coord(1,1))+
                      (coord(2,0)-coord(2,1))*(coord(2,0)-coord(2,1)));
    }
    else
      dserror("coord matrix of nodal positions has wrong dimensions here!");
  }
  else if (BEAMCONTACT::RigidsphereElement(*element))
    diameter = (*diameter_)[elecolid];

  // we want an additive extrusion of the radius, hence "2.0"
  if(additiveextrusion_)
    (*diameter_)[elecolid] = diameter + 2.0*extrusionvalue;
  else
    (*diameter_)[elecolid] = diameter * extrusionvalue;

  // first bounding box
  for(int dof=0; dof<coord.M(); dof++)
    (*allbboxes_)[dof][elecolid] = 0.5*(coord(dof,0) + coord(dof,1));

  bool periodicbbox = false;
  int shift = 0;
  if(numshifts>0)
  {
    for(int dof=0; dof<coord.M(); dof++)
    {
      for(int pos=0; pos<coord.N(); pos++)
      {
        if(coord(dof,pos)>=(*periodlength_)[dof]- 0.5*(*diameter_)[elecolid] || coord(dof,pos)<=0.5*(*diameter_)[elecolid])
        {
          shift++;
          // shift downwards
          if(coord(dof,pos)>=(*periodlength_)[dof]- 0.5*(*diameter_)[elecolid])
          {
            (*allbboxes_)[3*shift+dof][elecolid]       = (*allbboxes_)[dof][elecolid]-(*periodlength_)[dof];
            (*allbboxes_)[3*shift+(dof+1)%3][elecolid] = (*allbboxes_)[(dof+1)%3][elecolid];
            (*allbboxes_)[3*shift+(dof+2)%3][elecolid] = (*allbboxes_)[(dof+2)%3][elecolid];
            // no need to check the second position if pos==0
            periodicbbox = true;
            break;
          }
          // shift upwards
          else if(coord(dof,pos)<=0.5*(*diameter_)[elecolid])
          {
            (*allbboxes_)[3*shift+dof][elecolid]       = (*allbboxes_)[dof][elecolid] + (*periodlength_)[dof];
            (*allbboxes_)[3*shift+(dof+1)%3][elecolid] = (*allbboxes_)[(dof+1)%3][elecolid];
            (*allbboxes_)[3*shift+(dof+2)%3][elecolid] = (*allbboxes_)[(dof+2)%3][elecolid];
            // no need to check the second position if pos==0
            periodicbbox = true;
            break;
          }
        }
      }
      if(periodicbbox)
        break;
    }
    (*numshifts_)[elecolid] = shift; // store number of shifts
  }

  if(periodicBC_ && shift<3)
    for(int i=allbboxes_->NumVectors()-1; i>(shift+1)*3-1; i--)
      (*allbboxes_)[i][elecolid] = -1e9;

  // last entry: element GID
  (*allbboxes_)[allbboxes_->NumVectors()-1][elecolid] = elegid;
  return;
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
bool Beam3ContactOctTree::locateAll()
{
#ifdef MEASURETIME
  double t_octree = Teuchos::Time::wallTime();
#endif

  numcrit1_ = 0;
  numcrit2_ = 0;

  // Convert Epetra_MultiVector allbboxes_ to vector(std::vector<double>)
//  std::cout<<allbboxes_->MyLength()<<", "<<allbboxes_->NumVectors()<<std::endl;
  std::vector<std::vector<double> > allbboxesstdvec(allbboxes_->MyLength(), std::vector<double>(allbboxes_->NumVectors(),0.0));
  EpetraMultiVecToStdVec(*allbboxes_,allbboxesstdvec);

  // initial tree depth value (will be incremented with each recursive call of locateBox()
  int treedepth = 0;

  // Parameters and Initialization
  std::vector<std::vector<int> > bboxesinoctants;
  bboxesinoctants.clear();
  octreelimits_.clear();

  // explanation: bbox2octant[i][j]: bounding box i has a vector of octants j, to which the bounding box belongs
  std::vector<std::vector<int> > bbox2octant(allbboxes_->MyLength(), std::vector<int>());

  // Recursively construct octree; Proc 0 only (parallel computing impossible)
  // max depth of OctreeMap
  int maxdepthlocal = 0;
  int bboxlengthlocal = 0;
  // max depth of octree
  int maxnumoctlocal =  0;
  int maxnumoctglobal = 0;
  int maxdepthglobal = 0;
  int bboxlengthglobal = 0;


  // calculate root box limits
  rootbox_ = GetRootBox();
  if(searchdis_.Comm().MyPID()==0)
  {
    locateBox(allbboxesstdvec, rootbox_, octreelimits_, bboxesinoctants, bbox2octant, treedepth);

    for(int i=0; i<(int)bbox2octant.size(); i++)
      if((int)bbox2octant[i].size()>maxnumoctlocal)
        maxnumoctlocal = (int)bbox2octant[i].size();

    bboxlengthlocal = (int)bboxesinoctants.size();
    for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
      if((int)bboxesinoctants[i].size()>maxdepthlocal)
        maxdepthlocal = (int)bboxesinoctants[i].size();
  }

  // communicate bbox2octant_ to all Procs
  searchdis_.Comm().MaxAll(&maxnumoctlocal, &maxnumoctglobal, 1);
  bbox2octant_ = Teuchos::rcp(new Epetra_MultiVector(*(searchdis_.ElementColMap()),maxnumoctglobal,true));

  // fill epetra vector
  if(!searchdis_.Comm().MyPID())
  {
    bbox2octant_->PutScalar(-9.0);
    StdVecToEpetraMultiVec(bbox2octant,*bbox2octant_);
  }

  Epetra_MultiVector bbox2octantrow(*(searchdis_.ElementRowMap()),maxnumoctglobal,true);
  CommunicateMultiVector(bbox2octantrow, *bbox2octant_);

  searchdis_.Comm().MaxAll(&maxdepthlocal, &maxdepthglobal, 1);
  searchdis_.Comm().MaxAll(&bboxlengthlocal, &bboxlengthglobal, 1);

  /* build temporary, fully overlapping map and row map for octree
   * Note: maxdepthglobal does not occur for a converging Newton iteration. Yet, in some cases, when
   * encountering divergence for the Newton scheme, this might happen.
   * In biopolymer network simulations, this setting is not unlikely and unavoidable. a maximum depth of
   * 0 means, there are no bounding boxes/elements in any octants. Hence, we will not detect any contact and
   * therefore skip the rest of the octree algorithm.*/
  if(maxdepthglobal>0)
  {
    // create octree maps
    std::vector<int> gids;
    for (int i=0 ; i<bboxlengthglobal; i++ )
      gids.push_back(i);
    // crosslinker column and row map
    Epetra_Map octtreerowmap((int)gids.size(), 0, discret_.Comm());
    Epetra_Map octtreemap(-1, (int)gids.size(), &gids[0], 0, discret_.Comm());

    // build Epetra_MultiVectors which hold the BBs of the OctreeMap; for communication
    bboxesinoctants_ = Teuchos::rcp(new Epetra_MultiVector(octtreemap,maxdepthglobal),true);
    Epetra_MultiVector bboxinoctrow(octtreerowmap,maxdepthglobal, true);

    // fill bboxinoct for Proc 0
    if(searchdis_.Comm().MyPID()==0)
    {
      bboxesinoctants_->PutScalar(-9.0);
//      for (int i=0 ; i<(int)bboxesinoctants.size(); i++ )
//        for(int j=0; j<(int)bboxesinoctants[i].size(); j++)
//          (*bboxesinoctants_)[j][i] = bboxesinoctants[i][j];
      StdVecToEpetraMultiVec(bboxesinoctants,*bboxesinoctants_);
    }

    // Communication
    CommunicateMultiVector(bboxinoctrow, *bboxesinoctants_);

#ifdef OCTREEDEBUG
    std::ostringstream filename;
    if(!discret_.Comm().MyPID())
    {
      filename << "BBinOct_"<<searchdis_.Comm().NumProc()<<"_procs.dat";
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "w");
      std::stringstream myfile;
      for (int u=0; u<bboxesinoctants_->MyLength(); u++)
      {
        for (int v=0; v<bboxesinoctants_->NumVectors(); v++)
          myfile << std::scientific<<(*bboxesinoctants_)[v][u] <<" ";
        myfile <<std::endl;
      }
      fprintf(fp, myfile.str().c_str());
      fclose(fp);
    }
#endif
#ifdef MEASURETIME
     if(!searchdis_.Comm().MyPID())
       std::cout << "\nOctree building time:\t\t" << Teuchos::Time::wallTime() - t_octree<< " seconds" << std::endl;
#endif
     return true;
  }
  else
  {
    return false;
  }
}// end of method locateAll


/*-----------------------------------------------------------------------------------------*
 |  locateBox function (private);                                               meier 02/11|
 *----------------------------------------------------------------------------------------*/
void Beam3ContactOctTree::locateBox(std::vector<std::vector<double> >& allbboxesstdvec,
                                    LINALG::Matrix<6,1>&               lim,
                                    std::vector<LINALG::Matrix<6,1> >& octreelimits,
                                    std::vector<std::vector<int> >&    bboxesinoctants,
                                    std::vector<std::vector<int> >&    bbox2octant,
                                    int& treedepth)
{
  double extrusionvalue = GetBoundingBoxExtrusionValue();

  // edge length vector of the suboctants
  LINALG::Matrix<3,1> newedgelength;
  // vector with limits of the 8 sub octants
  std::vector<LINALG::Matrix<6,1> > suboctlimits;
  CreateSubOctants(lim,newedgelength,suboctlimits);

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
  Teuchos::RCP<LINALG::Matrix<3,1> > octcenter = Teuchos::null;

  //Goes through all suboctants
  for( int oct=0; oct<8; oct++)
  {
    // Define temporary vector of same size as current allbboxesstdvec
    std::vector<std::vector<double> > bboxsubset;
    bboxsubset.clear();

    // we need the octant centers when applying cylindrical or spherical bounding boxes.
    if(boundingbox_==Beam3ContactOctTree::cyloriented || boundingbox_==Beam3ContactOctTree::spherical)
    {
      octcenter = Teuchos::rcp(new LINALG::Matrix<3,1>);
      for(int i=0; i<(int)octcenter->M(); i++)
        (*octcenter)(i) = 0.5*(suboctlimits[oct](2*i)+suboctlimits[oct](2*i+1));
    }

    int numshifts = 0;

    for( int i=0; i<(int)allbboxesstdvec.size(); i++)
    {
      int lid = searchdis_.ElementColMap()->LID((int)allbboxesstdvec[i][allbboxesstdvec[i].size()-1]);
      if(periodicBC_)
        numshifts = (*numshifts_)[lid];

      for(int shift=0; shift<numshifts+1; shift++)
      {
        // Processes columns indices 1 to 6
        // 2)loop over the limits of the current octant and check if the current bounding box lies within this octant.
        // 3)Then, check component-wise and leave after first "hit"
        bool assignbboxtooctant = false;
        switch(boundingbox_)
        {
          case Beam3ContactOctTree::axisaligned:
          {
            if(AABBIsInThisOctant(suboctlimits[oct], allbboxesstdvec[i], shift))
              assignbboxtooctant = true;
          }
          break;
          case Beam3ContactOctTree::cyloriented:
          {
//            std::cout<<"Octant-Lim:"<<suboctlimits[oct]<<std::endl;
            if(COBBIsInThisOctant(*octcenter, newedgelength,allbboxesstdvec[i],extrusionvalue,lid,shift))
              assignbboxtooctant = true;
          }
          break;
          case Beam3ContactOctTree::spherical:
          {
            if(SPBBIsInThisOctant(*octcenter, newedgelength, allbboxesstdvec[i], lid, shift))
              assignbboxtooctant = true;
          }
          break;
          default: dserror("No or an invalid Octree type was chosen. Check your input file!");
          break;
        }
        if(assignbboxtooctant)
        {
          bboxsubset.push_back(allbboxesstdvec[i]);
          break;
        }
      }
    }// end of for-loop which goes through all elements of input

    // current tree depth
    int currtreedepth = treedepth+1;
    // Check for further recursion by checking number of boxes in octant (first criterion)....................
    int N = (int)bboxsubset.size();

    //If to divide further, let LocateBox call itself with updated inputs
    if (N > minbboxesinoctant_ && currtreedepth < maxtreedepth_-1)
      locateBox(bboxsubset, suboctlimits[oct], octreelimits, bboxesinoctants, bbox2octant, currtreedepth);
    else
    {
      // no further discretization of the volume because either the maximal tree depth or the minimal number of bounding
      // boxes per octant has been reached
      // this vector holds the IDs of the bounding boxes in this octant
      if(N>0)
      {
        std::vector<int> boxids;
        boxids.clear();
        //Push back Limits of suboctants to octreelimits
        octreelimits.push_back(suboctlimits[oct]);

        for (int m = 0; m < (int)bboxsubset.size(); m++)
        {
          // note: the Bounding Box ID is the last column entry of the m-th entry vector bboxsubset
          int bboxgid = (int)bboxsubset[m][(int)bboxsubset[m].size()-1];
          boxids.push_back(bboxgid);
          // assign current octant number to the bounding box
          int collid = searchdis_.ElementColMap()->LID(bboxgid);
          //std::cout<<"Collid = "<<collid<<std::endl;
          bbox2octant[collid].push_back((int)bboxesinoctants.size());
        }
        // add bounding box IDs of this octant to the global vector
        bboxesinoctants.push_back(boxids);

//        std::cout<<"Octant "<<bboxesinoctants.size()<<", treedepth = "<<currtreedepth<<": ";
//        for(int i=0; i<(int)bboxesinoctants[bboxesinoctants.size()-1].size(); i++)
//          std::cout<<bboxesinoctants[bboxesinoctants.size()-1][i]<<" ";
//        std::cout<<", num of bboxes: "<<bboxsubset.size()<<std::endl;
      }
    }
  }// end of loop which goes through all suboctants
  return;
} // end of method locateBox

/*----------------------------------------------------------------------------------*
 |  Calculate limits of the root box   (public)                        mueller 11/11|
 *----------------------------------------------------------------------------------*/
LINALG::Matrix<6,1> Beam3ContactOctTree::GetRootBox()
{
  int entriesperbbox = 0;
  if(boundingbox_ == Beam3ContactOctTree::spherical)
    entriesperbbox = 3;
  else
    entriesperbbox = 6;

  LINALG::Matrix<6,1> lim(true);
  // determine globally extremal coordinates and use them as root box.
  // initialize
  lim(0) = 1.0e9;
  lim(1) = -1.0e9;
  // loop over allbboxes_ and determine the extremes

  // number of shifts due to periodic BCs (if present, else 0)
  int numshifts = 0;
  for(int i=0; i<allbboxes_->MyLength(); i++)
  {
    if(periodicBC_)
      numshifts = (int)(*numshifts_)[i];
    for(int j=0; j<entriesperbbox*numshifts+entriesperbbox; j++)
    {
      if((*allbboxes_)[j][i]<lim(0)) // mins
        lim(0) = (*allbboxes_)[j][i];
      if((*allbboxes_)[j][i]>lim(1)) // maxes
        lim(1) = (*allbboxes_)[j][i];
    }
  }
  lim(2) = lim(0);
  lim(4) = lim(0);
  lim(3) = lim(1);
  lim(5) = lim(1);
  // safety factor
  lim.Scale(1.5);

  return lim;
}
/*-----------------------------------------------------------------------------------*
 |  Create limits of sub octants      (private)                        mueller 02/15|
 |----------------------------------------------------------------------------------*/
void Beam3ContactOctTree::CreateSubOctants(LINALG::Matrix<6,1>&               parentoctlimits,
                                           LINALG::Matrix<3,1>&               suboctedgelength,
                                           std::vector<LINALG::Matrix<6,1> >& suboctlimits)
{
  // Center of parent octant
  LINALG::Matrix<3,1> parentcenter;

  for(int i=0; i<(int)parentcenter.M(); i++)
  {
    parentcenter(i) = 0.5*(parentoctlimits(2*i)+parentoctlimits(2*i+1));
    suboctedgelength(i) = 0.5*fabs(parentoctlimits(2*i+1)-(parentoctlimits)(2*i));
  }

  // safety: only cubic octants for now.
  if(fabs(suboctedgelength(0)-suboctedgelength(1))>1e-8 || fabs(suboctedgelength(0)-suboctedgelength(2))>1e-8)
    dserror("Check GetRootBox()! Octants must be cubic!");

  suboctlimits.clear();
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<2; j++)
    {
      for(int k=0; k<2; k++)
      {
        LINALG::Matrix<6,1> sublim;
        sublim(0) = parentcenter(0) + (i-1)*suboctedgelength(0);
        sublim(1) = parentcenter(0) +  i   *suboctedgelength(0);
        sublim(2) = parentcenter(1) + (j-1)*suboctedgelength(1);
        sublim(3) = parentcenter(1) +  j   *suboctedgelength(1);
        sublim(4) = parentcenter(2) + (k-1)*suboctedgelength(2);
        sublim(5) = parentcenter(2) +  k   *suboctedgelength(2);

        suboctlimits.push_back(sublim);
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------------------*
 |  Check if AABB in octant           (private)                        mueller 02/15|
 |----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::AABBIsInThisOctant(LINALG::Matrix<6,1>& suboctlimits,
                                             std::vector<double>& bboxcoords,
                                             int &                shift)
{
  if(!((suboctlimits(0) >= bboxcoords[6*shift+1]) ||
      (bboxcoords[6*shift] >= suboctlimits(1))    ||
      (suboctlimits(2) >= bboxcoords[6*shift+3])  ||
      (bboxcoords[6*shift+2] >= suboctlimits(3))  ||
      (suboctlimits(4) >= bboxcoords[6*shift+5])  ||
      (bboxcoords[6*shift+4] >= suboctlimits(5))))
    return true;
  else
    return false;
}

/*-----------------------------------------------------------------------------------*
 |  Check if COBB in octant           (private)                        mueller 02/15|
 |----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::COBBIsInThisOctant(LINALG::Matrix<3,1>& octcenter,
                                             LINALG::Matrix<3,1>& newedgelength,
                                             std::vector<double>& bboxcoords,
                                             double&              extrusionvalue,
                                             int&                 lid,
                                             int &                shift)
{
  /* General idea:
   * Two checks ensure whether a bounding box lies within the given octant
   * 1) Check the two end points of the bounding box
   * 2) Check the linear connection between the end points for an intersection with the octant faces.
   * @1) Evaluate the closest intersection point of a line going through the octant center and end points.
   * If the distance to the end point is smaller than the distance between octant center and the intersection with the nearest
   * octant face plane-> intersection
   * @2) Calculate the intersection points of the line connecting the two end point of the bounding box with the octant face planes.
   * If the signed distance parameter lambda_k lies within the boundaries of the bounding box -> intersection.
   */
  // increase the radius of the box in order to have some additional buffer volume.
  double boxradius = 0.0;
  if(additiveextrusion_)
    boxradius = extrusionvalue + 0.5*(*diameter_)[lid];
  else
    boxradius = 0.5*extrusionvalue*(*diameter_)[lid];

  // criterion I: distance
  // check for both bounding box end points
  LINALG::Matrix<3,1> dir;
  for(int j=0; j<2; j++)
  {
    // std::vector to LINALG::Matrix (slower but nicer to read)
    LINALG::Matrix<3,1> bboxendpoint;
    bboxendpoint(0) = bboxcoords[6*shift+3*j];
    bboxendpoint(1) = bboxcoords[6*shift+3*j+1];
    bboxendpoint(2) = bboxcoords[6*shift+3*j+2];

    dir = bboxendpoint;
    dir -= octcenter;
    double dist = dir.Norm2();
    dir.Scale(1.0/dist);

    // find shortest distance to face planes
    LINALG::Matrix<3,1> lambda;
    for(int k=0; k<(int)lambda.M(); k++)
    {
      double normal = 0.0;
      if(bboxendpoint(k)>octcenter(k))
        normal = 1.0;
      else if(bboxendpoint(k)<octcenter(k))
        normal = -1.0;
      else // no component in direction k -> no intersection -> sinularity -> lambda ->infty -> insert high enough dummy value
      {
        lambda(k) = 1.0e9;
        break;
      }
      lambda(k) = 0.5*(newedgelength(k)*normal)/(dir(k));
    }
    // closest intersection point with face plane
    double lambda_min = lambda.MinValue();

    //only one of the endpoint must lie in the octant
    if(lambda_min+boxradius>=dist)
      return true;
  }

  // Criterion II: some part between the two given end points of the bbox may lie within the octant
  // determine unit ofrientation for line equation
  LINALG::Matrix<3,1> bboxorient;
  for(int k=0; k<(int)bboxorient.M(); k++)
    bboxorient(k) = bboxcoords[6*shift+k+3]-bboxcoords[6*shift+k];
  double bboxlength = bboxorient.Norm2();
//    std::cout<<"bboxlentgth = "<<bboxlength<<std::endl;
  bboxorient.Scale(1.0/bboxorient.Norm2());
  // component-wise check if the intersection point of the line lies outside or inside the octant (using an outward-facing normal)
  for(int k=0; k<(int)bboxorient.M(); k++)
  {
    // no intersection will occur since bbox is not slanted towards the face
    if(bboxorient(k) == 0.0)
      continue;

    // retrieve the correct value for the outward-facing normal
//      std::cout<<"bboxcoords[6*shift+k,"<<6*shift+k<<"] = "<<bboxcoords[6*shift+k]<<std::endl;
    double normal = 0.0;
    if(bboxcoords[6*shift+k]>=octcenter(k))
      normal = 1.0;
    else if(bboxcoords[6*shift+k]<octcenter(k))
      normal = -1.0;

    // signed distance from the first endpoint to the intersection point
    double lambda_k = (octcenter(k) + 0.5*newedgelength(k)*normal - bboxcoords[6*shift+k])/bboxorient(k);
    // calculate intersection point coords and subsequently the distance from octcenter to the intersecion
    //(only for promising candidates)
    if(lambda_k>-boxradius && lambda_k<bboxlength+boxradius)
    {
      LINALG::Matrix<3,1> dist;
      for(int l=0; l<(int)dist.M(); l++)
        dist(l) = bboxcoords[6*shift+l] + lambda_k*bboxorient(l)-octcenter(l);
      // bbox is in octant if one intersection point lies on the octant surface.
      if(dist.Norm2()<=0.5*newedgelength(k)*sqrt(3.0)+boxradius)
      {
//        if((int)(*numshifts_)[lid]>0)
//        {
//          LINALG::Matrix<3,1> aux(dist);
//          aux += octcenter;
//          std::cout<<"BBOX "<<lid<<": ";
//          for(int k=0; k<(int)bboxcoords.size(); k++)
//            if(bboxcoords[k]>-1e9)
//              std::cout<<bboxcoords[k]<<" ";
//          std::cout<<std::endl;
//          std::cout<<"; \nISEC : "<<aux(0)<<" "<<aux(1)<<" "<<aux(2)<<"\nLambda = "<<lambda_k<<", normal = "<<normal<<std::endl;
//          dserror("CRIT II");
//        }
        return true;
      }
    }
  }
  return false;
}

/*-----------------------------------------------------------------------------------*
 |  Check if SPBB in octant           (private)                        mueller 02/15|
 |----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::SPBBIsInThisOctant(LINALG::Matrix<3,1>& octcenter,
                                             LINALG::Matrix<3,1>& newedgelength,
                                             std::vector<double>& bboxcoords,
                                             int &                lid,
                                             int &                shift)
{
  double maxdistx = 0.5*(newedgelength(0)+(*diameter_)[lid]);
  double maxdisty = 0.5*(newedgelength(1)+(*diameter_)[lid]);
  double maxdistz = 0.5*(newedgelength(2)+(*diameter_)[lid]);

  double distx = fabs(octcenter(0)-bboxcoords[3*shift]);
  double disty = fabs(octcenter(1)-bboxcoords[3*shift+1]);
  double distz = fabs(octcenter(2)-bboxcoords[3*shift+2]);
  // check: component-wise distance
  if(distx <= maxdistx && disty <= maxdisty && distz <= maxdistz)
    return true;
  else
    return false;
}

/*-----------------------------------------------------------------------------------*
 |  Bounding Box Intersection function (private)                          meier 01/11|
 |  Intersects Bounding Boxes in same line of map OctreeMap                          |
 |  Gives back vector of intersection pairs                                          |
 *----------------------------------------------------------------------------------*/
void Beam3ContactOctTree::BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
                                                  std::vector<std::vector<DRT::Element*> >& contactpairelements)
{
#ifdef MEASURETIME
  double t_search = Teuchos::Time::wallTime();
#endif
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

          //Here we have introduced an criteria in order to sort neighboring contact element pairs out.
          //This method is based on the assumption, that two elements sharing the same node should not get into contact. In contrary to
          //a former criterion based on filament numbers (which forbids self contact), this method works for arbitrary element types and still allows for self contact!!!
          if(BEAMCONTACT::ElementsShareNode(*element1,*element2))
            considerpair = false;
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
            case Beam3ContactOctTree::spherical:
              intersection = IntersectionSPBB(bboxIDs);
            break;
            default: dserror("No or an invalid Octree type was chosen. Check your input file!");
            break;
          }

          if (intersection)
          {
            // note: creation of unique "first" entries in map, attention: IDs identical to crosslinker GIDs!!
            int mapfirst = (bboxIDs[0] + 1)*basisnodes_ + bboxIDs[1];
            contactpairmap.insert ( std::pair<int, std::vector<int> > (mapfirst, bboxIDs));
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
      //std::cout << std::scientific << (*it).first <<"  "<< ((*it).second)[0]<<" "<< ((*it).second)[1]<<std::endl;
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
    std::vector<DRT::Element*> temp_vec(2);
    temp_vec[0]=tempele1;
    temp_vec[1]=tempele2;
    contactpairelements.push_back(temp_vec);
  }
  //if(!discret_.Comm().MyPID())
    //std::cout<<"number of boxes: "<<counter<<std::endl;

#ifdef MEASURETIME
  double isectimelocal = Teuchos::Time::wallTime() - t_search;
  double isectimeglobal = 0.0;

  searchdis_.Comm().MaxAll(&isectimelocal, &isectimeglobal, 1);
  discret_.Comm().Barrier();
  if(!searchdis_.Comm().MyPID())
    std::cout << "Intersection time:\t\t" << isectimeglobal << " seconds"<<std::endl;
#endif

  return;
}//end of method BoundingBoxIntersection()

/*-----------------------------------------------------------------------------------*
 |  Axis Aligned Bounding Box Intersection function when both bounding boxes         |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionAABB(const std::vector<int>& bboxIDs, Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    dserror("Not in use!");

  // translate box / element GIDs to ElementColMap()-LIDs
  // note: GID and ColumnMap LID are usually the same except for crosslinker elements from statistical mechanics
  int entry1 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int entry2 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  //Initialization....................
  double a_xmin, a_xmax, a_ymin, a_ymax, a_zmin, a_zmax;
  double b_xmin, b_xmax, b_ymin, b_ymax, b_zmin, b_zmax;

  int I = 0;
  int J = 0;

  if(periodicBC_)
  {
    I = (int)(*numshifts_)[entry1];
    J = (int)(*numshifts_)[entry2];
  }
  // note: n shifts means n+1 segments
  for(int i=0; i<I+1; i++)
  {
    for(int j=0; j<J+1; j++)
    {
      //Intersection Test
      a_xmin=(*allbboxes_)[i*6][entry1];     a_xmax=(*allbboxes_)[i*6+1][entry1];
      a_ymin=(*allbboxes_)[i*6+2][entry1];   a_ymax=(*allbboxes_)[i*6+3][entry1];
      a_zmin=(*allbboxes_)[i*6+4][entry1];   a_zmax=(*allbboxes_)[i*6+5][entry1];

      b_xmin=(*allbboxes_)[j*6][entry2];     b_xmax=(*allbboxes_)[j*6+1][entry2];
      b_ymin=(*allbboxes_)[j*6+2][entry2];   b_ymax=(*allbboxes_)[j*6+3][entry2];
      b_zmin=(*allbboxes_)[j*6+4][entry2];   b_zmax=(*allbboxes_)[j*6+5][entry2];

      // if intersection exists, return true
      if (!((a_xmin >= b_xmax || b_xmin >= a_xmax) || (a_ymin >= b_ymax || b_ymin >= a_ymax) || (a_zmin >= b_zmax || b_zmin >= a_zmax)))
        return true;
    }
  }
  return false;
}

/*-----------------------------------------------------------------------------------*
 |  Cylindrical Oriented Bounding Box Intersection function when both bounding boxes |
 |  represent actual finite elements  (private)                         mueller 11/11|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionCOBB(const std::vector<int>& bboxIDs, Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    dserror("Not in use!");

  /* intersection test by calculating the distance between the two bounding box center lines
   * and comparing it to the respective diameters of the beams*/
  bool intersection = false;


  int bboxid0 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int bboxid1 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  // In case of a hypothetical BB, simply take the last beam element's diameter (does the job for now)
  double bbox0diameter = 0.0;
  double bbox1diameter = 0.0;
  bbox0diameter = (*diameter_)[bboxid0];
  bbox1diameter = (*diameter_)[bboxid1];

  double radiusextrusion = extrusionvalue_->at(1);
  if(additiveextrusion_)
  {
    bbox0diameter += 2.0*radiusextrusion;
    bbox1diameter += 2.0*radiusextrusion;
  }
  else
  {
    bbox1diameter *= radiusextrusion;
    bbox0diameter *= radiusextrusion;
  }

  //Distance at which intersection happens
  double distancelimit = 0.5*(bbox0diameter+bbox1diameter);

  LINALG::TMatrix<double,3,1> t1(true);
  LINALG::TMatrix<double,3,1> t2(true);
  LINALG::TMatrix<double,3,1> r1_a;
  LINALG::TMatrix<double,3,1> r1_b;
  LINALG::TMatrix<double,3,1> r2_a;
  LINALG::TMatrix<double,3,1> r2_b;

  int numshifts_bbox0 = 0;
  int numshifts_bbox1 = 0;

  if(periodicBC_)
  {
    numshifts_bbox0 = (int)(*numshifts_)[bboxid0];
    numshifts_bbox1 = (int)(*numshifts_)[bboxid1];
  }

//  if(numshifts_bbox0>0 || numshifts_bbox1>0)
//    std::cout<<"PAIR "<<bboxid0<<", "<<bboxid1<<"; numshifts: "<<numshifts_bbox0<<", "<<numshifts_bbox1<<std::endl;

//  double t_start = Teuchos::Time::wallTime();
  // intersect
  for(int i=0; i<numshifts_bbox0+1; i++)
  {
    for(int j=0; j<numshifts_bbox1+1; j++)
    {
      for(int k=0; k<(int)r1_a.M(); k++)
      {
        // first bbox
        r1_a(k) = (*allbboxes_)[6*i+k][bboxid0];
        r1_b(k) = (*allbboxes_)[6*i+k+3][bboxid0];
        // second bbox
        r2_a(k) = (*allbboxes_)[6*j+k][bboxid1];
        r2_b(k) = (*allbboxes_)[6*j+k+3][bboxid1];
      }

//      double t_02 = Teuchos::Time::wallTime();
      t1 = FADUTILS::DiffVector(r1_b,r1_a);
      t2 = FADUTILS::DiffVector(r2_b,r2_a);
      double angle = BEAMCONTACT::CalcAngle(t1,t2);
//      std::cout<<"(i,j) = "<<i<<", "<<j<<",  Angle Calc : "<< Teuchos::Time::wallTime()-t_02<<std::endl;

//      double t_03 = Teuchos::Time::wallTime();
      // non-parallel case
      if(angle>ANGLETOL)
      {
        std::pair<double,double> closestpoints(std::make_pair(0.0,0.0));
        bool etaset=false;
        intersection = BEAMCONTACT::IntersectArbitraryCylinders(r1_a,r1_b,r2_a,r2_b,distancelimit,closestpoints,etaset);
//        std::cout<<"(i,j) = "<<i<<", "<<j<<",  Isec nonpar: "<< Teuchos::Time::wallTime()-t_03<<std::endl;
      }
      else // parallel case
      {
        intersection = BEAMCONTACT::IntersectParallelCylinders(r1_a,r1_b,r2_a,r2_b,distancelimit);
//        std::cout<<"(i,j) = "<<i<<", "<<j<<",  Isec par   : "<< Teuchos::Time::wallTime()-t_03<<std::endl;
      }
//      if(numshifts_bbox0>0 || numshifts_bbox1>0)
//        std::cout<<"  (i,j) = "<<i<<", "<<j<<": intersection = "<<intersection<<std::endl;
      if(intersection)
        break;
    }
    // return upon first found intersection, the remainder does not matter anyway+
    if(intersection)
      break;
  }

//  std::cout<<"Intersection total: "<< Teuchos::Time::wallTime()-t_start<<std::endl;
  return intersection;
}

/*-----------------------------------------------------------------------------------*
 |  Spherical Bounding Box Intersection function when both bounding boxes           |
 |  for linkers                                       (private)        mueller 01/12|
 *----------------------------------------------------------------------------------*/
bool Beam3ContactOctTree::IntersectionSPBB(const std::vector<int>& bboxIDs, Teuchos::RCP<Epetra_SerialDenseMatrix> bboxlimits)
{
  // Why bboxlimits separately: The idea is that we can use this method to check whether a hypothetical bounding box (i.e. without an element)
  // can be tested for intersection. Hence, we store the limits of this bounding box into bboxlimits if needed.
  // Since the hypothetical bounding box stands for a crosslinker to be set, we just need the exact dimensions of the element
  if(bboxlimits!=Teuchos::null)
    dserror("Not in use!");

  int bboxid0 = searchdis_.ElementColMap()->LID(bboxIDs[0]);
  int bboxid1 = searchdis_.ElementColMap()->LID(bboxIDs[1]);

  int numshifts_bbox0 = 0;
  int numshifts_bbox1 = 0;

  if(periodicBC_)
  {
    numshifts_bbox0 = (*numshifts_)[bboxid0];
    numshifts_bbox1 = (*numshifts_)[bboxid1];
  }

  for(int i=0; i<numshifts_bbox0+1; i++)
  {
    for(int j=0; j<numshifts_bbox1+1;j++)
    {
      // calculate distance between sphere centers
      double centerdist = sqrt(((*allbboxes_)[0][bboxid0]-(*allbboxes_)[0][bboxid1])*((*allbboxes_)[0][bboxid0]-(*allbboxes_)[0][bboxid1])+
                               ((*allbboxes_)[1][bboxid0]-(*allbboxes_)[1][bboxid1])*((*allbboxes_)[1][bboxid0]-(*allbboxes_)[1][bboxid1])+
                               ((*allbboxes_)[2][bboxid0]-(*allbboxes_)[2][bboxid1])*((*allbboxes_)[2][bboxid0]-(*allbboxes_)[2][bboxid1]));

      double bboxradius0 = 0.5*(*diameter_)[bboxid0];
      double bboxradius1 = 0.5*(*diameter_)[bboxid1];

      // intersect spheres
      if(centerdist <= bboxradius0 + bboxradius1)
        return true; // spheres have overlapping volumes
    }
  }
  return false; // no overlap
}

/*-----------------------------------------------------------------------*
 | communicate Vector to all Processors                    mueller 11/11 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::CommunicateVector(Epetra_Vector& InVec, Epetra_Vector& OutVec, bool zerofy, bool doexport, bool doimport)
{
  /* zerofy InVec at the beginning of each search except for Proc 0
   * for subsequent export and reimport. This way, we guarantee redundant information
   * on all processors. */

  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_.Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::CommunicateMultiVector(Epetra_MultiVector& InVec, Epetra_MultiVector& OutVec, bool zerofy, bool doexport, bool doimport)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec.Map(), InVec.Map());
  Epetra_Import importer(OutVec.Map(), InVec.Map());
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_.Comm().MyPID()!=0 && zerofy)
      OutVec.PutScalar(0.0);
    InVec.Export(OutVec, exporter, Add);
  }
  if(doimport)
    OutVec.Import(InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | Calc. max and min values of node coordinates              meier 05/14 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::CalcCornerPos(DRT::Element* element,
                                        std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
                                        LINALG::SerialDenseMatrix& coord)
{

  LINALG::Matrix<3,1> coord_max(true);
  LINALG::Matrix<3,1> coord_min(true);

  for(int i=0; i<element->NumNode(); i++)
  {
    int gid = element->Nodes()[i]->Id();
    LINALG::Matrix<3,1> coord_aux = currentpositions[gid];
    for(int j=0; j<3; j++)
    {
      if (coord_aux(j)<coord_min(j))
        coord_min(j)=coord_aux(j);
      if (coord_aux(j)>coord_max(j))
        coord_max(j)=coord_aux(j);
    }
  }

  /*  We don't fill real node coordinates in the matrix coord, but only the min/max x-,y-,z- values over
      all nodes of the considered elements. These coordinates coord_max and coord_min can thus be interpreted
      as fictitious nodes which are sufficient in order to create a correct bounding box. */

  for(int j=0; j<3; j++)
  {
    coord(j,0) = coord_min(j);
    coord(j,1) = coord_max(j);
  }

  return;
}

/*-----------------------------------------------------------------------*
 | unshift coordinates               mueller 02/15 |
 *-----------------------------------------------------------------------*/
void Beam3ContactOctTree::UndoEffectOfPeriodicBoundaryCondition(Epetra_SerialDenseMatrix& coord,
                                                                std::vector<int>&         cut,
                                                                int&                      numshifts)
{
  if(coord.M() != 3 || coord.N() != 2)
    dserror("coord must have the dimension M()==3, N()==2!");
  if((int)cut.size()!=3)
    dserror("cut is of wrong size %i!", (int)cut.size());
  // By definition, we shift the second position!
  for (int dof=0; dof<coord.M(); dof++)
  {
    if (fabs(coord(dof,1)-periodlength_->at(dof)-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
    {
      numshifts++;
      cut[dof] = 1.0;
      coord(dof,1) -= periodlength_->at(dof);
    }
    if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
    {
      numshifts++;
      cut[dof] = 2.0;
      coord(dof,1) += periodlength_->at(dof);
    }
  }
  return;
}
/*-----------------------------------------------------------------------*
 | Get bounding-box-specific extrusion value               mueller 02/15 |
 *-----------------------------------------------------------------------*/
double Beam3ContactOctTree::GetBoundingBoxExtrusionValue()
{
  double extrusionvalue = -1.0;
  switch(boundingbox_)
  {
    case Beam3ContactOctTree::axisaligned:
      extrusionvalue = extrusionvalue_->at(0);
      break;
    case Beam3ContactOctTree::cyloriented:
      extrusionvalue = std::max(extrusionvalue_->at(0), extrusionvalue_->at(1));
      break;
    case Beam3ContactOctTree::spherical:
      extrusionvalue = extrusionvalue_->at(0);
      break;
    default:
      dserror("No bounding box type chosen!");
  }
  if(extrusionvalue<0.0)
    dserror("Check bounding box extrusion value %d < 0.0!", extrusionvalue);
  return extrusionvalue;
}
