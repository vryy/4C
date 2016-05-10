/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact.cpp
\brief One beam contact pair (two beam elements) consisting of several contact segments

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3contact.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "beam3contact_tangentsmoothing.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure/strtimint_impl.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_inpar/inpar_statmech.H"

#include "Teuchos_TimeMonitor.hpp"
//TODO: Abfangen, dass Kontaktpunkte am Elementübergang zweimal ausgewertet werden!!!

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
CONTACT::Beam3contact<numnodes, numnodalvalues>::Beam3contact( const DRT::Discretization& pdiscret,
                                                               const DRT::Discretization& cdiscret,
                                                               const std::map<int,int>& dofoffsetmap,
                                                               DRT::Element* element1,
                                                               DRT::Element* element2,
                                                               Teuchos::ParameterList& beamcontactparams):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
bcparams_(beamcontactparams),
iter_(0),
numstep_(0),
R1_(BEAMCONTACT::CalcEleRadius(element1)),
R2_(BEAMCONTACT::CalcEleRadius(element2)),
maxactivegap_(GetMaxActiveDist()),
maxsegdist1_(0.0),
maxsegdist2_(0.0),
numseg1_(1),
numseg2_(1),
boundarynode1_(std::make_pair(0,0)),
boundarynode2_(std::make_pair(0,0))
{
  for (int i=0;i<3*numnodes*numnodalvalues;i++)
  {
    ele1pos_(i)=0.0;
    ele2pos_(i)=0.0;
  }
  for (int i=0;i<3*numnodes;i++)
  {
    nodaltangentssmooth1_(i)=0.0;
    nodaltangentssmooth2_(i)=0.0;
  }

  int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(bcparams_,"BEAMS_SMOOTHING");
  if (smoothing == INPAR::BEAMCONTACT::bsm_cpp)
  {
    const DRT::ElementType & eot1 = element1_->ElementType();
    if(eot1 != DRT::ELEMENTS::Beam3Type::Instance() and eot1 != DRT::ELEMENTS::Beam3rType::Instance())
      dserror("Tangent smoothing only implemented for beams of type beam3 and beam3r!");

    //For both elements the 2 direct neighbor elements are determined and saved in the B3CNeighbor-Class
    //variables neighbors1_ and neighbors2_.
    {
      neighbors1_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element1);
      neighbors2_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element2);
    }
  }

  //In case we want to apply a segment-based integration at the endpoints of the physical beam (in order to avoid
  //strong discontinuities in the integrand) we have to check, if a master beam element node coincides with a beams endpoint!
  bool determine_neighbors=false;
  bool endpointpenalty = DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_ENDPOINTPENALTY");
  if(endpointpenalty) determine_neighbors=true;

  #ifdef ENDPOINTSEGMENTATION
    determine_neighbors=true;
  #endif

  if(determine_neighbors)
  {
    neighbors1_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element1);
    neighbors2_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element2);

    bool leftboundarynode1 = false;
    bool rightboundarynode1 = false;

    if(neighbors1_->GetLeftNeighbor()==NULL)
      leftboundarynode1=true;

    if(neighbors1_->GetRightNeighbor()==NULL)
      rightboundarynode1=true;

    boundarynode1_=std::make_pair(leftboundarynode1,rightboundarynode1);

    bool leftboundarynode2 = false;
    bool rightboundarynode2 = false;

    if(neighbors2_->GetLeftNeighbor()==NULL)
      leftboundarynode2=true;

    if(neighbors2_->GetRightNeighbor()==NULL)
      rightboundarynode2=true;

    boundarynode2_=std::make_pair(leftboundarynode2,rightboundarynode2);
  }

  //TODO: calculate the real length for initially curved elements (e.g. beam3eb_anisotrop)
  //Calculate initial length of beam elements (approximation for initially curved elements!)
  LINALG::TMatrix<double,3,1> lvec1(true);
  LINALG::TMatrix<double,3,1> lvec2(true);
  for(int i=0;i<3;i++)
  {
    lvec1(i)=(element1_->Nodes())[0]->X()[i]-(element1_->Nodes())[1]->X()[i];
    lvec2(i)=(element2_->Nodes())[0]->X()[i]-(element2_->Nodes())[1]->X()[i];
  }

  double l1=lvec1.Norm2();
  double l2=lvec2.Norm2();

  if (element1->ElementType() != element2->ElementType())
    dserror("The class beam3contact only works for contact pairs of the same beam element type!");

  if (element1->Id() >= element2->Id())
    dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");

  int penaltylaw = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(beamcontactparams,"BEAMS_PENALTYLAW");
  if(penaltylaw != INPAR::BEAMCONTACT::pl_lp and penaltylaw != INPAR::BEAMCONTACT::pl_qp)
  {
    if(beamcontactparams.get<double>("BEAMS_PENREGPARAM_F0",-1.0)==-1.0 or beamcontactparams.get<double>("BEAMS_PENREGPARAM_G0",-1.0)==-1.0 or beamcontactparams.get<double>("BEAMS_PENREGPARAM_C0",-1.0)==-1.0)
      dserror("Regularized penalty law chosen, but not all regularization parameters are set!");
  }

  cpvariables_.resize(0);
  gpvariables_.resize(0);
  epvariables_.resize(0);

  if(DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcontactparams,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_uzawa)
    dserror("Uzawa is not implemented for beam3contact elements so far!");

  if (DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_DAMPING")!=INPAR::BEAMCONTACT::bd_no)
    dserror("Damping is not implemented for beam3contact elements so far!");

  if (bcparams_.get<double>("BEAMS_GAPSHIFTPARAM",0.0)!=0.0 and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW")!=INPAR::BEAMCONTACT::pl_lpqp)
    dserror("BEAMS_GAPSHIFTPARAM only possible for penalty law LinPosQuadPen!");

  double perpshiftangle1 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE1")/180.0*M_PI;
  double parshiftangle2 = bcparams_.get<double>("BEAMS_PARSHIFTANGLE2")/180.0*M_PI;

  if(parshiftangle2<=perpshiftangle1)
    dserror("No angle overlap between large-angle and small-angle contact!");

  bool beamsdebug = DRT::INPUT::IntegralValue<int>(beamcontactparams,"BEAMS_DEBUG");

  //Check, if a unique closest point solution can be guaranteed for angles alpha > BEAMS_PERPSHIFTANGLE1
  if((perpshiftangle1<acos(1.0-2*MAXCROSSSECTIONTOCURVATURE)) and !beamsdebug)
    dserror("Choose larger shifting angle BEAMS_PERPSHIFTANGLE1 in order to enable a unique CPP!");

  double segangle = bcparams_.get<double>("BEAMS_SEGANGLE",-1.0)/180.0*M_PI;

  if(bcparams_.get<double>("BEAMS_SEGANGLE",-1.0)<0.0)
    dserror("Input variable BEAMS_SEGANGLE has to be defined!");

  double safetyfac=1.5;
  //Determine bound for search of large-angle contact pairs
  deltalargeangle_=perpshiftangle1-safetyfac*2*segangle;

  //In case of a negative value of deltalargeangle_ all pairs have to be evaluate by the large-angle contact formulation
  if(deltalargeangle_<=0)
    deltalargeangle_=0;

  //Determine bound for search of small-angle contact pairs
  deltasmallangle_ = parshiftangle2+safetyfac*2*segangle;

  //Check, if we have enough gauss points in order to find every contact point!!
  //Calculate maximal length distance between two gauss points (the factor 1.5 takes into account the not evenly distributed locations
  //of the Gauss points -> this does hold for a number of Gauss points <= 10!!!)
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::BEAMCONTACTGAUSSRULE);
  int intintervals = bcparams_.get<int>("BEAMS_NUMINTEGRATIONINTERVAL");

  double deltal1=1.5*l1/(intintervals*gausspoints.nquad);

  if(l2+1.0e-8<l1/intintervals)
    dserror("Length of second (master) beam has to be larger than length of one integration interval on first (slave) beam!");

  if(gausspoints.nquad>10)
    dserror("So far, not more than 10 Gauss points are allowed!");

  //TODO We have not considered the factor of 4 occurring in the formula of maximal Gauss point distance, therefore
  //we have an additional safety factor here...
  //TODO!!!!
  if((deltal1>R1_/sin(parshiftangle2)) and !beamsdebug)
    dserror("Not enough Gauss points crossing of beams possible!!!");

  return;
}
/*----------------------------------------------------------------------*
 |  end: constructor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::Evaluate( LINALG::SparseMatrix& stiffmatrix,
                                                                Epetra_Vector& fint,
                                                                const double& pp,
                                                                std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap,
                                                                Teuchos::ParameterList& timeintparams,
                                                                bool fdcheck)
{
  //**********************************************************************
  // Evaluation of contact forces and stiffness
  //**********************************************************************
  // (1) Closest Point Projection (CPP)
  //     -> find closest point where contact forces are evaluated
  // (2) Compute some auxiliary quantities
  //     -> normal vector, gap, shape functions, contact flag,
  //     -> linearizations of all geometric quantities
  // (3) Compute contact forces and stiffness
  //     -> stiffness terms are directly assembled to global matrix
  //     -> contact forces are only returned as global vector
  // (4) Perform some finite difference checks
  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
  //***************Get some parameters in the beginning*******************

#ifdef FDCHECK
  if(fdcheck==false)
    FDCheck(stiffmatrix,fint,pp,contactpairmap,timeintparams,fdcheck);
#endif

  //All updates that have to be done in every iteration have do be done here,
  //since most of the elements leave directly after the closest point projection!
  SetClassVariables(timeintparams);

  //Subdevide the two elements in segments with linear approximation
  std::vector<LINALG::TMatrix<double,3,1> > endpoints1(0);
  std::vector<LINALG::TMatrix<double,3,1> > endpoints2(0);

  //TODO: remove 0 and 1: So far the number 0 and 1 are used in order to distinguish
  //between element 1 and element 2. However, this is only necessary for debugging purposes
  //and can be removed later!
  maxsegdist1_=CreateSegments(element1_, endpoints1, numseg1_,0);
  maxsegdist2_=CreateSegments(element2_, endpoints2, numseg2_,1);

  //Make pairs of close segments: Most of the pairs are already sorted out
  //at this point and don't have to be considered further in the following CPP
  //Additionally, we store the relative orientation of the pairs

  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> > closelargeanglesegments;
  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> > closesmallanglesegments;
  std::vector<std::pair<int,int> > closeendpointsegments(0);
  closelargeanglesegments.clear();
  closesmallanglesegments.clear();

  bool endpoint_penalty=DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_ENDPOINTPENALTY");

  //Sub-devision of contact elements in search segments or not?
  #ifndef NOSEGMENTATION
    GetCloseSegments(endpoints1,endpoints2,closesmallanglesegments,closelargeanglesegments,closeendpointsegments,maxactivegap_);
  #else
    LINALG::TMatrix<double,3,1> segmentdata(true);
    segmentdata(0)=0.0;  //segment angle
    segmentdata(1)=0.0;    //eta1_seg
    segmentdata(2)=0.0;    //eta2_seg
    closesmallanglesegments[std::make_pair(0,0)]=segmentdata;
    closelargeanglesegments[std::make_pair(0,0)]=segmentdata;

    if(endpoint_penalty)
    {
      if(boundarynode1_.first or boundarynode1_.second or boundarynode2_.first or boundarynode2_.second)
        closeendpointsegments.push_back(std::make_pair(0,0));
    }
  #endif

  //**********************************************************************
  // (1) Closest Point Projection for all close large angle segments(CPP)
  //**********************************************************************

  //Treat large-angle contact pairs if existing
  if(closelargeanglesegments.size()>0)
  {
    //Get active large angle pairs (valid closest point projections) and create vector of cpvariables_
    GetActiveLargeAnglePairs(endpoints1,endpoints2,closelargeanglesegments,pp);

    //Evaluate contact contribution of large-angle-contact (residual and stiffness) for all closest points found before
    EvaluateActiveLargeAnglePairs(stiffmatrix,fint);
  }

  //Treat small angle contact pairs if existing
  if(closesmallanglesegments.size()>0)
  {
    #ifndef ENDPOINTSEGMENTATION
      //Get active small angle pairs (valid Gauss points) and create vector of gpvariables_
      GetActiveSmallAnglePairs(closesmallanglesegments);

      //Evaluate contact contribution of small-angle-contact (residual and stiffness) for all closest points found before
      EvaluateActiveSmallAnglePairs(stiffmatrix,fint);
    #else
      //In case of endpoint segmentation some additional quantities have to be transferred between the methods
      //GetActiveSmallAnglePairs() and EvaluateActiveSmallAnglePairs().
      std::pair<int,int> iminmax = std::make_pair(0,0);
      std::pair<bool,bool> leftrightsolutionwithinsegment = std::make_pair(false,false);
      std::pair<double,double> eta1_leftrightboundary = std::make_pair(0.0, 0.0);

      //Get active small angle pairs (valid Gauss points) and create vector of gpvariables_
      GetActiveSmallAnglePairs(closesmallanglesegments,&iminmax,&leftrightsolutionwithinsegment,&eta1_leftrightboundary);

      //Evaluate contact contribution of small-angle-contact (residual and stiffness) for all closest points found before
      EvaluateActiveSmallAnglePairs(stiffmatrix,fint,&iminmax,&leftrightsolutionwithinsegment,&eta1_leftrightboundary);
    #endif
  }

  if(endpoint_penalty)
  {
    //Treat endpoint contact pairs if existing
    if(closeendpointsegments.size()>0)
    {
      //Get active endpoint pairs and create vector of epvariables_
      GetActiveEndPointPairs(closeendpointsegments,pp);

      //Evaluate contact contribution of endpoint-contact (residual and stiffness) for all closest points found before
      EvaluateActiveEndPointPairs(stiffmatrix,fint);
    }
  }

  return (true);
}
/*----------------------------------------------------------------------*
 |  end: Evaluate the element
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active large angle pairs                             meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetActiveLargeAnglePairs(std::vector<LINALG::TMatrix<double,3,1> >& endpoints1,
                                                                               std::vector<LINALG::TMatrix<double,3,1> >& endpoints2,
                                                                               std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closelargeanglesegments,
                                                                               const double pp)
{
  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >::iterator iter;

    for (iter=closelargeanglesegments.begin(); iter != closelargeanglesegments.end(); ++iter)
    {
      LINALG::TMatrix<double,3,1> segmentdata = iter->second;
      std::pair<int,int> leftpoint_ids = iter->first;
      int nseg1 = endpoints1.size()-1;
      int nseg2 = endpoints2.size()-1;
      int segid1 = leftpoint_ids.first;
      int segid2 = leftpoint_ids.second;
      double l1 = 2.0/nseg1;
      double l2 = 2.0/nseg2;
      double eta_left1 = -1.0+segid1*l1;
      double eta_left2 = -1.0+segid2*l2;

      std::pair<TYPE,TYPE> closestpoint(0.0,0.0);
      bool validpairfound=false;

      //The method ClosestPointProjection() only delivers a valid solution (validpairfound=true), if eta1 \in [eta_left1,eta_left1+l1],
      //eta2 \in [eta_left2,eta_left2+l2] and gap<maxactivegap_!
      validpairfound=ClosestPointProjection(eta_left1,eta_left2,l1,l2,segmentdata,closestpoint,segid1,segid2);

      //With the following block we sort out identical contact points that occure more than once
      //within this element pair -> this is possible, when the contact point lies on the boundary
      //between two segments!
      bool allready_found=false;

      for(int i=0;i<(int)cpvariables_.size();i++)
      {
        double eta1_eval = FADUTILS::CastToDouble(cpvariables_[i]->GetCP().first);
        double eta2_eval = FADUTILS::CastToDouble(cpvariables_[i]->GetCP().second);

        if(fabs(eta1_eval-FADUTILS::CastToDouble(closestpoint.first))<XIETARESOLUTIONFAC*XIETAITERATIVEDISPTOL and fabs(eta2_eval-FADUTILS::CastToDouble(closestpoint.second))<XIETARESOLUTIONFAC*XIETAITERATIVEDISPTOL)
          allready_found=true;
      }

      if(validpairfound and !allready_found)
      {
        std::pair<int,int> integration_ids = std::make_pair(-2,-2);
        cpvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }
}
/*----------------------------------------------------------------------*
 |  end: Get active large angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active large angle pairs                        meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateActiveLargeAnglePairs(LINALG::SparseMatrix& stiffmatrix,
                                                                                    Epetra_Vector& fint)
{
  for (int numcp=0;numcp<(int)cpvariables_.size();numcp++)
  {
    //**********************************************************************
    // (2) Compute some auxiliary quantities
    //**********************************************************************

    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    TYPE eta1=cpvariables_[numcp]->GetCP().first;
    TYPE eta2=cpvariables_[numcp]->GetCP().second;

    #ifdef AUTOMATICDIFF
      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
    #endif

    // update shape functions and their derivatives
    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    ComputeNormal(r1, r2, r1_xi, r2_xi, cpvariables_[numcp],0);

    // call function to compute penalty force
    CalcPenaltyLaw(cpvariables_[numcp]);

    // get shift angles from input file
    double perpshiftangle1 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE1")/180.0*M_PI;
    double perpshiftangle2 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE2")/180.0*M_PI;

    // call function to compute scale factor of penalty parameter
    CalcPerpPenaltyScaleFac(cpvariables_[numcp],r1_xi,r2_xi,perpshiftangle1,perpshiftangle2);

    //In case of large-angle-contact, the length specific energy and the 'real' energy are identical
    double lengthspec_energy=FADUTILS::CastToDouble(cpvariables_[numcp]->GetEnergy());
    cpvariables_[numcp]->SetIntegratedEnergy(lengthspec_energy);

//    std::cout << "cpvariables_[numcp]->GetNormal(): " << cpvariables_[numcp]->GetNormal() << std::endl;
//    std::cout << "numcp: " << numcp << std::endl;
//    std::cout << "xi: " << cpvariables_[numcp]->GetCP().first.val() << std::endl;
//    std::cout << "eta: " << cpvariables_[numcp]->GetCP().second.val() << std::endl;
//    std::cout << "gap: " << cpvariables_[numcp]->GetGap().val() << std::endl;
//    std::cout << "angle: " << cpvariables_[numcp]->GetAngle()/M_PI*180.0 << std::endl;
//    std::cout << "r1_xi: " << r1_xi << std::endl;
//    std::cout << "r2_xi: " << r2_xi << std::endl;
//    std::cout << "|r1_xi|: " << r1_xi.Norm2() << std::endl;
//    std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
//    std::cout << "r1_xi*r2_xi: " << FADUTILS::ScalarProduct(r1_xi,r2_xi) << std::endl;
//    std::cout << "cpvariables_[numcp]->Getfp(): " << cpvariables_[numcp]->Getfp() << std::endl;

    // call function to compute contact contribution to residual vector
    EvaluateFcContact(&fint, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, cpvariables_[numcp], 1.0,true,false,false,false);

    // call function to compute contact contribution to stiffness matrix
    EvaluateStiffcContact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi,cpvariables_[numcp],1.0,true,false,false,false);

  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active large angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active small angle pairs                             meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetActiveSmallAnglePairs( std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closesmallanglesegments,
                                                                                std::pair<int,int>* iminmax,
                                                                                std::pair<bool,bool>* leftrightsolutionwithinsegment,
                                                                                std::pair<double,double>* eta1_leftrightboundary)
{
  //lengths in parameter space of created segments
  double l1=2.0/numseg1_;
  double l2=2.0/numseg2_;

  int numpairs = closesmallanglesegments.size();
  std::vector<std::pair<double, double > > inversepairs(numpairs,std::make_pair(0.0,0.0));
  int pairiter=0;
  int intintervals = bcparams_.get<int>("BEAMS_NUMINTEGRATIONINTERVAL");

  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >::iterator iter;

  for (iter=closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
  {
    std::pair<int,int> segment_ids = iter->first;
    int segid1 = segment_ids.first;
    int segid2 = segment_ids.second;
    double eta1_segleft = -1.0+segid1*l1;
    double eta2_segleft = -1.0+segid2*l2;

    inversepairs[numpairs-pairiter-1]=std::make_pair(eta1_segleft,eta2_segleft);
    pairiter++;
  }
  int imin = 0;
  int imax = intintervals -1;

  #ifdef ENDPOINTSEGMENTATION
  double eta1_leftboundary = -1.0;
  double eta1_rightboundary = 1.0;
  double leftintervallength = 2.0/intintervals;
  double rightintervallength = 2.0/intintervals;
  bool leftsolutionwithinsegment=false;
  bool rightsolutionwithinsegment=false;

  if(boundarynode2_.first)
  {
    for (iter=closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
    {
      std::pair<int,int> segment_ids = iter->first;
      int segid2 = segment_ids.second;

      if(segid2==0)
      {
        int segid1 = segment_ids.first;
        double eta1_segleft = -1.0+segid1*l1;
        double eta2_segleft = -1.0;
        double eta1_boundary_trial = 0.0;
        bool dummy;

        bool solutionwithinsegment=false;

        double gap_dummy=0.0;
        double alpha_dummy=0.0;

        #ifndef CHANGEENDPOINTPROJECTION
          solutionwithinsegment=PointToLineProjection(eta2_segleft,eta1_segleft,l1,eta1_boundary_trial,gap_dummy,alpha_dummy,dummy,true,true);
        #else
          solutionwithinsegment=PointToLineProjection(eta2_segleft,eta1_segleft,l1,eta1_boundary_trial,gap_dummy,alpha_dummy,dummy,true,true,true);
        #endif

        if(solutionwithinsegment)
        {
          //Determine if the projection eta1_boundary_trial is a left boundary of the integration segment or a right boundary of the integration segment
          //This is done in the following way: First, we determine the tangent of the master boundary node in a way, such that the tangent points into the
          //elements interior. Then, we determine the tangent on slave beam at the projection point eta1_boundary_trial. This tangent automatically points
          //into positive eta1-direction=integration direction. Thus, if the scalar product of these two tangents is positive, the master element evolves
          //in positive eta1-direction and consequently, eta1_boundary_trial is the left boundary of the integration segment. If the scalar product is negative,
          //eta1_boundary_trial is the right boundary of the integration segment
          LINALG::TMatrix<TYPE, 3, 1> inward_tangent_master=r_xi(eta2_segleft, element2_);
          LINALG::TMatrix<TYPE, 3, 1> tangent_slave=r_xi(eta1_boundary_trial, element1_);
          double orientation=FADUTILS::CastToDouble(FADUTILS::ScalarProduct(inward_tangent_master,tangent_slave));
          if(orientation>0)//left boundary
          {
            leftsolutionwithinsegment=true;
            eta1_leftboundary=eta1_boundary_trial;
            //determine ID of integration interval in which the point eta1_leftboundary lies
            imin = BEAMCONTACT::GetIntervalId(eta1_leftboundary,intintervals,true);
            //get length of segmented integration interval
            leftintervallength=-1.0+(imin+1)*2.0/intintervals-eta1_leftboundary;
            break;
          }
          else if(orientation<0)//right boundary
          {
            rightsolutionwithinsegment=true;
            eta1_rightboundary=eta1_boundary_trial;
            //determine ID of integration interval in which the point eta1_leftboundary lies
            imax = BEAMCONTACT::GetIntervalId(eta1_rightboundary,intintervals,false);
            //get length of segmented integration interval
            rightintervallength=eta1_rightboundary-(-1.0+imax*2.0/intintervals);
            break;
          }
          else //This can only happen, if both beams are exactly perpendicular AND the master beam endpoint projects perpendicular on the slave beam!
            dserror("The very unlikely case orientation==0 is not implemented so far!");
        }
      }
    }
  }
  if(boundarynode2_.second)
  {
    for (iter=closesmallanglesegments.begin(); iter != closesmallanglesegments.end(); ++iter)
    {
      std::pair<int,int> segment_ids = iter->first;
      int segid2 = segment_ids.second;

      if(segid2==numseg2_-1)
      {
        int segid1 = segment_ids.first;
        double eta1_segleft = -1.0+segid1*l1;
        double eta2_segright = 1.0;
        double eta1_boundary_trial = 0.0;
        bool dummy;

        bool solutionwithinsegment=false;

        double gap_dummy=0.0;
        double alpha_dummy=0.0;

        #ifndef CHANGEENDPOINTPROJECTION
          solutionwithinsegment=PointToLineProjection(eta2_segright,eta1_segleft,l1,eta1_boundary_trial,gap_dummy,alpha_dummy,dummy,true,true);
        #else
          solutionwithinsegment=PointToLineProjection(eta2_segright,eta1_segleft,l1,eta1_boundary_trial,gap_dummy,alpha_dummy,dummy,true,true,true);
        #endif

        if(solutionwithinsegment)
        {
          //Determine if the projection eta1_boundary_trial is a left boundary of the integration segment or a right boundary of the integration segment
          //This is done in the following way: First, we determine the tangent of the master boundary node in a way, such that the tangent points into the
          //elements interior. Then, we determine the tangent on slave beam at the projection point eta1_boundary_trial. This tangent automatically points
          //into positive eta1-direction=integration direction. Thus, if the scalar product of these two tangents is positive, the master element evolves
          //in positive eta1-direction and consequently, eta1_boundary_trial is the left boundary of the integration segment. If the scalar product is negative,
          //eta1_boundary_trial is the right boundary of the integration segment
          LINALG::TMatrix<TYPE, 3, 1> inward_tangent_master=r_xi(eta2_segright, element2_);
          //Scale tangent of right element node (eta2=1.0) in order to get inward tangent!
          inward_tangent_master.Scale(-1.0);
          LINALG::TMatrix<TYPE, 3, 1> tangent_slave=r_xi(eta1_boundary_trial, element1_);
          double orientation=FADUTILS::CastToDouble(FADUTILS::ScalarProduct(inward_tangent_master,tangent_slave));
          if(orientation>0)//left boundary
          {
            if(leftsolutionwithinsegment)
              dserror("Something went wrong here: both boundary nodes of the master beam (discretized by one finite element?!?) are prejected as left boundary of the integration segment!");

            leftsolutionwithinsegment=true;
            eta1_leftboundary=eta1_boundary_trial;
            //determine ID of integration interval in which the point eta1_leftboundary lies
            imin = BEAMCONTACT::GetIntervalId(eta1_leftboundary,intintervals,true);
            //get length of segmented integration interval
            leftintervallength=-1.0+(imin+1)*2.0/intintervals-eta1_leftboundary;
            break;
          }
          else if(orientation<0)//right boundary
          {
            if(rightsolutionwithinsegment)
              dserror("Something went wrong here: both boundary nodes of the master beam (discretized by one finite element?!?) are prejected as right boundary of the integration segment!");

            rightsolutionwithinsegment=true;
            eta1_rightboundary=eta1_boundary_trial;
            //determine ID of integration interval in which the point eta1_leftboundary lies
            imax = BEAMCONTACT::GetIntervalId(eta1_rightboundary,intintervals,false);
            //get length of segmented integration interval
            rightintervallength=eta1_rightboundary-(-1.0+imax*2.0/intintervals);
            break;
          }
          else //This can only happen, if both beams are exactly perpendicular AND the master beam endpoint projects perpendicular on the slave beam!
            dserror("The very unlikely case orientation==0 is not implemented so far!");
        }
      }
    }
  }
  if(leftsolutionwithinsegment and rightsolutionwithinsegment and imin==imax)
    dserror("It is not possible to cut an integration interval from both sides, choose a larger value intintervals!");
  #endif

  // gaussian points
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::BEAMCONTACTGAUSSRULE);

  //loop over all integration intervals
  for(int interval=imin;interval<=imax;interval++)
  {
    //Calculate parameter bounds of considered intergration interval
    double eta1_min = -1.0+interval*2.0/intintervals;
    double eta1_max = -1.0+(interval+1)*2.0/intintervals;

    //Get jacobi factor of considered interval
    TYPE jacobi_interval=1.0;

    //standard case of equidistant intervals
    #ifndef ENDPOINTSEGMENTATION
      //map from segment coordinate xi to element coordinate eta
    jacobi_interval=1.0/intintervals;
    //case of smaller integration intervals due to segmentation at the beams endpoints
    #else
      if(interval==imin and leftsolutionwithinsegment)
      {
        jacobi_interval= leftintervallength/2.0;
      }
      else if(interval==imax and rightsolutionwithinsegment)
      {
        jacobi_interval= rightintervallength/2.0;
        //std::cout << "rightintervallength: " << rightintervallength << std::endl;
      }
      else
      {
        jacobi_interval=1.0/intintervals;
      }
    #endif

    std::vector<std::pair<double, double > > curintsegpairs(0);
    int size = inversepairs.size();

    //All segment pairs for which the segment on the slave beam 1 intersects with the considered integration interval are filtered
    //out and stored in the vector curintsegpairs. These pairs are relevant for the integration procedure on the current interval.
    for(int k=0;k<size;k++)
    {
      double eta1_segleft = (inversepairs[size-1-k]).first;
      double eta1_segright = eta1_segleft + l1;
      //Since the vector inversepairs is sorted with respect to the location of the slave segment (the slave segment with
      //the lowest bounding parameter coordinates eta1_segleft and eta1_segright lie on the last position of the vector
      //inversepairs), it is sufficient to start with the last element and leave the k-loop as soon as we have found the first
      //segment pair without intersection. This procedure only works, if we delete a segment pair as soon as we realize that it will
      //not be relevant for the next integration interval anymore, see comment at (*).
      if(eta1_segleft<eta1_max+1.0e-10)
      {
        //store relevant pairs in new vector
        curintsegpairs.push_back(inversepairs[size-1-k]);

        //(*) In case eta1_segright (the largest parameter coordinate lying within the slave segment) is smaller than
        //eta1_max(the upper bound of the integration interval), the considered segment will not be relevant for the next
        //integration interval at i+1 and can be deleted.
        if(eta1_segright< eta1_max-1.0e-10)
        {
          inversepairs.pop_back();
        }
      }
      //In case we have no relevant segment pair, we will leave the loop already after the first iteration!
      else
      {
        break;
      }
    }

    //If segments exist, evaluate the corresponding Gauss points
    if(curintsegpairs.size()>0)
    {
      //loop over Gauss point of considered integration interval
      for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
      {
        //integration points in parameter space and weights
        const double xi = gausspoints.qxg[numgp][0];

        //Get Gauss point coordinate at slave element
        double eta1_slave=0.0;

        //standard case of equidistant intervals
        #ifndef ENDPOINTSEGMENTATION
          //map from segment coordinate xi to element coordinate eta
          eta1_slave=eta1_min + (1.0+xi)/intintervals;
        //case of smaller integration intervals due to segmentation at the beams endpoints
        #else
          if(interval==imin and leftsolutionwithinsegment)
          {
            eta1_slave=eta1_leftboundary +(1.0+xi)/2.0*leftintervallength;
          }
          else if(interval==imax and rightsolutionwithinsegment)
          {
            eta1_slave=eta1_min +(1.0+xi)/2.0*rightintervallength;
          }
          else
          {
            eta1_slave=eta1_min + (1.0+xi)/intintervals;
          }
        #endif

        for(int k=0;k<(int)curintsegpairs.size();k++)
        {
          double eta1_segleft = (curintsegpairs[k]).first;
          double eta1_segright = eta1_segleft + l1;

          //TODO: This procedure can also be made more efficient by deleting all segments of curintsegpairs which
          //are not relevant for the following Gauss points anymore (see intersection of integration intervals and segment pairs)
          if(BEAMCONTACT::WithinInterval(eta1_slave,eta1_segleft,eta1_segright))
          {
            double eta2_segleft = (curintsegpairs[k]).second;
            double eta2_master = 0.0;
            bool pairactive=false;

            double gap_dummy=0.0;
            double alpha_dummy=0.0;

            bool solutionwithinsegment=PointToLineProjection(eta1_slave,eta2_segleft,l2,eta2_master,gap_dummy,alpha_dummy,pairactive,true);

            if(solutionwithinsegment)
            {

              if(pairactive)
              {
                TYPE eta1=eta1_slave;
                TYPE eta2=eta2_master;
                int leftpoint_id1=BEAMCONTACT::GetSegmentId(eta1_slave,numseg1_);
                int leftpoint_id2=BEAMCONTACT::GetSegmentId(eta2_master,numseg2_);
                std::pair<TYPE,TYPE> closestpoint(std::make_pair(eta1,eta2));
                std::pair<int,int> integration_ids = std::make_pair(numgp,interval);
                std::pair<int,int> leftpoint_ids = std::make_pair(leftpoint_id1,leftpoint_id2);
                TYPE jacobi = GetJacobi(element1_)*jacobi_interval;

                const double parallel_pp = bcparams_.get<double>("BEAMS_BTBLINEPENALTYPARAM");

                if(parallel_pp<0.0)
                  dserror("BEAMS_BTBLINEPENALTYPARAM not set!");

                //Create data container for each Gauss point (in case of small-angle contact the number of the Gauss point [numgp] and the number of the
                //integration interval [interval] are stored in the pair segids_ of the class beamcontactvaribles!)
                gpvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,parallel_pp,jacobi)));
              }
              //We can leave the k-loop as soon as we have found a valid projection for the given Gauss point eta1_slave
              break;
            }
          }
        }
      }//for (int numgp=0; numgp<gausspoints.nquad; ++numgp)
    }//if(curintsegpairs.size()>0)
  }//for(int interval=imin;interval<=imax;interval++)

#ifdef ENDPOINTSEGMENTATION
  if(iminmax==NULL or leftrightsolutionwithinsegment==NULL or eta1_leftrightboundary==NULL)
    dserror("In case of ENDPOINTSEGMENTATION no NUll pointer should be handeld in!!!");

  *iminmax = std::make_pair(imin,imax);
  *leftrightsolutionwithinsegment = std::make_pair(leftsolutionwithinsegment,rightsolutionwithinsegment);
  *eta1_leftrightboundary = std::make_pair(eta1_leftboundary, eta1_rightboundary);
#endif

}
/*----------------------------------------------------------------------*
 |  end: Get active small angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active small angle pairs                        meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateActiveSmallAnglePairs(LINALG::SparseMatrix& stiffmatrix,
                                                                                    Epetra_Vector& fint,
                                                                                    std::pair<int,int>* iminmax,
                                                                                    std::pair<bool,bool>* leftrightsolutionwithinsegment,
                                                                                    std::pair<double,double>* eta1_leftrightboundary)
{
  //Compute linearizations of integration interval boundaries if necessary
  #ifdef ENDPOINTSEGMENTATION

  if(iminmax==NULL or leftrightsolutionwithinsegment==NULL or eta1_leftrightboundary==NULL)
    dserror("In case of ENDPOINTSEGMENTATION no NUll pointer should be handeld in!!!");

  int imin=(*iminmax).first;
  int imax=(*iminmax).second;

  bool leftsolutionwithinsegment=(*leftrightsolutionwithinsegment).first;
  bool rightsolutionwithinsegment=(*leftrightsolutionwithinsegment).second;

  double eta1_leftboundary =(*eta1_leftrightboundary).first;
  double eta1_rightboundary =(*eta1_leftrightboundary).second;

  LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1> delta_xi_R(true);
  LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1> delta_xi_L(true);

  if(leftsolutionwithinsegment)
  {
    TYPE eta1_bound=eta1_leftboundary;
    TYPE eta2=-1.0;
    ComputeLinXiBound(delta_xi_L,eta1_bound,eta2);
  }
  if(rightsolutionwithinsegment)
  {
    TYPE eta1_bound=eta1_rightboundary;
    TYPE eta2=1.0;
    ComputeLinXiBound(delta_xi_R,eta1_bound,eta2);
  }
  #endif

  // gaussian points
  DRT::UTILS::IntegrationPoints1D gausspoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::BEAMCONTACTGAUSSRULE);

  //Evaluate all active Gauss points
  for (int numgptot=0;numgptot<(int)gpvariables_.size();numgptot++)
  {
    std::pair<TYPE,TYPE> closestpoint = gpvariables_[numgptot]->GetCP();
    TYPE eta1=closestpoint.first;
    TYPE eta2=closestpoint.second;

    #ifdef AUTOMATICDIFF
      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1,eta2);
    #endif

    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    // update shape functions and their derivatives
    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    ComputeNormal(r1, r2, r1_xi, r2_xi, gpvariables_[numgptot],1);

    // call function to compute penalty force
    CalcPenaltyLaw(gpvariables_[numgptot]);

    // get shift angles from input file
    double parshiftangle1 = bcparams_.get<double>("BEAMS_PARSHIFTANGLE1")/180.0*M_PI;
    double parshiftangle2 = bcparams_.get<double>("BEAMS_PARSHIFTANGLE2")/180.0*M_PI;

    // call function to compute scale factor of penalty parameter
    CalcParPenaltyScaleFac(gpvariables_[numgptot],r1_xi,r2_xi,parshiftangle1,parshiftangle2);

//    std::cout << "gpvariables_[numgp]->GetNormal(): " << gpvariables_[numgp]->GetNormal() << std::endl;
//    std::cout << "numgptot: " << numgptot << std::endl;
//    std::cout << "xi: " << gpvariables_[numgptot]->GetCP().first.val() << std::endl;
//    std::cout << "eta: " << gpvariables_[numgptot]->GetCP().second.val() << std::endl;
//    std::cout << "gap: " << gpvariables_[numgptot]->GetGap().val() << std::endl;
//    std::cout << "angle: " << gpvariables_[numgptot]->GetAngle()/M_PI*180.0 << std::endl;
//    std::cout << "r1_xi: " << r1_xi << std::endl;
//    std::cout << "r2_xi: " << r2_xi << std::endl;
//    std::cout << "|r1_xi|: " << r1_xi.Norm2() << std::endl;
//    std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
//    std::cout << "r1_xi*r2_xi: " << FADUTILS::ScalarProduct(r1_xi,r2_xi) << std::endl;
//    std::cout << "gpvariables_[numgptot]->Getfp(): " << gpvariables_[numgp]->Getfp() << std::endl;

    //Determine the integration-segment-local Gauss point-ID of the considered gpvariable
    int numgploc = gpvariables_[numgptot]->GetIntIds().first;

    double weight = gausspoints.qwgt[numgploc];
    TYPE jacobi = gpvariables_[numgptot]->GetJacobi();

    #ifdef ENDPOINTSEGMENTATION
      int numinterval = gpvariables_[numgptot]->GetIntIds().second;
    #endif

    //TODO: Here we apply an element jacobian that is constant along the beam element. This works only for initially straight elements!
    //Furthermore we assume, that the element is subdivided in a total of intintervals integration intervals of equal length!
    //The intfac has NOT to be of TYPE FAD in order to deal with non-constant jacobis (in case of ENDPOINTSEGMENTATION) since we explicitly
    //consider the linearization of the jacobi in EvaluateStiffcContactIntSeg()!
    double intfac = FADUTILS::CastToDouble(jacobi)*weight;

    //Convert the length specific energy into a 'real' energy
    //while the length specific energy is used for later calculation, the real (or over the length integrated) energy
    //is a pure output variable and can therefore be of type double!
    double lengthspec_energy=FADUTILS::CastToDouble(gpvariables_[numgptot]->GetEnergy());
    double integrated_energy=lengthspec_energy*intfac;
    gpvariables_[numgptot]->SetIntegratedEnergy(integrated_energy);

    // call function to compute contact contribution to residual vector
    EvaluateFcContact(&fint, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, gpvariables_[numgptot],intfac,false,true,false,false);

    #ifndef ENDPOINTSEGMENTATION
    // call function to compute contact contribution to stiffness matrix
    EvaluateStiffcContact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi,gpvariables_[numgptot],intfac,false,true,false,false);
    #else
    TYPE jacobi_interval = jacobi/GetJacobi(element1_);
    // In case of segment-based integration, we apply a special FAD linearization technique
    // Case 1: segmentation on left side of integration interval
    if(leftsolutionwithinsegment and numinterval==imin)
    {
      //We need the linearization of the mapping from the element parameter space to the integration interval parameter
      //space: xi_ele=xi_left*(1.0-xi_local)/2.0+xi_right*(1.0+xi_local)/2.0.
      //-> d(xi_ele)/d(xi_left)=(1.0-xi_local)/2.0 and d(xi_ele)/d(xi_right)=(1.0+xi_local)/2.0

      double d_xi_ele_d_xi_left=(1.0-gausspoints.qxg[numgploc][0])/2.0;
      EvaluateStiffcContactIntSeg(stiffmatrix, delta_xi_L, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, gpvariables_[numgptot],intfac,d_xi_ele_d_xi_left,-jacobi_interval);
    }
    // Case 2: segmentation on right side of integration interval
    else if(rightsolutionwithinsegment and numinterval==imax)
    {
      double d_xi_ele_d_xi_right=(1.0+gausspoints.qxg[numgploc][0])/2.0;
      EvaluateStiffcContactIntSeg(stiffmatrix, delta_xi_R, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, gpvariables_[numgptot],intfac,d_xi_ele_d_xi_right,jacobi_interval);
    }
    // Case 3: No segmentation necessary
    else
    {
      EvaluateStiffcContact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi,gpvariables_[numgptot],intfac,false,true,false,false);
    }
    #endif
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active small angle pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get active endpoint pairs                                meier 12/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetActiveEndPointPairs(std::vector<std::pair<int,int> >& closeendpointsegments, const double pp)
{
  for(int i=0;i<(int)closeendpointsegments.size();i++)
  {
    int segid1=closeendpointsegments[i].first;
    int segid2=closeendpointsegments[i].second;

    //lengths in parameter space of created segments
    double l1=2.0/numseg1_;
    double l2=2.0/numseg2_;
    //element parameter coordinates of left segment boundary
    double eta1_segleft = -1.0+segid1*l1;
    double eta2_segleft = -1.0+segid2*l2;

    if(segid1==0 and boundarynode1_.first)
    {
      //given parameter coordinate
      double eta1=-1.0;
      //searched parameter coordinate
      double eta2=0.0;
      bool pairactive=false;

      double gap_dummy=0.0;
      double alpha_dummy=0.0;

      bool solutionwithinsegment=PointToLineProjection(eta1,eta2_segleft,l2,eta2,gap_dummy,alpha_dummy,pairactive,false);

      if(solutionwithinsegment and pairactive)
      {
          std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

          std::pair<int,int> integration_ids = std::make_pair(1,0);
          std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

          //Create data container for each end point
          //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
          //has been considered as endpoint, b has the same meaning for element 2
          epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }

    #ifndef ONLYLEFTENDPOINTCONTACT
    if(segid1==numseg1_-1 and boundarynode1_.second)
    {
      //given parameter coordinate
      double eta1=1.0;
      //searched parameter coordinate
      double eta2=0.0;
      bool pairactive=false;

      double gap_dummy=0.0;
      double alpha_dummy=0.0;

      bool solutionwithinsegment=PointToLineProjection(eta1,eta2_segleft,l2,eta2,gap_dummy,alpha_dummy,pairactive,false);

      if(solutionwithinsegment and pairactive)
      {
          std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

          std::pair<int,int> integration_ids = std::make_pair(1,0);
          std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

          //Create data container for each end point
          //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
          //has been considered as endpoint, b has the same meaning for element 2
          epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }
    #endif

    if(segid2==0 and boundarynode2_.first)
    {
      //given parameter coordinate
      double eta2=-1.0;
      //searched parameter coordinate
      double eta1=0.0;
      bool pairactive=false;

      double gap_dummy=0.0;
      double alpha_dummy=0.0;

      bool solutionwithinsegment=PointToLineProjection(eta2,eta1_segleft,l1,eta1,gap_dummy,alpha_dummy,pairactive,false,true);

      if(solutionwithinsegment and pairactive)
      {
          std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

          std::pair<int,int> integration_ids = std::make_pair(0,1);
          std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

          //Create data container for each end point
          //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
          //has been considered as endpoint, b has the same meaning for element 2
          epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }

    #ifndef ONLYLEFTENDPOINTCONTACT
    if(segid2==numseg2_-1 and boundarynode2_.second)
    {
      //given parameter coordinate
      double eta2=1.0;
      //searched parameter coordinate
      double eta1=0.0;
      bool pairactive=false;

      double gap_dummy=0.0;
      double alpha_dummy=0.0;

      bool solutionwithinsegment=PointToLineProjection(eta2,eta1_segleft,l1,eta1,gap_dummy,alpha_dummy,pairactive,false,true);

      if(solutionwithinsegment and pairactive)
      {
          std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

          std::pair<int,int> integration_ids = std::make_pair(0,1);
          std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

          //Create data container for each end point
          //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
          //has been considered as endpoint, b has the same meaning for element 2
          epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }
    #endif

    if((segid1==0 and boundarynode1_.first) and (segid2==0 and boundarynode2_.first))
    {
      double eta1=-1.0;
      double eta2=-1.0;
      LINALG::Matrix<3,1> deltanodalpos(true);
      for(int i=0;i<3;i++)
      {
        deltanodalpos(i)=FADUTILS::CastToDouble(ele2pos_(i)-ele1pos_(i));
      }

      double gap = deltanodalpos.Norm2()-R1_-R2_;
      if(CheckContactStatus(gap) or CheckDampingStatus(gap))
      {
        std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

        std::pair<int,int> integration_ids = std::make_pair(1,1);
        std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

        //Create data container for each end point
        //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
        //has been considered as endpoint, b has the same meaning for element 2
        epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }

    #ifndef ONLYLEFTENDPOINTCONTACT
    if((segid1==0 and boundarynode1_.first) and (segid2==numseg2_-1 and boundarynode2_.second))
    {
      double eta1=-1.0;
      double eta2=1.0;
      LINALG::Matrix<3,1> deltanodalpos(true);
      for(int i=0;i<3;i++)
      {
        deltanodalpos(i)=FADUTILS::CastToDouble(ele2pos_(6+i)-ele1pos_(i));
      }

      double gap = deltanodalpos.Norm2()-R1_-R2_;
      if(CheckContactStatus(gap) or CheckDampingStatus(gap))
      {
        std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

        std::pair<int,int> integration_ids = std::make_pair(1,1);
        std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

        //Create data container for each end point
        //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
        //has been considered as endpoint, b has the same meaning for element 2
        epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }

    if((segid1==numseg1_-1 and boundarynode1_.second) and (segid2==0 and boundarynode2_.first))
    {
      double eta1=1.0;
      double eta2=-1.0;
      LINALG::Matrix<3,1> deltanodalpos(true);
      for(int i=0;i<3;i++)
      {
        deltanodalpos(i)=FADUTILS::CastToDouble(ele2pos_(i)-ele1pos_(6+i));
      }

      double gap = deltanodalpos.Norm2()-R1_-R2_;
      if(CheckContactStatus(gap) or CheckDampingStatus(gap))
      {
        std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

        std::pair<int,int> integration_ids = std::make_pair(1,1);
        std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

        //Create data container for each end point
        //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
        //has been considered as endpoint, b has the same meaning for element 2
        epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }

    if((segid1==numseg1_-1 and boundarynode1_.second) and (segid2==numseg2_-1 and boundarynode2_.second))
    {
      double eta1=1.0;
      double eta2=1.0;
      LINALG::Matrix<3,1> deltanodalpos(true);
      for(int i=0;i<3;i++)
      {
        deltanodalpos(i)=FADUTILS::CastToDouble(ele2pos_(6+i)-ele1pos_(6+i));
      }

      double gap = deltanodalpos.Norm2()-R1_-R2_;
      if(CheckContactStatus(gap) or CheckDampingStatus(gap))
      {
        std::pair<TYPE,TYPE> closestpoint(std::make_pair((TYPE)eta1,(TYPE)eta2));

        std::pair<int,int> integration_ids = std::make_pair(1,1);
        std::pair<int,int> leftpoint_ids = std::make_pair(segid1,segid2);;

        //Create data container for each end point
        //in case of end-point-contact the variable integration_ids contains two bool values (a,b): a \in {0,1} contains the information, if a node of element 1
        //has been considered as endpoint, b has the same meaning for element 2
        epvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,integration_ids,pp,1.0)));
      }
    }
    #endif
  }//for(int i=0;i<closeendpointsegments.size();i++)

}
/*----------------------------------------------------------------------*
 |  end: Get active endpoint pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate active endpoint pairs                           meier 12/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateActiveEndPointPairs(LINALG::SparseMatrix& stiffmatrix,
                                                                                    Epetra_Vector& fint)
{
  for (int numep=0;numep<(int)epvariables_.size();numep++)
  {
    //**********************************************************************
    // (2) Compute some auxiliary quantities
    //**********************************************************************

    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    TYPE eta1=epvariables_[numep]->GetCP().first;
    TYPE eta2=epvariables_[numep]->GetCP().second;

    #ifdef AUTOMATICDIFF
      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
    #endif

    // update shape functions and their derivatives
    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    ComputeNormal(r1, r2, r1_xi, r2_xi, epvariables_[numep],2);

    // call function to compute penalty force
    CalcPenaltyLaw(epvariables_[numep]);

    epvariables_[numep]->SetPPfac(1.0);
    epvariables_[numep]->SetDPPfac(0.0);

    //In case of endpoint-contact, the length specific energy and the 'real' energy are identical
    double lengthspec_energy=FADUTILS::CastToDouble(epvariables_[numep]->GetEnergy());
    epvariables_[numep]->SetIntegratedEnergy(lengthspec_energy);

//    std::cout << "epvariables_[numep]->GetNormal(): " << epvariables_[numep]->GetNormal() << std::endl;
//    std::cout << "numep: " << numep << std::endl;
//    std::cout << "xi: " << epvariables_[numep]->GetCP().first.val() << std::endl;
//    std::cout << "eta: " << epvariables_[numep]->GetCP().second.val() << std::endl;
//    std::cout << "gap: " << epvariables_[numep]->GetGap().val() << std::endl;
//    std::cout << "angle: " << epvariables_[numep]->GetAngle()/M_PI*180.0 << std::endl;
//    std::cout << "r1_xi: " << r1_xi << std::endl;
//    std::cout << "r2_xi: " << r2_xi << std::endl;
//    std::cout << "|r1_xi|: " << r1_xi.Norm2() << std::endl;
//    std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
//    std::cout << "r1_xi*r2_xi: " << FADUTILS::ScalarProduct(r1_xi,r2_xi) << std::endl;
//    std::cout << "epvariables_[numep]->Getfp(): " << epvariables_[numep]->Getfp() << std::endl;

    bool fixedendpointxi = epvariables_[numep]->GetIntIds().first;
    bool fixedendpointeta = epvariables_[numep]->GetIntIds().second;

    // call function to compute contact contribution to residual vector
    EvaluateFcContact(&fint, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, epvariables_[numep],1.0,false,false,fixedendpointxi,fixedendpointeta);

    // call function to compute contact contribution to stiffness matrix
    EvaluateStiffcContact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi,epvariables_[numep],1.0,false,false,fixedendpointxi,fixedendpointeta);
  }
}
/*----------------------------------------------------------------------*
 |  end: Evaluate active endpoint pairs
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate scalar contact force                           meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::CalcPenaltyLaw(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables)
{
    //First parameter for contact force regularization
    double g0 = bcparams_.get<double>("BEAMS_PENREGPARAM_G0",-1.0);
    TYPE fp=0.0;
    TYPE dfp=0.0;
    TYPE e=0.0;
    double pp=variables->GetPP();
    TYPE gap=variables->GetGap();

    if(!CheckContactStatus(FADUTILS::CastToDouble(gap)))
      return;

    switch (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW"))
    {
      case INPAR::BEAMCONTACT::pl_lp:               //linear penalty force law
      {
        fp= - pp*gap;
        dfp=-pp;
        e=-1.0/2.0*pp*gap*gap;

        break;
      }
      case INPAR::BEAMCONTACT::pl_qp:               //quadratic penalty force law
      {
        fp=pp*gap*gap;
        dfp=2*pp*gap;
        e=1.0/3.0*pp*gap*gap*gap;

        break;
      }
      case INPAR::BEAMCONTACT::pl_lnqp:             //quadratic regularization for negative gaps
      {
        if(g0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

        if(gap>-g0)
        {
          //std::cout << "Regularized Penalty!" << std::endl;
          fp=pp/(2.0*g0)*gap*gap;
          dfp=pp/(g0)*gap;
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=-pp*(gap+g0/2.0);
          dfp=-pp;
        }

        break;
      }
      case INPAR::BEAMCONTACT::pl_lpqp:             //quadratic regularization for positiv gaps
      {
        if(g0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

        //Parameter to shift penalty law
        double gbar=bcparams_.get<double>("BEAMS_GAPSHIFTPARAM",0.0);
        gap=gap+gbar;

        double f0=g0*pp/2.0;
        double factor_a=pp/(g0) -f0/(g0*g0); //=pp/(2*g0)
        double factor_b=-pp;
        double factor_c=f0;
        if(gap>0)
        {
          //std::cout << "Regularized Penalty!" << std::endl;
          fp=factor_a*gap*gap+factor_b*gap+factor_c;
          dfp=2*factor_a*gap+factor_b;
          e=-pp*g0*g0/6.0+(pp/(6.0*g0)*gap*gap*gap-pp/2.0*gap*gap+pp*g0/2.0*gap);
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=f0 - pp*gap;
          dfp=-pp;
          e=-(pp*g0*g0/6.0 + pp/2.0*gap*gap - pp*g0/2.0*gap);
        }

        break;
      }
      case INPAR::BEAMCONTACT::pl_lpcp:             //cubic regularization for positive gaps
      {
        if(g0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

        //Third parameter for contact force regularization
        double c0 = bcparams_.get<double>("BEAMS_PENREGPARAM_C0",-1.0);
        if(c0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_C0!");

        //k \in ~[1;3] delivers sensible results representing a parable without turning point
        //k \in ~[3;6] delivers a parable with turning point and consequentely also small negative contact forces ~0.1*f0
        //k=2.0 is  identical to the quadratic regularization for positive gaps!
        double k=c0;
        double f0=pp*g0/k;
        double factor_a=-pp/(g0*g0) +2*f0/(g0*g0*g0);
        double factor_b=2*pp/(g0) -3*f0/(g0*g0);
        double factor_c=-pp;
        double factor_d=f0;
        if(gap>0.0)
        {
          //std::cout << "Regularized Penalty!" << std::endl;
          fp=factor_a*gap*gap*gap+factor_b*gap*gap+factor_c*gap+factor_d;
          dfp=3*factor_a*gap*gap+2*factor_b*gap+factor_c;
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=f0 - pp*gap;
          dfp=-pp;
        }

        break;
      }
      case INPAR::BEAMCONTACT::pl_lpdqp:            //double quadratic regularization for positiv gaps
      {
        if(g0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

        //Third parameter for contact force regularization
        double c0 = bcparams_.get<double>("BEAMS_PENREGPARAM_C0",-1.0);
        if(c0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_C0!");

        //Second parameter for contact force regularization
        double f0 = bcparams_.get<double>("BEAMS_PENREGPARAM_F0",-1.0);
        if(f0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_F0!");

        double k=c0;        //transition between first and second quadratic regularization part: k \in [0;2.0]
        double g1=k*f0/pp;
        double c_tilde=f0;
        double b_tilde=-pp;
        double a_bar=(2*f0-pp*g1)/(2*g0*(g0-g1));
        double b_bar=-2*g0*a_bar;
        double c_bar=-g0*g0*a_bar -g0*b_bar;
        double a_tilde=(2*g1*a_bar+b_bar-b_tilde)/(2*g1);

        if(gap>g1)
        {
          //std::cout << "Regularized Penalty: g1 < gap < g0!" << std::endl;
          fp=a_bar*gap*gap+b_bar*gap+c_bar;
          dfp=2*a_bar*gap+b_bar;
        }
        else if(gap>0)
        {
          //std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
          fp=a_tilde*gap*gap+b_tilde*gap+c_tilde;
          dfp=2*a_tilde*gap+b_tilde;
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=f0 - pp*gap;
          dfp=-pp;
        }

        break;
      }
      case INPAR::BEAMCONTACT::pl_lpep:             //exponential regularization for positiv gaps. Here g0 represents the cut off radius!
      {
        if(g0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

        //Second parameter for contact force regularization
        double f0 = bcparams_.get<double>("BEAMS_PENREGPARAM_F0",-1.0);
        if(f0==-1.0)
          dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_F0!");

        if(gap>0)
        {
          //std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
          fp=f0*exp(-pp*gap/f0);
          dfp=-pp*exp(-pp*gap/f0);
          if(f0*exp(-pp*g0/f0)>0.01*f0)
          {
            std::cout << "      Warning - g0: " << g0 << " f0*exp(-pp*g0/f0): " << f0*exp(-pp*g0/f0) << "-> Choose higher cut-off radius g0!"<< std::endl;
          }
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=f0 - pp*gap;
          dfp=-pp;
        }

        break;
      }
    }

    variables->Setfp(fp);
    variables->Setdfp(dfp);
    variables->SetEnergy(e);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate scalar contact force
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate angle-dependent perp-penalty scale factor       meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::CalcPerpPenaltyScaleFac(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables,
                                                                              LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                              LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                              const double shiftangle1,
                                                                              const double shiftangle2)
{
  //Penalty scale factor that reduces the penalty parameter for small angles
  TYPE ppfac = 1.0;
  TYPE dppfac = 0.0;

  if(shiftangle1>M_PI/2.0 and shiftangle2>M_PI/2.0)
  {
    ppfac = 0.0;
    dppfac = 0.0;
  }
  else
  {
    TYPE s = fabs(FADUTILS::ScalarProduct(r1_xi,r2_xi)/(FADUTILS::VectorNorm<3>(r1_xi)*FADUTILS::VectorNorm<3>(r2_xi)));
    double s1 = cos(shiftangle1);
    double s2 = cos(shiftangle2);

    if(shiftangle1<0.0 or shiftangle1>M_PI/2.0 or shiftangle2<0.0 or shiftangle2>M_PI/2.0 or shiftangle1>=shiftangle2)
      dserror("Invalid choice of shift angles!");

    if(FADUTILS::CastToDouble(s)>s1)
      ppfac=0.0;
    else if(FADUTILS::CastToDouble(s)>s2)
    {
      #ifndef CONSISTENTTRANSITION
        ppfac=0.5*(cos(M_PI*(s-s2)/(s1-s2))+1.0);
        dppfac=-0.5*M_PI/(s1-s2)*sin(M_PI*(s-s2)/(s1-s2));
      #else
        if(CONSISTENTTRANSITION==1)
        {
          TYPE simple_fac=0.5*(cos(M_PI*(s-s2)/(s1-s2))+1.0);
          TYPE d_simple_fac=-0.5*M_PI/(s1-s2)*sin(M_PI*(s-s2)/(s1-s2));
          ppfac=simple_fac*simple_fac;
          dppfac=2*simple_fac*d_simple_fac;
        }
        else if(CONSISTENTTRANSITION==2)
        {
          TYPE simple_fac=0.5*(-cos(M_PI*(s-s2)/(s1-s2))+1.0);
          TYPE d_simple_fac=0.5*M_PI/(s1-s2)*sin(M_PI*(s-s2)/(s1-s2));
          ppfac=1-simple_fac*simple_fac;
          dppfac=-2*simple_fac*d_simple_fac;
        }
        else
          dserror("Inadmissible value of CONSISTENTTRANSITION, only the values 1 and 2 are allowed!");
      #endif
    }
  }

  //set class variable
  cpvariables->SetPPfac(ppfac);
  cpvariables->SetDPPfac(dppfac);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate angle-dependent perp-penalty scale factor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate angle-dependent par-penalty scale factor       Meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::CalcParPenaltyScaleFac(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > gpvariables,
                                                                              LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                              LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                              const double shiftangle1,
                                                                              const double shiftangle2)
{
  //Penalty scale factor that reduces the penalty parameter for small angles
  TYPE ppfac = 1.0;
  TYPE dppfac = 0.0;

  if(shiftangle1>M_PI/2.0 and shiftangle2>M_PI/2.0)
  {
    ppfac = 1.0;
    dppfac = 0.0;
  }
  else
  {
    TYPE s = fabs(FADUTILS::ScalarProduct(r1_xi,r2_xi)/(FADUTILS::VectorNorm<3>(r1_xi)*FADUTILS::VectorNorm<3>(r2_xi)));
    double s1 = cos(shiftangle1);
    double s2 = cos(shiftangle2);

    if(shiftangle1<0.0 or shiftangle1>M_PI/2.0 or shiftangle2<0.0 or shiftangle2>M_PI/2.0 or shiftangle1>=shiftangle2)
      dserror("Invalid choice of shift angles!");

    if(FADUTILS::CastToDouble(s)>s1)
      ppfac=1.0;
    else if(FADUTILS::CastToDouble(s)>s2)
    {
      #ifndef CONSISTENTTRANSITION
        ppfac=0.5*(-cos(M_PI*(s-s2)/(s1-s2))+1.0);
        dppfac=0.5*M_PI/(s1-s2)*sin(M_PI*(s-s2)/(s1-s2));
      #else
        TYPE simple_fac=0.5*(-cos(M_PI*(s-s2)/(s1-s2))+1.0);
        TYPE d_simple_fac=0.5*M_PI/(s1-s2)*sin(M_PI*(s-s2)/(s1-s2));
        ppfac=simple_fac*simple_fac;
        dppfac=2*simple_fac*d_simple_fac;
      #endif
    }
    else
      ppfac=0.0;
  }

  //set class variable
  gpvariables->SetPPfac(ppfac);
  gpvariables->SetDPPfac(dppfac);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate angle-dependent par-penalty scale factor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Subdivide elements into segments for CPP                 meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
double CONTACT::Beam3contact<numnodes, numnodalvalues>::CreateSegments(DRT::Element* ele,
                                                                       std::vector<LINALG::TMatrix<double,3,1> >& endpoints_final,
                                                                       int& numsegment,
                                                                       int i)
{
  //endpoints of the segments
  std::vector<LINALG::TMatrix<double,3,1> > endpoints((int)MAXNUMSEG+1,LINALG::TMatrix<double,3,1>(true));
  double segangle = bcparams_.get<double>("BEAMS_SEGANGLE")/180.0*M_PI;

  numsegment=1;
  double deltaxi=2.0;

  if(i==0)
  {
    numsegment=INITSEG1;
    deltaxi=2.0/INITSEG1;
  }

  if(i==1)
  {
    numsegment=INITSEG2;
    deltaxi=2.0/INITSEG2;
  }

  double xi1(0.0);
  double xi2(0.0);
  LINALG::TMatrix<double,3,1> r1(true);
  LINALG::TMatrix<double,3,1> t1(true);
  LINALG::TMatrix<double,3,1> r2(true);
  LINALG::TMatrix<double,3,1> t2(true);
  LINALG::TMatrix<double,3,1> rm(true);
  double l=0.0;
  double segdist=0.0;
  double maxsegdist=0.0;
  bool moresegments = true;

  while (moresegments)
  {
    //We have to zero maxsegdist for each new segment distribution, otherwise we would get a larger value of a former rougher distribution!
    maxsegdist=0.0;
    moresegments = false;
    for(int i=0;i<numsegment;i++)
    {
      if (numsegment>(int)MAXNUMSEG)
              dserror("Not more segments than MAXNUMSEG per element possible! Increase MAXNUMSEG or apply finer discretization!");

      xi1= 0.0;
      xi2= 0.0;
      xi1= -1.0 +i/((double)numsegment)*2.0;
      xi2= -1.0 +(i+1)/((double)numsegment)*2.0;  //The cast to double is necessary here to avoid integer round-off
      LINALG::TMatrix<TYPE,3,1> auxmatrix(true);

      auxmatrix=r(xi1,ele);
      r1=FADUTILS::CastToDouble<TYPE,3,1>(auxmatrix);
      auxmatrix=r(xi2,ele);
      r2=FADUTILS::CastToDouble<TYPE,3,1>(auxmatrix);
      auxmatrix=r_xi(xi1,ele);
      t1=FADUTILS::CastToDouble<TYPE,3,1>(auxmatrix);
      auxmatrix=r_xi(xi2,ele);
      t2=FADUTILS::CastToDouble<TYPE,3,1>(auxmatrix);
      auxmatrix=r((xi1+xi2)/2.0,ele);
      rm=FADUTILS::CastToDouble<TYPE,3,1>(auxmatrix);

      endpoints[i]=r1;
      endpoints[i+1]=r2;
      l=FADUTILS::VectorNorm<3>(FADUTILS::DiffVector(r1,r2));
      //TODO: adapt this tolerance if necessary!!!
      segdist=1.0*l/2.0*tan(segangle);

      if(segdist>maxsegdist)
        maxsegdist=segdist;
      if(!CheckSegment(r1,t1,r2,t2,rm,segdist))
        moresegments = true;
    }

    deltaxi=deltaxi/2;
    numsegment=numsegment*2;
  }
  numsegment=numsegment/2;

  #ifdef NOSEGMENTATION
    if(numsegment>1)
      dserror("Choose higher SEGANGLE since in the case NOSEGMENTATION only one segment per element is allowed!");
  #endif

  //std::cout << "numsegment: " << numsegment << std::endl;
  endpoints_final.resize(numsegment+1);

  for (int i=0;i<numsegment+1;i++)
  {
    endpoints_final[i]=endpoints[i];
  }

  return maxsegdist;
}
/*----------------------------------------------------------------------*
 |  end: Subdivide elements into segments for CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Max. distance at which a contact force becomes active     meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
double CONTACT::Beam3contact<numnodes, numnodalvalues>::GetMaxActiveDist()
{
  double maxactivedist = 0.0;
  int penaltylaw = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW");

  switch (penaltylaw)
  {
    case INPAR::BEAMCONTACT::pl_lp:
    case INPAR::BEAMCONTACT::pl_qp:
    case INPAR::BEAMCONTACT::pl_lnqp:
    {
      maxactivedist = 0.0;
      break;
    }
    case INPAR::BEAMCONTACT::pl_lpqp:
    {
      double g0 = bcparams_.get<double>("BEAMS_PENREGPARAM_G0",-1.0);
      if(g0==-1.0)
        dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

      //Parameter to shift penalty law
      double gbar=bcparams_.get<double>("BEAMS_GAPSHIFTPARAM",0.0);

      maxactivedist = g0 - gbar;

      break;
    }
    case INPAR::BEAMCONTACT::pl_lpcp:
    case INPAR::BEAMCONTACT::pl_lpdqp:
    case INPAR::BEAMCONTACT::pl_lpep:
    {
      maxactivedist = bcparams_.get<double>("BEAMS_PENREGPARAM_G0",-1.0);
      if(maxactivedist==-1.0)
        dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");
      break;
    }
  }
  if (DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_DAMPING")!=INPAR::BEAMCONTACT::bd_no)
  {
    double gd1 = bcparams_.get<double>("BEAMS_DAMPREGPARAM1",-1000.0);
    if(gd1==-1000.0)
      dserror("Damping parameter BEAMS_DAMPINGPARAM, BEAMS_DAMPREGPARAM1 and BEAMS_DAMPREGPARAM2 have to be chosen!");
    if (gd1>maxactivedist)
      maxactivedist=gd1;
  }

  return maxactivedist;
}

/*----------------------------------------------------------------------*
 |  Check, if segments are fine enough                       meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::CheckSegment( LINALG::TMatrix<double,3,1>& r1,
                                                                    LINALG::TMatrix<double,3,1>& t1,
                                                                    LINALG::TMatrix<double,3,1>& r2,
                                                                    LINALG::TMatrix<double,3,1>& t2,
                                                                    LINALG::TMatrix<double,3,1>& rm,
                                                                    double& segdist)
{
  LINALG::TMatrix<double,3,1> t_lin(true);
  LINALG::TMatrix<double,3,1> rm_lin(true);
  double angle1(0.0);
  double angle2(0.0);
  double dist(0.0);
  double segangle = bcparams_.get<double>("BEAMS_SEGANGLE")/180.0*M_PI;

  //Calculate tangent and midpint of linear nodal interpolation
  for(int i=0;i<3;i++)
  {
    t_lin(i)=r2(i)-r1(i);
    rm_lin(i)=(r2(i)+r1(i))/2.0;
  }

  LINALG::TMatrix<double,3,1> diffvec(true);
  diffvec = FADUTILS::DiffVector(rm_lin, rm);
  dist = (double)FADUTILS::VectorNorm<3>(diffvec);
  angle1 = (double)BEAMCONTACT::CalcAngle(t1,t_lin);
  angle2 = (double)BEAMCONTACT::CalcAngle(t2,t_lin);

  if(fabs(angle1)<segangle and fabs(angle2)<segangle) //segment distribution is fine enough
  {
    if(fabs(dist)>segdist)
      dserror("Value of segdist too large, approximation as circle segment not possible!");

    return true;
  }
  else                                                                       //we still need more segments
    return false;
}

/*----------------------------------------------------------------------*
 |  Find segments close to each other                        meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetCloseSegments( const std::vector<LINALG::TMatrix<double,3,1> >& endpoints1,
                                                                        const std::vector<LINALG::TMatrix<double,3,1> >& endpoints2,
                                                                        std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closesmallanglesegments,
                                                                        std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closelargeanglesegments,
                                                                        std::vector<std::pair<int,int> >& closeendpointsegments,
                                                                        double maxactivedist)
{

  LINALG::TMatrix<double,3,1> t1(true);
  LINALG::TMatrix<double,3,1> t2(true);
  LINALG::TMatrix<double,3,1> r1_a(true);
  LINALG::TMatrix<double,3,1> r1_b(true);
  LINALG::TMatrix<double,3,1> r2_a(true);
  LINALG::TMatrix<double,3,1> r2_b(true);
  double angle(0.0);

  bool endpoint_penalty=DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_ENDPOINTPENALTY");

  //Safety factor for determination of close segments
  double safetyfac = 1.1;
  //Distance at which intersection happens
  double distancelimit=safetyfac*(maxsegdist1_+maxsegdist2_+maxactivedist+R1_+R2_);

  int numseg1 = (int)endpoints1.size()-1;
  int numseg2 = (int)endpoints2.size()-1;

  //TODO: This check is implemented in a brute force way. However, this should be efficient enough as long
  //as the number of segments per element remains small!
  for (int i=0;i<numseg1;i++)
  {
    r1_a=endpoints1[i];
    r1_b=endpoints1[i+1];
    t1=FADUTILS::DiffVector(r1_b,r1_a);
    for (int j=0;j<numseg2;j++)
    {
      r2_a=endpoints2[j];
      r2_b=endpoints2[j+1];
      t2=FADUTILS::DiffVector(r2_b,r2_a);

      angle=BEAMCONTACT::CalcAngle(t1,t2);

      //*******1) intersection between two parallel cylinders*********************************************************
      if(fabs(angle)<ANGLETOL)
      {
        if(BEAMCONTACT::IntersectParallelCylinders(r1_a,r1_b,r2_a,r2_b,distancelimit))
        {
          LINALG::TMatrix<double,3,1> segmentdata(true);
          segmentdata(0)=angle;  //segment angle
          segmentdata(1)=1000.0;    //eta1_seg
          segmentdata(2)=1000.0;    //eta2_seg

          //Add new small angle pair
          if(fabs(angle)<=deltasmallangle_)
            closesmallanglesegments[std::make_pair(i,j)]=segmentdata;
          if(fabs(angle)>=deltalargeangle_)
            closelargeanglesegments[std::make_pair(i,j)]=segmentdata;

          //If the element lies on the boundary of a physical beam, we sort out the corresponding boundary segments
          if(endpoint_penalty)
          {
            if((i==0 and boundarynode1_.first) or (i==numseg1-1 and boundarynode1_.second) or (j==0 and boundarynode2_.first) or (j==numseg2-1 and boundarynode2_.second))
              closeendpointsegments.push_back(std::make_pair(i,j));
          }
        }

      }
      //*******2) intersection between two arbitrary oriented cylinders************************************************
      else
      {
        std::pair<double,double> closestpoints(std::make_pair(0.0,0.0));
        bool etaset=false;
        if(BEAMCONTACT::IntersectArbitraryCylinders(r1_a,r1_b,r2_a,r2_b,distancelimit,closestpoints,etaset))
        {
          LINALG::TMatrix<double,3,1> segmentdata(true);
          segmentdata(0)=angle;                  //segment angle

          if(etaset)
          {
            segmentdata(1)=closestpoints.first;    //eta1_seg
            segmentdata(2)=closestpoints.second;   //eta2_seg
          }
          else
          {
            segmentdata(1)=1000.0;   //eta1_seg
            segmentdata(2)=1000.0;   //eta2_seg
          }
          //Add new small angle pair
          if(fabs(angle)<=deltasmallangle_)
            closesmallanglesegments[std::make_pair(i,j)]=segmentdata;
          if(fabs(angle)>=deltalargeangle_)
            closelargeanglesegments[std::make_pair(i,j)]=segmentdata;

          //If the element lies on the boundary of a physical beam, we sort out the corresponding boundary segments
          if(endpoint_penalty)
          {
            if((i==0 and boundarynode1_.first) or (i==numseg1-1 and boundarynode1_.second) or (j==0 and boundarynode2_.first) or (j==numseg2-1 and boundarynode2_.second))
              closeendpointsegments.push_back(std::make_pair(i,j));
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Closest point projection                                  meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::ClosestPointProjection( double& eta_left1,
                                                                              double& eta_left2,
                                                                              double& l1,
                                                                              double& l2,
                                                                              LINALG::TMatrix<double,3,1>& segmentdata,
                                                                              std::pair<TYPE,TYPE>& solutionpoints,
                                                                              int segid1, int segid2)
{

  std::vector<std::pair<double,double> > startingpoints(0);
  bool validpairfound=false;
  double gap=0.0;
  double eta_right1=eta_left1+l1;
  double eta_right2=eta_left2+l2;

  double etalocal1 = segmentdata(1);
  double etalocal2 = segmentdata(2);

  if(fabs(etalocal1)<=1.0 and fabs(etalocal2)<=1.0)
  {
    startingpoints.push_back(std::make_pair(eta_left1+0.5*l1*(1+etalocal1),eta_left2+0.5*l2*(1+etalocal2))); //cp of linear segment approximation as starting point
  }

  startingpoints.push_back(std::make_pair(eta_left1+0.5*l1,eta_left2+0.5*l2));//segment midpoint as starting point

  //Other combinations of (etalocal1, etalocal2 \in {-1;0;1}) for each segment -> 8 additional combinations besides (0,0)
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      if(!(i==0 and j==0)) //we already have the segment midpoint combination (0,0)
        startingpoints.push_back(std::make_pair(eta_left1+i*0.5*l1,eta_left2+j*0.5*l2));
    }
  }

  for (int numstartpoint=0;numstartpoint<(int)startingpoints.size();numstartpoint++)
  {
    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    //Tangent and derivatives for tangent field smoothing (only for Reissner beams)
    LINALG::TMatrix<TYPE,3,1> t1(true);
    LINALG::TMatrix<TYPE,3,1> t1_xi(true);
    LINALG::TMatrix<TYPE,3,1> t2(true);
    LINALG::TMatrix<TYPE,3,1> t2_xi(true);

    // initialize function f and Jacobian df for Newton iteration
    LINALG::TMatrix<TYPE,2,1> f(true);
    LINALG::TMatrix<TYPE,2,2> df(true);
    LINALG::TMatrix<TYPE,2,2> dfinv(true);

    // initial scalar residual (L2-norm of f)
    double residual = 0.0;
    double lastresidual = 0.0;
    double residual0 = 0.0;
    int iter=0;

    TYPE eta1=startingpoints[numstartpoint].first;
    TYPE eta2=startingpoints[numstartpoint].second;
    double eta1_old=FADUTILS::CastToDouble(eta1);
    double eta2_old=FADUTILS::CastToDouble(eta2);
    bool converged = false;
    bool elementscolinear = false;

    #ifdef FADCHECKS
      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
    #endif

    //std::cout << "numstartpoint: " << numstartpoint << std::endl;

    //**********************************************************************
    // local Newton iteration
    //**********************************************************************
    for (int i=0;i<BEAMCONTACTMAXITER;++i)
    {
      //store residual of last iteration
      lastresidual=residual;
      iter++;

      // reset shape function variables to zero
      N1.Clear();
      N2.Clear();
      N1_xi.Clear();
      N2_xi.Clear();
      N1_xixi.Clear();
      N2_xixi.Clear();

      // update shape functions and their derivatives
      GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
      // update coordinates and derivatives of contact points
      ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

      // use delta_r = r1-r2 as auxiliary quantity
      delta_r=FADUTILS::DiffVector(r1,r2);

      // compute norm of difference vector to scale the equations
      // (this yields better conditioning)
      // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
      // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
      double norm_delta_r = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(delta_r));
      gap=norm_delta_r-R1_-R2_;

      // The closer the beams get, the smaller is norm_delta_r, but
      // norm_delta_r is not allowed to be too small, else numerical problems occur.
      // It can happen quite often that the centerlines of two beam elements of the same physical beam
      // cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they will be sorted out later anyways.
      //std::cout << "norm_delta_r: " << norm_delta_r << std::endl;
      if (norm_delta_r < NORMTOL)
      {
        // this exludes pairs with IDs i and i+2, i.e. contact with the next but one element
        if (FADUTILS::CastToDouble(FADUTILS::Norm(eta1)) <= 1.0 and FADUTILS::CastToDouble(FADUTILS::Norm(eta2)) <= 1.0)
        {
          dserror("Beam axis identical, choose smaller time step!");
        }
        else
        {
          break;
        }
      }

      int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(bcparams_,"BEAMS_SMOOTHING");
      if (smoothing != INPAR::BEAMCONTACT::bsm_none) //smoothed case
      {
        // Evaluate nodal tangents in each case. However, they are used only if smoothing=INPAR::BEAMCONTACT::bsm_cpp
        CONTACT::B3TANGENTSMOOTHING::ComputeTangentsAndDerivs<numnodes,numnodalvalues>(t1, t1_xi, nodaltangentssmooth1_, N1, N1_xi);
        CONTACT::B3TANGENTSMOOTHING::ComputeTangentsAndDerivs<numnodes,numnodalvalues>(t2, t2_xi, nodaltangentssmooth2_, N2, N2_xi);
      }

      // evaluate f at current eta1, eta2
      EvaluateOrthogonalityCondition(f, delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);

      TYPE jacobi1 = 1.0;
      TYPE jacobi2 = 1.0;
      jacobi1 =  GetJacobi(element1_);
      jacobi2 =  GetJacobi(element2_);

      // compute the scalar residuum
      // The residual is scaled with 1/element_length since an absolute
      // residual norm is used as local CPP convergence criteria and r_xi scales with the element_length
       residual = sqrt((double)FADUTILS::CastToDouble((TYPE)(f(0)*f(0)/(jacobi1*jacobi1) + f(1)*f(1)/(jacobi2*jacobi2))));

//      std::cout << "iter: " << iter << std::endl;
//      std::cout << "residual: " << residual << std::endl;
//      std::cout << "eta1: " << eta1 << std::endl;
//      std::cout << "eta2: " << eta2 << std::endl;
//      std::cout << "angle: " << BEAMCONTACT::CalcAngle(r1_xi,r2_xi) << std::endl;
//      std::cout << "iter: " << iter << std::endl;
//      std::cout << "residual: " << residual << std::endl;

      if(iter==1)
        residual0=residual;

      // check if Newton iteration has converged
      #ifndef RELBEAMCONTACTTOL
        if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL and fabs(eta1_old-FADUTILS::CastToDouble(eta1))<XIETAITERATIVEDISPTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
        {
          converged=true;
          break;
        }
      #else
        if(residual0>1.0e-6)
        {
          if (FADUTILS::CastToDouble(residual/residual0) < RELBEAMCONTACTTOL and fabs(eta1_old-FADUTILS::CastToDouble(eta1))<XIETAITERATIVEDISPTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
          {
            converged=true;
            break;
          }
        }
        else
        {
          if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL and fabs(eta1_old-FADUTILS::CastToDouble(eta1))<XIETAITERATIVEDISPTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
          {
            converged=true;
            break;
          }
        }
      #endif

      // evaluate Jacobian of f at current eta1, eta2
      // Note: Parallel elements can not be handled with this beam contact formulation;
      EvaluateLinOrthogonalityCondition(df, dfinv, delta_r, norm_delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, t1, t2, t1_xi, t2_xi, elementscolinear);

      #ifdef FADCHECKS
        std::cout << "f: " << f << std::endl;
        std::cout << "df: " << df << std::endl;
        BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
        FADCheckLinOrthogonalityCondition(delta_r, norm_delta_r, r1_xi,r2_xi, t1, t2);
      #endif

      if (elementscolinear)
      {
        break;
      }

      eta1_old=FADUTILS::CastToDouble(eta1);
      eta2_old=FADUTILS::CastToDouble(eta2);

      //update element coordinates of contact point
      eta1 += -dfinv(0,0)*f(0) - dfinv(0,1)*f(1);
      eta2 += -dfinv(1,0)*f(0) - dfinv(1,1)*f(1);

      #ifdef FADCHECKS
        BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
      #endif

      //std::cout << "eta1: " << eta1 << std::endl;
      //std::cout << "eta2: " << eta2 << std::endl;
    }//for (int i=0;i<BEAMCONTACTMAXITER;++i)
      //**********************************************************************

    // Newton iteration unconverged after BEAMCONTACTMAXITER
    if (!converged)
    {
      //Initialize g_min with a very large value, at which no active conact should occur!
      double g_min=1000*R2_;
      if(CheckContactStatus(g_min) or CheckDampingStatus(g_min))
        dserror("Are sure that contact should be active at such large gaps?");

      double alpha_g_min=0.0;
      //In case no valid point-to-line solution is found in () (pointtolinesolfound=false)
      //it is assumed, that the distance between the segments is large enough such that no contact can occur.
      bool pointtolinesolfound=false;
      double eta1_min=0.0;
      double eta2_min=0.0;

      CheckUnconvergedSegmentPair(eta_left1,eta_left2,l1,l2,eta1_min,eta2_min,g_min,alpha_g_min,pointtolinesolfound);

      //Check, if we have found a valid point-to-line projection and if the solution is not a boundary minimum
      if(pointtolinesolfound and eta1_min>(eta_left1-1.0e-10) and eta1_min<(eta_left1+l1+1.0e-10) and eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {

        double perpshiftangle1 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE1")/180.0*M_PI;
        //Here, we apply the conservative estimate that the closest-point gap is by 0.1*R2_ smaller than g_min
        double g_min_estimate=g_min-0.1*R2_;

        //TODO
        if((CheckContactStatus(g_min_estimate) or CheckDampingStatus(g_min_estimate)) and fabs(alpha_g_min)>=perpshiftangle1)
        {
          std::cout << std::endl << "Serious Warning!!!!! Local CPP not converged: CP-Approximation applied!" << std::endl;
          std::cout << "element1_->Id(): " << element1_->Id() << std::endl;
          std::cout << "element2_->Id(): " << element2_->Id() << std::endl;
          std::cout << "R2_: " << R2_ << std::endl;
          std::cout << "g_min: " << g_min << std::endl;
          std::cout << "alpha_g_min: " << alpha_g_min/M_PI*180 << "°"<< std::endl;
          std::cout << "numstartpoint: " << numstartpoint << std::endl;
          std::cout << "iter: " << iter << std::endl;
          std::cout << "residual0: " << residual0 << std::endl;
          std::cout << "lastresidual: " << lastresidual << std::endl;
          std::cout << "residual: " << residual << std::endl;
          std::cout << "eta1_min: " << FADUTILS::CastToDouble(eta1_min) << std::endl;
          std::cout << "eta1: " << FADUTILS::CastToDouble(eta1) << std::endl;
          std::cout << "eta1_old: " << eta1_old << std::endl;
          std::cout << "eta2_min: " << FADUTILS::CastToDouble(eta2_min) << std::endl;
          std::cout << "eta2: " << FADUTILS::CastToDouble(eta2) << std::endl;
          std::cout << "eta2_old: " << eta2_old << std::endl;

          //We need here the original elements of the problem discretization in order to read out time-dependent
          //element quantities (such as the current curvature) since element1_ and element2_ are pure copies of
          //these elements generated once in the beginning of the simulation and which do therefore not contain
          //the current values of such element quantities
          DRT::Element* element1 = pdiscret_.lColElement(pdiscret_.ElementColMap()->LID(element1_->Id()));
          DRT::Element* element2 = pdiscret_.lColElement(pdiscret_.ElementColMap()->LID(element2_->Id()));

          const DRT::ElementType & eot = element1->ElementType();
          if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
          {
            const DRT::ELEMENTS::Beam3eb* beam3ebelement1 = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element1);
            double kappamax1 = beam3ebelement1->GetKappaMax();
            const DRT::ELEMENTS::Beam3eb* beam3ebelement2 = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element2);
            double kappamax2 = beam3ebelement2->GetKappaMax();

            std::cout << "kappamax1: " << kappamax1 << std::endl;
            std::cout << "kappamax2: " << kappamax2 << std::endl << std::endl;

          }

          //Apply Point-To-Line solution as approximation for CPP or...
          #ifdef CPP_APPROX
            eta1=eta1_min;
            eta2=eta2_min;
            gap=g_min;
            double angle=alpha_g_min;

            if(fabs(eta1-1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta1+1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta2-1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta2+1.0)<1.1*XIETAITERATIVEDISPTOL)
                dserror("|eta1|=1 or |eta2|=1, danger of multiple gauss point evaluation!");

            double perpshiftangle1 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE1")/180.0*M_PI;

            if((CheckContactStatus(gap) or CheckDampingStatus(gap)) and angle >=perpshiftangle1)
              validpairfound=true;

            solutionpoints.first=FADUTILS::CastToDouble(eta1);
            solutionpoints.second=FADUTILS::CastToDouble(eta2);

            break;

          //... or abort simulation
          #else
            dserror("CPP is not converged, eventhough the corresponding closest point is active! "
                    "Decrease the value of SEGANGLE or increase your shifting angles!");
          #endif
        }
      }

      eta1 = 1e+12;
      eta2 = 1e+12;
    }
    else
    {

      //if we have already found a converged solution with valid closest points eta1 \in [eta_left1;eta_right1] and
      //eta2 \in [eta_left2;eta_right2], we can finish here and don't have to apply more starting points
      if( eta_left1 - XIETAITERATIVEDISPTOL <= eta1 and eta1 <= eta_right1 + XIETAITERATIVEDISPTOL and eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {
        if(fabs(eta1-1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta1+1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta2-1.0)<1.1*XIETAITERATIVEDISPTOL or fabs(eta2+1.0)<1.1*XIETAITERATIVEDISPTOL)
            dserror("|eta1|=1 or |eta2|=1, danger of multiple gauss point evaluation!");

        if(FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(r1_xi))<1.0e-8 or FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(r2_xi))<1.0e-8)
          dserror("Tangent vector of zero length, choose smaller time step!");

        double angle=fabs(BEAMCONTACT::CalcAngle(FADUTILS::CastToDouble<TYPE,3,1>(r1_xi),FADUTILS::CastToDouble<TYPE,3,1>(r2_xi)));

        double perpshiftangle1 = bcparams_.get<double>("BEAMS_PERPSHIFTANGLE1")/180.0*M_PI;

        if((CheckContactStatus(gap) or CheckDampingStatus(gap)) and angle >=perpshiftangle1)
          validpairfound=true;

        //Here, we perform an additional security check: If a unique CCP solution exists, the Newton scheme should find it
        //with the first starting point. Otherwise, the problem may be ill-conditioned!
        if(validpairfound and numstartpoint!=0)
          dserror("Valid CCP solution has not been found with the first starting point. Choose smaller value of SEGANGLE!");

        solutionpoints.first=FADUTILS::CastToDouble(eta1);
        solutionpoints.second=FADUTILS::CastToDouble(eta2);

        break;
      }
    }
  }//for (int numstartpoint=0;numstartpoint<startingpoints.size();numstartpoint++)

  return validpairfound;

}
/*----------------------------------------------------------------------*
|  end: Closest point projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Closest Point-To-Line Projection                         meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::PointToLineProjection(double& eta1_slave,
                                                                            double& eta_left2,
                                                                            double& l2,
                                                                            double& eta2_master,
                                                                            double& gap,
                                                                            double& alpha,
                                                                            bool& pairactive,
                                                                            bool smallanglepair,
                                                                            bool invertpairs,
                                                                            bool orthogonalprojection)
{

  /* Attention: With the paramters invertpairs and orthogonalprojection, 4 different types of projections can be realized.
   * The parameter invertpairs=true simply changes the meaning of eta1 and eta2 within this method, i.e. in this case the given
   * parameter coordinate is eta2 of beam2 and the searched parameter is eta1 on beam1 while in the standard case eta1 is given
   * and eta2 is searched. The parameter orthogonalprojection changes the projection method: In the standard case (orthogonalprojection=false),
   * the orthogonality condition is satisfied at the projection side of the searched parameter (i.e. on beam2 if invertpairs=false or beam1 if inverpairs=true)
   * while for orthogonalprojection=true the orthogonality condition is satisfied at the projection side of the given parameter
   * (i.e. on beam1 if invertpairs=false or beam2 if inverpairs=true). This leads to the following four possible projections:
   * |invertpairs|orthogonalprojection|given parameter|searched parameter|orthogonality condition satisfied on|
   * |   false   |       false        |      eta1     |      eta2        |              beam2                 | -> default case
   * |   true    |       false        |      eta2     |      eta1        |              beam1                 |
   * |   false   |       true         |      eta1     |      eta2        |              beam1                 |
   * |   true    |       true         |      eta2     |      eta1        |              beam2                 |
   */
  std::vector<double> startingpoints(3,0.0);
  double gap_test=0.0;
  double eta_right2=eta_left2+l2;

  //We start with the segment midpoint, since there it is most likely to find a closest point solution within the segment
  startingpoints[0]=eta_left2+0.5*l2; //segment midpoint as starting point
  startingpoints[1]=eta_left2; //left segment point as starting point
  startingpoints[2]=eta_left2+l2; //right segment point as starting point

  for (int numstartpoint=0;numstartpoint<(int)startingpoints.size();numstartpoint++)
  {

    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    // initialize function f and Jacobian df for Newton iteration
    TYPE f=0.0;
    TYPE df=0.0;
    LINALG::TMatrix<TYPE,2,2> dfinv(true);

    // initial scalar residual (L2-norm of f)
    double residual = 0.0;
    double lastresidual = 0.0;
    double residual0 = 0.0;
    int iter=0;

    TYPE eta1=eta1_slave;
    TYPE eta2=startingpoints[numstartpoint];
    double eta2_old=FADUTILS::CastToDouble(eta2);

    bool converged = false;

    #ifdef FADCHECKS
      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
    #endif

    //**********************************************************************
    // local Newton iteration
    //**********************************************************************
    for (int i=0;i<BEAMCONTACTMAXITER;++i)
    {
      //store residual of last iteration
      lastresidual=residual;
      iter++;

      // reset shape function variables to zero
      N1.Clear();
      N2.Clear();
      N1_xi.Clear();
      N2_xi.Clear();
      N1_xixi.Clear();
      N2_xixi.Clear();

      bool inversion_possible=false;
      bool endpointpenalty = DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_ENDPOINTPENALTY");
      if(endpointpenalty) inversion_possible=true;

      #ifdef ENDPOINTSEGMENTATION
      inversion_possible=true;
      #endif

      if(inversion_possible)
      {
        //In the case of ENDPOINTSEGMENTATION or ENDPOINTPENALTY it can be necessary to make an invere projection (from the master beam onto the slave beam). In this
        //case, the local variables (e.g. r1, r1_xi...) inside PointToLineProjection() with index 1 represent the master beam which has the global index 2. In order
        //to get the right nodal positions ele2pos_ for the local variables r1, r1_xi, r1_xixi, we have to invert the arguments of the function call ComputeCoordsAndDerivs()!
        if(!invertpairs)
        {
          //update shape functions and their derivatives
          GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
          //update coordinates and derivatives of contact points
          ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
        }
        else
        {
          //update shape functions and their derivatives
          GetShapeFunctions(N2, N1, N2_xi, N1_xi, N2_xixi, N1_xixi, eta2, eta1);
          //update coordinates and derivatives of contact points
          ComputeCoordsAndDerivs(r2, r1, r2_xi, r1_xi, r2_xixi, r1_xixi, N2, N1, N2_xi, N1_xi, N2_xixi, N1_xixi);
        }
      }
      else
      {
        // update shape functions and their derivatives
        GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
        // update coordinates and derivatives of contact points
        ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
      }

      // use delta_r = r1-r2 as auxiliary quantity
      delta_r=FADUTILS::DiffVector(r1,r2);

      // compute norm of difference vector to scale the equations
      // (this yields better conditioning)
      // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
      // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
      double norm_delta_r = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(delta_r));
      gap_test=norm_delta_r-R1_-R2_;

      // The closer the beams get, the smaller is norm_delta_r, but
      // norm_delta_r is not allowed to be too small, else numerical problems occur.
      // It can happen quite often that the centerlines of two beam elements of the same physical beam
      // cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they will be sorted out later anyways.
      //std::cout << "norm_delta_r: " << norm_delta_r << std::endl;
      if (norm_delta_r < NORMTOL)
      {
        // this exludes pairs with IDs i and i+2, i.e. contact with the next but one element
        if (FADUTILS::CastToDouble(FADUTILS::Norm(eta1)) + FADUTILS::CastToDouble(FADUTILS::Norm(eta2)) < NEIGHBORTOL)
        {
          dserror("Beam axis identical, choose smaller time step!");
        }
      }

      // evaluate f at current eta1, eta2
      EvaluatePTLOrthogonalityCondition(f, delta_r, norm_delta_r, r1_xi, r2_xi, orthogonalprojection);

      //The residual will be scaled with the length of the element whose tangent appears in the orthogonality conditions.
      //Which of the two elements is relevant in this context depends of the parameters orthogonalprojection and invertpairs:
      TYPE jacobi = 1.0;

      if(invertpairs==false and orthogonalprojection==false)//default case
      {
        jacobi =  GetJacobi(element2_);
      }
      else if(invertpairs==false and orthogonalprojection==true)
      {
        jacobi =  GetJacobi(element1_);
      }
      else if(invertpairs==true and orthogonalprojection==false)
      {
        jacobi =  GetJacobi(element1_);
      }
      else if(invertpairs==true and orthogonalprojection==true)
      {
        jacobi =  GetJacobi(element2_);
      }
      jacobi =  GetJacobi(element2_);

      // compute the scalar residuum
      // The residual is scaled with 1/element_length since r_xi scales with the element_length
      residual = fabs((double)FADUTILS::CastToDouble((TYPE)(f/jacobi)));

//      std::cout << "iter: " << iter << std::endl;
//      std::cout << "residual: " << residual << std::endl;
//      std::cout << "eta1: " << eta1.val() << std::endl;
//      std::cout << "eta2: " << eta2.val() << std::endl;
//      std::cout << "r1: " << r1 << std::endl;
//      std::cout << "r2: " << r2 << std::endl;
//      std::cout << "r1_xi: " << r1_xi << std::endl;
//      std::cout << "r2_xi: " << r2_xi << std::endl;
//      std::cout << "r1_xixi: " << r1_xixi << std::endl;
//      std::cout << "r2_xixi: " << r2_xixi << std::endl;
//      std::cout << "ele1pos_: " << ele1pos_ << std::endl;
//      std::cout << "ele2pos_: " << ele2pos_ << std::endl;
//      std::cout << "angle: " << BEAMCONTACT::CalcAngle(FADUTILS::CastToDouble<TYPE,3,1>(r1_xi),FADUTILS::CastToDouble<TYPE,3,1>(r2_xi))/M_PI*180.0 << std::endl;

      if(iter==1)
        residual0=residual;

      // check if Newton iteration has converged
      #ifndef RELBEAMCONTACTTOL
        if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
        { converged=true;
          break;
        }
      #else
        if(residual0>1.0e-6)
        {
          if (FADUTILS::CastToDouble(residual/residual0) < RELBEAMCONTACTTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
          { converged=true;
            break;
          }
        }
        else
        {
          if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL and fabs(eta2_old-FADUTILS::CastToDouble(eta2))<XIETAITERATIVEDISPTOL)
          { converged=true;
            break;
          }
        }
      #endif

      // evaluate Jacobian of f at current eta1, eta2
      // Note: It has to be checked, if the linearization is equal to zero;
      bool validlinearization = EvaluateLinPTLOrthogonalityCondition(df, delta_r, norm_delta_r, r1_xi, r2_xi, r2_xixi, orthogonalprojection);

      if(!validlinearization)
        dserror("Linearization of point to line projection is zero, choose tighter search boxes!");

      #ifdef FADCHECKS
        std::cout << "f: " << f << std::endl;
        std::cout << "df: " << df << std::endl;
        BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
        FADCheckLinOrthogonalityCondition(delta_r, norm_delta_r, r1_xi,r2_xi, t1, t2);
      #endif

      eta2_old=FADUTILS::CastToDouble(eta2);

      //update master element coordinate of contact point
      eta2 += -f/df;

      #ifdef FADCHECKS
        BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
      #endif

    }//for (int i=0;i<BEAMCONTACTMAXITER;++i)
      //**********************************************************************

    // Newton iteration unconverged after BEAMCONTACTMAXITER
    if (!converged)
    {
      std::cout << "iter: " << iter << std::endl;
      std::cout << "residual0: " << residual0 << std::endl;
      std::cout << "lastresidual: " << lastresidual << std::endl;
      std::cout << "residual: " << residual << std::endl;
      std::cout << "eta2: " << FADUTILS::CastToDouble(eta2) << std::endl;
      std::cout << "eta2_old: " << eta2_old << std::endl;


      //TODO:
      dserror("Local Newton loop unconverged. Adapt segangle or the shift angles for small-anlge contact!");

      eta1 = 1e+12;
      eta2 = 1e+12;
    }
    else
    {

      //if we have already found a converged solution with valid closest point eta2 \in [eta_left2;eta_right2], we can finish here and don't have to apply more starting points
      if(eta_left2 - XIETAITERATIVEDISPTOL <= eta2 and eta2 <= eta_right2 + XIETAITERATIVEDISPTOL)
      {

        if(fabs(eta2-1.0)<1.1*XIETAITERATIVEDISPTOL)
        {
          bool throw_error=true;

          //There is no danger of multiple evaluation, if the considered node is a boundary node of the beam:
          if((invertpairs==false and boundarynode2_.second) or (invertpairs==true and boundarynode1_.second))
            throw_error=false;

          if(throw_error)
          {
            std::cout << "ID1: " << element1_->Id() << std::endl;
            std::cout << "ID2: " << element2_->Id() << std::endl;
            std::cout << "eta1: " << eta1 << std::endl;
            std::cout << "eta2: " << eta2 << std::endl;
            //TODO: In some cases a warning is sufficient, but in general we need the dserror("");
            std::cout << "Serious Warning!!!!! eta2=1, danger of multiple gauss point evaluation! " << std::endl;
            //dserror("eta2=1, danger of multiple gauss point evaluation!");
          }
        }

        if(fabs(eta2+1.0)<1.1*XIETAITERATIVEDISPTOL)
        {
          bool throw_error=true;

          //There is no danger of multiple evaluation, if the considered node is a boundary node of the beam:
          if((invertpairs==false and boundarynode2_.first) or (invertpairs==true and boundarynode1_.first))
            throw_error=false;

          if(throw_error)
            dserror("eta2=-1, danger of multiple gauss point evaluation!");
        }

        if(FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(r1_xi))<1.0e-8 or FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(r2_xi))<1.0e-8)
          dserror("Tangent vector of zero length, choose smaller time step!");

        bool relevant_angle=true;
        double angle=fabs(BEAMCONTACT::CalcAngle(FADUTILS::CastToDouble<TYPE,3,1>(r1_xi),FADUTILS::CastToDouble<TYPE,3,1>(r2_xi)));
        if(smallanglepair)
        {
          double parshiftangle2 = bcparams_.get<double>("BEAMS_PARSHIFTANGLE2")/180.0*M_PI;

          if(angle > parshiftangle2)
            relevant_angle=false;
        }
        if((CheckContactStatus(gap_test) or CheckDampingStatus(gap_test)) and relevant_angle)
          pairactive=true;

        eta2_master=FADUTILS::CastToDouble(eta2);

        //Here, we perform an additional security check: If a unique CCP solution exists, the Newton scheme should find it
        //with the first starting point. Otherwise, the problem may be ill-conditioned!
        if(pairactive and numstartpoint!=0)
          dserror("Valid Point-To-Line solution has not been found with the first starting point. Choose smaller value of SEGANGLE!");

        gap=gap_test;
        alpha=angle;

        //eta2 \in [eta_left2;eta_right2] --> true
        return true;
      }
    }
  }//for (int numstartpoint=0;numstartpoint<startingpoints.size();numstartpoint++)

  //no eta2 \in [eta_left2;eta_right2] --> false
  return false;

}
/*----------------------------------------------------------------------*
|  end: Point-To-Line Projection
 *----------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------*
|  Determine minimal distance and contact angle for unconverged segment pair     meier 05/15|
 *------------------------------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::CheckUnconvergedSegmentPair( double& eta_left1,
                                                                                   double& eta_left2,
                                                                                   double& l1,
                                                                                   double& l2,
                                                                                   double& eta1_min,
                                                                                   double& eta2_min,
                                                                                   double& g_min,
                                                                                   double& alpha_g_min,
                                                                                   bool& pointtolinesolfound)
{

  //Calculate initial length of slave element
  LINALG::TMatrix<double,3,1> lengthvec1(true);
  for(int i=0;i<3;i++)
  {
    lengthvec1(i)=(element1_->Nodes())[0]->X()[i]-(element1_->Nodes())[1]->X()[i];
  }
  //length1 = physical length; l1=length in parameter space
  double length1=lengthvec1.Norm2();

  int n=1;
  //subdivide the slave segment by n+1 test points until the distance between the
  //test points is smaller than half of the cross-section radius
  while(l1/2*length1/n>R2_/2)
  {
    n=2*n;
  }

  double eta1_closestpoint=0.0;
  double eta2_closestpoint=0.0;

  for(int i=0;i<n+1;i++)
  {
    double eta1_slave=eta_left1+i*l1/n;
    double eta2_segleft=eta_left2;
    double eta2_master = 0.0;
    bool pairactive=false;
    double gap=0.0;
    double alpha=0.0;

    bool solutionwithinsegment=PointToLineProjection(eta1_slave,eta2_segleft,l2,eta2_master,gap,alpha,pairactive,true);

    if(solutionwithinsegment)
    {
      pointtolinesolfound=true;
      if(gap<g_min)
      {
        g_min=gap;
        alpha_g_min=alpha;
        eta1_closestpoint=eta1_slave;
        eta2_closestpoint=eta2_master;
      } //search also for the second-smallest gap
    }
  }

  bool cp_at_right_neighbor=false;

  //if we have a boundary minimum on the left, we also investigate the left neighbor point
  if(fabs(eta1_closestpoint-eta_left1)<1.0e-10)
  {
    double eta1_slave=eta_left1-l1/n;
    double eta2_segleft=eta_left2;
    double eta2_master = 0.0;
    bool pairactive=false;
    double gap=0.0;
    double alpha=0.0;

    bool solutionwithinsegment=PointToLineProjection(eta1_slave,eta2_segleft,l2,eta2_master,gap,alpha,pairactive,true);

    if(solutionwithinsegment)
    {
      if(gap<g_min)
      {
        cp_at_right_neighbor=true;
        g_min=gap;
        alpha_g_min=alpha;
        eta1_closestpoint=eta1_slave;
        eta2_closestpoint=eta2_master;
      } //search also for the second-smallest gap
    }
  }

  //if we have a boundary minimum on the right, we also investigate the right neighbor point
  if(fabs(eta1_closestpoint-(eta_left1+l1))<1.0e-10)
  {
    double eta1_slave=eta_left1+(n+1)*l1/n;
    double eta2_segleft=eta_left2;
    double eta2_master = 0.0;
    bool pairactive=false;
    double gap=0.0;
    double alpha=0.0;

    bool solutionwithinsegment=PointToLineProjection(eta1_slave,eta2_segleft,l2,eta2_master,gap,alpha,pairactive,true);

    if(solutionwithinsegment)
    {
      if(gap<g_min)
      {
        if(cp_at_right_neighbor)
          dserror("This should not happen, that we have a local minimum on the right and on the left neighbor!");

        g_min=gap;
        alpha_g_min=alpha;
        eta1_closestpoint=eta1_slave;
        eta2_closestpoint=eta2_master;
      } //search also for the second-smallest gap
    }
  }


  eta1_min=eta1_closestpoint;
  eta2_min=eta2_closestpoint;

  return;
}
/*------------------------------------------------------------------------------------------*
|  end: Determine minimal distance and contact angle for unconverged segment pair
 *------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateFcContact(Epetra_Vector* fint,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r1,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r2,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                                                                        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
                                                                        const double& intfac,
                                                                        bool cpp,
                                                                        bool gp,
                                                                        bool fixedendpointxi,
                                                                        bool fixedendpointeta,
                                                                        LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD,
                                                                        LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD)
{
  //Check for sensible combinations:
  if((cpp and (gp or fixedendpointxi or fixedendpointeta)) or (gp and (fixedendpointxi or fixedendpointeta)) )
    dserror("This is no possible combination of the parameters cpp, gp, fixedendpointxi and fixedendpointeta!");

  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // temporary vectors for contact forces, DOF-GIDs and owning procs
  LINALG::TMatrix<TYPE, dim1, 1> fc1(true);
  LINALG::TMatrix<TYPE, dim2, 1> fc2(true);
  Epetra_SerialDenseVector fcontact1(dim1);
  Epetra_SerialDenseVector fcontact2(dim2);

  //TODO: Introduce this quantities as class variables?
  std::vector<int>  lm1(dim1);
  std::vector<int>  lm2(dim2);
  std::vector<int>  lmowner1(dim1);
  std::vector<int>  lmowner2(dim2);

  // flag indicating assembly
  bool DoNotAssemble = true;

  // node ids of both elements
  const int* node_ids1 = element1_->NodeIds();
  const int* node_ids2 = element2_->NodeIds();

  for (int i=0;i<numnodes;++i)
  {
    // get node pointer and dof ids
    DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    // compute force vector Fc1 and prepare assembly
    for (int j=0;j<3*numnodalvalues;++j)
    {
      lm1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmowner1[3*numnodalvalues*i+j] = node->Owner();
    }
  }

  for (int i=0;i<numnodes;++i)
  {
    // get node pointer and dof ids
    DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    // compute force vector Fc1 and prepare assembly
    for (int j=0;j<3*numnodalvalues;++j)
    {
      lm2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmowner2[3*numnodalvalues*i+j] = node->Owner();
    }
  }

  TYPE gap = variables->GetGap();
  LINALG::TMatrix<TYPE,3,1> normal = variables->GetNormal();
  TYPE fp = variables->Getfp();
  //The factor ppfac reduces the penalty parameter for the large-angle and small-angle formulation
  //in dependence of the current contact angle
  TYPE ppfac = variables->GetPPfac();

  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (CheckContactStatus(FADUTILS::CastToDouble(gap)))
  {
    DoNotAssemble=false;
    #ifndef CONSISTENTTRANSITION
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    //The variable intfac represents the integration factor containing the
    //Gauss weight and the jacobian. This factor is only necessary for the
    //small-angle formulation and is set to 1.0 otherwise!
    for (int i=0;i<dim1;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc1(i) +=  N1(j,i)*normal(j)*fp*ppfac*intfac;
      }
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    //The variable intfac represents the integration factor containing the
    //Gauss weight and the jacobian. This factor is only necessary for the
    //small-angle formulation and is set to 1.0 otherwise!
    for (int i=0;i<dim2;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc2(i) +=  -N2(j,i)*normal(j)*fp*ppfac*intfac;
      }
    }
    #else
    // initialize storage for linearizations
    LINALG::TMatrix<TYPE, dim1+dim2, 1>  delta_xi(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_eta(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_coscontactangle(true);

    LINALG::TMatrix<TYPE, 3, 1> delta_r = FADUTILS::DiffVector(r1,r2);
    TYPE norm_delta_r = FADUTILS::VectorNorm<3>(delta_r);
    LINALG::TMatrix<TYPE, 3, 1> normal = variables->GetNormal();
    TYPE fp = variables->Getfp();
    TYPE dfp = variables->Getdfp();
    TYPE dppfac = variables->GetDPPfac();
    TYPE e = variables->GetEnergy();

    // linearization of contact point
    if(cpp)  //in case of large-angle-contact (standard closest-point-projection), we need delta_xi and delta_eta.
    {
      ComputeLinXiAndLinEta(delta_xi,delta_eta,delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);
    }
    else if(gp or (fixedendpointxi and !fixedendpointeta)) //in case of small-angle-contact (xi remains fixed), we only need delta_eta, delta_xi remains zero (this does not hold in case of ENDPOINTSEGMENTATION)
    {                                                      //this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint eta not!
      ComputeLinEtaFixXi(delta_eta,delta_r,r2_xi,r2_xixi,N1,N2,N2_xi);
      #ifdef ENDPOINTSEGMENTATION
        dserror("The combination of ENDPOINTSEGMENTATION and CONSISTENTTRANSITION is not possible!");
      #endif
    }
    else if(fixedendpointeta and !fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi not...
    {
      ComputeLinXiFixEta(delta_xi,delta_r,r1_xi,r1_xixi,N2,N1,N1_xi);
    }
    else if(fixedendpointeta and fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi is fixed...
    {
      //..we need to do nothing since delta_xi and delta_eta have to remain zero!
    }

    // linearization of large-angle/small-angle scale factor
    ComputeLinCosContactAngle(delta_coscontactangle,delta_xi,delta_eta,r1_xi,r2_xi,r1_xixi,r2_xixi,N1_xi,N2_xi);
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    //The variable intfac represents the integration factor containing the
    //Gauss weight and the jacobian. This factor is only necessary for the
    //small-angle formulation and is set to 1.0 otherwise!
    for (int i=0;i<dim1;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc1(i) +=  (N1(j,i)*normal(j)*fp*ppfac)*intfac;
      }
      fc1(i) +=  (e*dppfac*delta_coscontactangle(i))*intfac;
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    //The variable intfac represents the integration factor containing the
    //Gauss weight and the jacobian. This factor is only necessary for the
    //small-angle formulation and is set to 1.0 otherwise!
    for (int i=0;i<dim2;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc2(i) +=  (-N2(j,i)*normal(j)*fp*ppfac)*intfac;
      }
      fc2(i) +=  (e*dppfac*delta_coscontactangle(dim1+i))*intfac;
    }
    #endif
  }

  //Quantities necessary for automatic differentiation
  #ifdef AUTOMATICDIFF
  if (fc1_FAD != NULL and fc2_FAD != NULL)
  {
    for (int i=0;i<dim1;++i)
    {
        (*fc1_FAD)(i) = fc1(i);
    }
    for (int i=0;i<dim2;++i)
    {
        (*fc2_FAD)(i) = fc2(i);
    }
  }
  #endif

  //**********************************************************************
  // assemble contact forces
  //**********************************************************************
  if (!DoNotAssemble and fint != NULL)
  {
    for (int i=0;i<dim1;++i)
    {
        fcontact1[i] = FADUTILS::CastToDouble(fc1(i));
    }
    for (int i=0;i<dim2;++i)
    {
        fcontact2[i] = FADUTILS::CastToDouble(fc2(i));
    }
    // assemble fc1 and fc2 into global contact force vector
    LINALG::Assemble(*fint,fcontact1,lm1,lmowner1);
    LINALG::Assemble(*fint,fcontact2,lm2,lmowner2);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute contact forces
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate contact stiffness                               meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateStiffcContact( LINALG::SparseMatrix& stiffmatrix,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r1,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r2,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                             const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
                                                                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi,
                                                                             Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
                                                                             const double& intfac,
                                                                             bool cpp,
                                                                             bool gp,
                                                                             bool fixedendpointxi,
                                                                             bool fixedendpointeta)
{
  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  LINALG::TMatrix<TYPE, dim1, dim1+dim2> stiffc1(true);
  LINALG::TMatrix<TYPE, dim2, dim1+dim2> stiffc2(true);
  LINALG::TMatrix<TYPE, dim1, dim1+dim2> stiffc1_FAD(true);
  LINALG::TMatrix<TYPE, dim2, dim1+dim2> stiffc2_FAD(true);
  Epetra_SerialDenseMatrix stiffcontact1(dim1,dim1+dim2);
  Epetra_SerialDenseMatrix stiffcontact2(dim2,dim1+dim2);
  std::vector<int>  lmrow1(dim1);
  std::vector<int>  lmrow2(dim2);
  std::vector<int>  lmrowowner1(dim1);
  std::vector<int>  lmrowowner2(dim2);
  std::vector<int>  lmcol1(dim1+dim2);
  std::vector<int>  lmcol2(dim1+dim2);

  // flag indicating assembly
  bool DoNotAssemble = true;
  TYPE gap = variables->GetGap();
  //The factor ppfac reduces the penalty parameter for the large-angle and small-angle formulation
  //in dependence of the current contact angle
  TYPE ppfac = variables->GetPPfac();
  TYPE dppfac = variables->GetDPPfac();

  //In order to accelerate convergence, we only apply the basic stiffness part in case of very large gaps!
  double basicstiffgap = bcparams_.get<double>("BEAMS_BASICSTIFFGAP",-1.0);
  bool completestiff = true;
  if (basicstiffgap!=-1.0)
  {
    if (basicstiffgap<0.0)
      dserror("The parameter BEAMS_BASICSTIFFGAP has to be positive!");
    else if (gap < -1.0*basicstiffgap)
    {
      completestiff=false;
    }
  }

  //Apply additional weighting of the basic stiffness term e.g. in the first iterations or when
  //the Newton scheme oscillates (no convergence after a certain number of iterations)
  double basicstiffweightfac = 1.0;
  #ifdef BASICSTIFFWEIGHT
  if (iter_<5)
  {
    basicstiffweightfac=BASICSTIFFWEIGHT;
  }
  #endif

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  if (CheckContactStatus(FADUTILS::CastToDouble(gap)))
  {
    DoNotAssemble = false;

    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();

    //TODO: Introduce this quantities as class variables?
    //********************************************************************
    // prepare assembly
    //********************************************************************
    // fill lmrow1 and lmrowowner1
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
        lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
      }
    }

    // fill lmrow2 and lmrowowner2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmrow2[3*numnodalvalues*i+j]=NodeDofGIDs[j];
        lmrowowner2[3*numnodalvalues*i+j]=node->Owner();
      }
    }

    // fill lmcol1 and lmcol2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      }
    }

    // fill lmcol1 and lmcol2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmcol1[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmcol2[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
      }
    }

    // initialize storage for linearizations
    LINALG::TMatrix<TYPE, dim1+dim2, 1>  delta_xi(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_eta(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_gap(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_gap_t(true);
    LINALG::TMatrix<TYPE, 3, dim1+dim2> delta_x1_minus_x2(true);
    LINALG::TMatrix<TYPE, 3, dim1+dim2> delta_n(true);
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_coscontactangle(true);

    LINALG::TMatrix<TYPE, 3, 1> delta_r = FADUTILS::DiffVector(r1,r2);
    TYPE norm_delta_r = FADUTILS::VectorNorm<3>(delta_r);
    LINALG::TMatrix<TYPE, 3, 1> normal = variables->GetNormal();
    TYPE fp = variables->Getfp();
    TYPE dfp = variables->Getdfp();

    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************

    // linearization of contact point
    if(cpp)  //in case of large-angle-contact (standard closest-point-projection), we need delta_xi and delta_eta.
    {
      ComputeLinXiAndLinEta(delta_xi,delta_eta,delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);
    }
    else if(gp or (fixedendpointxi and !fixedendpointeta)) //in case of small-angle-contact (xi remains fixed), we only need delta_eta, delta_xi remains zero (this does not hold in case of ENDPOINTSEGMENTATION)
    {                                                      //this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint eta not!
      ComputeLinEtaFixXi(delta_eta,delta_r,r2_xi,r2_xixi,N1,N2,N2_xi);
    }
    else if(fixedendpointeta and !fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi not...
    {
      ComputeLinXiFixEta(delta_xi,delta_r,r1_xi,r1_xixi,N2,N1,N1_xi);
    }
    else if(fixedendpointeta and fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi is fixed...
    {
      //..we need to do nothing since delta_xi and delta_eta have to remain zero!
    }

    // linearization of gap function which is equal to delta d
    ComputeLinGap(delta_gap,delta_xi,delta_eta,delta_r,norm_delta_r,r1_xi,r2_xi,N1,N2);

    // linearization of normal vector
    ComputeLinNormal(delta_n,delta_xi,delta_eta,delta_r,r1_xi,r2_xi,N1,N2);

    // linearization of large-angle/small-angle scale factor
    ComputeLinCosContactAngle(delta_coscontactangle,delta_xi,delta_eta,r1_xi,r2_xi,r1_xixi,r2_xixi,N1_xi,N2_xi);

    #ifdef FADCHECKS
    std::cout << "delta_xi: " << std::endl;
      for (int i=0;i<dim1+dim2;i++)
        std::cout << delta_xi(i).val() << "  ";
      std::cout << std::endl << "delta_eta: " << std::endl;
      for (int i=0;i<dim1+dim2;i++)
        std::cout << delta_eta(i).val() << "  ";
      std::cout << std::endl;
      FADCheckLinXiAndLinEta(delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);
    #endif

    //*************Begin of standard linearization of penalty contact forces**************
    //The full contact stiffness is only applied if the contact flag is true
    //and gap_ > -BEAMS_BASICSTIFFGAP. If gap_ < -BEAMS_BASICSTIFFGAP, only
    //the basic stiffness is applied.

    //********************************************************************
    // evaluate contact stiffness
    // (1) stiffc1 of first element
    //********************************************************************

    //********************************************************************
    // part I - basic stiffness
    //********************************************************************

    LINALG::TMatrix<TYPE,dim1,1> N1T_normal(true);
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<dim1;j++)
      {
        N1T_normal(j)+=N1(i,j)*normal(i);
      }
    }
    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        stiffc1(i,j) += basicstiffweightfac * N1T_normal(i) * ( ppfac * dfp * delta_gap(j) + dppfac * fp * delta_coscontactangle(j) );
      }
    }

    //The geometric part is only applied for gap_ < -BEAMS_BASICSTIFFGAP
    if (completestiff)
    {
      //********************************************************************
      // part II - geometric stiffness 1
      //********************************************************************
      for  (int i=0;i<3;i++)
      {
        for (int j=0;j<dim1;j++)
        {
          for (int k=0;k<dim1+dim2;k++)
          {
              stiffc1(j,k) += ppfac*fp*N1(i,j)*delta_n(i,k);
          }
        }
      }
      //********************************************************************
      // part III - geometric stiffness 2
      //********************************************************************
      LINALG::TMatrix<TYPE,dim1,1> N1xiT_normal(true);
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<dim1;j++)
        {
          N1xiT_normal(j) += N1_xi(i,j)*normal(i);
        }
      }

      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          stiffc1(i,j) += ppfac*fp*N1xiT_normal(i)*delta_xi(j);
        }
      }
    }
    //********************************************************************
    // evaluate contact stiffness
    // (2) stiffc2 of second element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************
    LINALG::TMatrix<TYPE,dim2,1> N2T_normal(true);
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<dim2;j++)
      {
        N2T_normal(j)+=N2(i,j)*normal(i);
      }
    }
    for (int i=0;i<dim2;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        stiffc2(i,j) += -basicstiffweightfac * N2T_normal(i) * ( ppfac * dfp * delta_gap(j) + dppfac * fp * delta_coscontactangle(j) );


      }
    }

    if (completestiff)
    {
      //********************************************************************
      // part II
      //********************************************************************
      for  (int i=0;i<3;i++)
      {
        for (int j=0;j<dim2;j++)
        {
          for (int k=0;k<dim1+dim2;k++)
          {
              stiffc2(j,k) += -ppfac*fp*N2(i,j)*delta_n(i,k);
          }
        }
      }
      //********************************************************************
      // part III
      //********************************************************************
      LINALG::TMatrix<TYPE,dim1,1> N2xiT_normal(true);
      for (int i=0;i<3;i++)
      {
        for (int j=0;j<dim2;j++)
        {
          N2xiT_normal(j) += N2_xi(i,j)*normal(i);
        }
      }

      for (int i=0;i<dim2;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          stiffc2(i,j) += -ppfac*fp*N2xiT_normal(i)*delta_eta(j);
        }
      }
    }
    //*************End of standard linearization of penalty contact forces****************

    stiffc1.Scale(intfac);
    stiffc2.Scale(intfac);

    // automatic differentiation for debugging
    #ifdef AUTOMATICDIFF
      LINALG::TMatrix<TYPE, dim1, 1> fc1_FAD(true);
      LINALG::TMatrix<TYPE, dim2, 1> fc2_FAD(true);
      EvaluateFcContact(NULL, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, variables, intfac, cpp, gp, fixedendpointxi, fixedendpointeta, &fc1_FAD, &fc2_FAD);

      if(cpp)  //in case of large-angle-contact (standard closest-point-projection), we need delta_xi and delta_eta.
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          for (int i=0;i<dim1;i++)
            stiffc1_FAD(i,j) = (fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc1_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
          for (int i=0;i<dim2;i++)
            stiffc2_FAD(i,j) = (fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc2_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
        }
      }
      else if(gp or (fixedendpointxi and !fixedendpointeta)) //in case of small-angle-contact (xi remains fixed), we only need delta_eta, delta_xi remains zero (this does not hold in case of ENDPOINTSEGMENTATION)
      {                                                      //this also holds in case of ENDPOINTPENALTY when the endpoint xi is fixed and the endpoint eta not!
        for (int j=0;j<dim1+dim2;j++)
        {
          for (int i=0;i<dim1;i++)
            stiffc1_FAD(i,j) = (fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
          for (int i=0;i<dim2;i++)
            stiffc2_FAD(i,j) = (fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
        }
      }
      else if(fixedendpointeta and !fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi not...
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          for (int i=0;i<dim1;i++)
            stiffc1_FAD(i,j) = (fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2)*delta_xi(j));
          for (int i=0;i<dim2;i++)
            stiffc2_FAD(i,j) = (fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2)*delta_xi(j));
        }
      }
      else if(fixedendpointeta and fixedendpointxi) //In case of ENDPOINTPENALTY when the endpoint eta is fixed and the endpoint xi is fixed...
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          for (int i=0;i<dim1;i++)
            stiffc1_FAD(i,j) = (fc1_FAD(i).dx(j));
          for (int i=0;i<dim2;i++)
            stiffc2_FAD(i,j) = (fc2_FAD(i).dx(j));
        }
      }

//      std::cout << "Pair: " << element1_->Id() << " / " << element2_->Id() << std::endl;

//      std::cout << "stiffc1: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
//          std::cout << stiffc1(i,j).val() << " ";
          if(fabs(stiffc1(i,j).val())>1.0e-7 and fabs((stiffc1(i,j).val()-stiffc1_FAD(i,j).val())/stiffc1(i,j).val())>1.0e-7)
          {
            //std::cout << std::endl << std::endl << "stiffc1(i,j).val(): " << stiffc1(i,j).val() << "   stiffc1_FAD(i,j).val(): " << stiffc1_FAD(i,j).val() << std::endl;
            //dserror("Error in linearization!");
          }
        }
//        std::cout << std::endl;
      }
//      std::cout << std::endl;
//      std::cout << "stiffc1_FAD: " << std::endl;
//      for (int i=0;i<dim1;i++)
//      {
//        for (int j=0;j<dim1+dim2;j++)
//        {
//          std::cout << stiffc1_FAD(i,j).val() << " ";
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
//      std::cout << "stiffc2: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
//          std::cout << stiffc2(i,j).val() << " ";
          if(fabs(stiffc2(i,j).val())>1.0e-7 and fabs((stiffc2(i,j).val()-stiffc2_FAD(i,j).val())/stiffc2(i,j).val())>1.0e-7)
          {
            //std::cout << std::endl << std::endl <<"stiffc2(i,j).val(): " << stiffc2(i,j).val() << "   stiffc2_FAD(i,j).val(): " << stiffc2_FAD(i,j).val() << std::endl;
            //dserror("Error in linearization!");
          }
        }
//        std::cout << std::endl;
      }
//      std::cout << std::endl;
//      std::cout << "stiffc2_FAD: " << std::endl;
//      for (int i=0;i<dim1;i++)
//      {
//        for (int j=0;j<dim1+dim2;j++)
//        {
//          std::cout << stiffc2_FAD(i,j).val() << " ";
//        }
//        std::cout << std::endl;
//      }
//      std::cout << std::endl;
    #endif
  }//if (CheckContactStatus(gap))

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  // change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.

  // now finally assemble stiffc1 and stiffc2
  if (!DoNotAssemble)
  {
    #ifndef AUTOMATICDIFF
      for (int j=0;j<dim1+dim2;j++)
      {
        for (int i=0;i<dim1;i++)
          stiffcontact1(i,j) = -FADUTILS::CastToDouble(stiffc1(i,j));
        for (int i=0;i<dim2;i++)
          stiffcontact2(i,j) = -FADUTILS::CastToDouble(stiffc2(i,j));
      }
    #else
      for (int j=0;j<dim1+dim2;j++)
      {
        for (int i=0;i<dim1;i++)
          stiffcontact1(i,j) = -FADUTILS::CastToDouble(stiffc1_FAD(i,j));
        for (int i=0;i<dim2;i++)
          stiffcontact2(i,j) = -FADUTILS::CastToDouble(stiffc2_FAD(i,j));
      }
    #endif

    stiffmatrix.Assemble(0,stiffcontact1,lmrow1,lmrowowner1,lmcol1);
    stiffmatrix.Assemble(0,stiffcontact2,lmrow2,lmrowowner2,lmcol2);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate contact stiffness
 *----------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------*
 |  FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION    meier 10/14|
 *------------------------------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateStiffcContactIntSeg( LINALG::SparseMatrix& stiffmatrix,
                                                                                   const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi_bound,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r1,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r2,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                                   const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                                   const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                                   const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                                   const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                                   const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                                                                                   Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables,
                                                                                   const double& intfac,
                                                                                   const double& d_xi_ele_d_xi_bound,
                                                                                   TYPE signed_jacobi_interval)
{
  #ifndef AUTOMATICDIFF
    dserror("This method only works with automatic differentiation!");
  #endif

  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  LINALG::TMatrix<TYPE, dim1, dim1+dim2> stiffc1_FAD(true);
  LINALG::TMatrix<TYPE, dim2, dim1+dim2> stiffc2_FAD(true);
  Epetra_SerialDenseMatrix stiffcontact1(dim1,dim1+dim2);
  Epetra_SerialDenseMatrix stiffcontact2(dim2,dim1+dim2);
  std::vector<int>  lmrow1(dim1);
  std::vector<int>  lmrow2(dim2);
  std::vector<int>  lmrowowner1(dim1);
  std::vector<int>  lmrowowner2(dim2);
  std::vector<int>  lmcol1(dim1+dim2);
  std::vector<int>  lmcol2(dim1+dim2);

  // flag indicating assembly
  bool DoNotAssemble = true;
  TYPE gap = cpvariables->GetGap();

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  if (CheckContactStatus(FADUTILS::CastToDouble(gap)))
  {
    DoNotAssemble = false;

    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();

    //TODO: Introduce this quantities as class variables?
    //********************************************************************
    // prepare assembly
    //********************************************************************
    // fill lmrow1 and lmrowowner1
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
        lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
      }
    }

    // fill lmrow2 and lmrowowner2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmrow2[3*numnodalvalues*i+j]=NodeDofGIDs[j];
        lmrowowner2[3*numnodalvalues*i+j]=node->Owner();
      }
    }

    // fill lmcol1 and lmcol2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      }
    }

    // fill lmcol1 and lmcol2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmcol1[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmcol2[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
      }
    }

    // initialize storage for linearizations
    LINALG::TMatrix<TYPE, dim1+dim2, 1> delta_eta(true);
    LINALG::TMatrix<TYPE, 3, 1> delta_r = FADUTILS::DiffVector(r1,r2);

    ComputeLinEtaFixXi(delta_eta,delta_r,r2_xi,r2_xixi,N1,N2,N2_xi);

    LINALG::TMatrix<TYPE, dim1, 1> fc1_FAD(true);
    LINALG::TMatrix<TYPE, dim2, 1> fc2_FAD(true);
    EvaluateFcContact(NULL, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, cpvariables, intfac, false, true, false, false, &fc1_FAD, &fc2_FAD);

    #ifdef AUTOMATICDIFF
    TYPE fac1(0.0);
    TYPE fac2(0.0);
    fac2=-FADUTILS::ScalarProduct(r2_xi, r2_xi) + FADUTILS::ScalarProduct(delta_r, r2_xixi);
    fac1=FADUTILS::ScalarProduct(r2_xi, r1_xi);
    for (int j=0;j<dim1+dim2;j++)
    {
      for (int i=0;i<dim1;i++)
        stiffc1_FAD(i,j) = fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2)*d_xi_ele_d_xi_bound*delta_xi_bound(j) +fc1_FAD(i).dx(dim1+dim2+1)*(delta_eta(j)-fac1/fac2*d_xi_ele_d_xi_bound*delta_xi_bound(j))+fc1_FAD(i).val()/(2.0*signed_jacobi_interval)*delta_xi_bound(j);
                        // d(f)/d(disp)    +d(f)/d(xi,GP)           *d(xi,GP)/d(disp)                      +d(f)/d(eta,GP)            *d(eta,GP)/d(disp)                                             +d(f)/d(xi,Bound)*d(xi,Bound)/d(disp)

      for (int i=0;i<dim2;i++)
        stiffc2_FAD(i,j) = fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2)*d_xi_ele_d_xi_bound*delta_xi_bound(j) +fc2_FAD(i).dx(dim1+dim2+1)*(delta_eta(j)-fac1/fac2*d_xi_ele_d_xi_bound*delta_xi_bound(j))+fc2_FAD(i).val()/(2.0*signed_jacobi_interval)*delta_xi_bound(j);
    }
    #endif

  }//if (CheckContactStatus(gap))

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  // change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.

  // now finally assemble stiffc1 and stiffc2
  if (!DoNotAssemble)
  {
    #ifndef AUTOMATICDIFF
      dserror("This method only works with AUTOMATICDIFF");
    #else
      for (int j=0;j<dim1+dim2;j++)
      {
        for (int i=0;i<dim1;i++)
          stiffcontact1(i,j) = -FADUTILS::CastToDouble(stiffc1_FAD(i,j));
        for (int i=0;i<dim2;i++)
          stiffcontact2(i,j) = -FADUTILS::CastToDouble(stiffc2_FAD(i,j));
      }
    #endif

    stiffmatrix.Assemble(0,stiffcontact1,lmrow1,lmrowowner1,lmcol1);
    stiffmatrix.Assemble(0,stiffcontact2,lmrow2,lmrowowner2,lmcol2);
  }

  return;
}
/*------------------------------------------------------------------------------------------*
 |  end: FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION
 *------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Linearizations of contact point                          meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinXiAndLinEta(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                            LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
{
  //**********************************************************************
  // we have to solve the following system of equations:
  //  _              _       _      _       _              _      _       _
  // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
  // |                | *  |         | =  |                 | * |        |
  // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
  //
  // this can be done easily because it is a linear 2x2-system.
  // we obtain the solution by inverting matrix L:
  //
  // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
  //
  //**********************************************************************

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  LINALG::TMatrix<TYPE,2,2> L(true);
  LINALG::TMatrix<TYPE,2,2> L_inv(true);
  LINALG::TMatrix<TYPE,2,dim1+dim2> B(true);
  LINALG::TMatrix<TYPE,2,dim1+dim2> D(true);

  // compute L elementwise
  L(0,0)=::FADUTILS::ScalarProduct(r1_xi, r1_xi) + ::FADUTILS::ScalarProduct(delta_r, r1_xixi);
  L(1,1)=-::FADUTILS::ScalarProduct(r2_xi, r2_xi) + ::FADUTILS::ScalarProduct(delta_r, r2_xixi);
  L(0,1)=-::FADUTILS::ScalarProduct(r2_xi, r1_xi);
  L(1,0)=-L(0,1);

  // invert L by hand
  TYPE det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
  if (FADUTILS::CastToDouble(FADUTILS::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
  L_inv(0,0) =  L(1,1) / det_L;
  L_inv(0,1) = -L(0,1) / det_L;
  L_inv(1,0) = -L(1,0) / det_L;
  L_inv(1,1) =  L(0,0) / det_L;

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1;j++)
    {
      B(0,j)+= -delta_r(i)*N1_xi(i,j) - r1_xi(i)*N1(i,j);
      B(1,j)+= - r2_xi(i)*N1(i,j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim2;j++)
    {
      B(0,j+dim1)+= r1_xi(i)*N2(i,j);
      B(1,j+dim1)+= -delta_r(i)*N2_xi(i,j) + r2_xi(i)*N2(i,j);
    }
  }

  // compute D = L^-1 * B
  D.Multiply(L_inv, B);

  // finally the linearizations / directional derivatives
  for (int i=0;i<dim1+dim2;i++)
  {
    delta_xi(i) = D(0,i);
    delta_eta(i) = D(1,i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Linearizations of contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Lin. of contact point coordinate eta with fixed xi        meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinEtaFixXi( LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                            const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                            const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
{

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  TYPE L=0.0;
  LINALG::TMatrix<TYPE,1,dim1+dim2> B(true);

  // compute L elementwise
  L=-FADUTILS::ScalarProduct(r2_xi, r2_xi) + FADUTILS::ScalarProduct(delta_r, r2_xixi);

//  std::cout << "r2_xi: " << r2_xi << std::endl;
//  std::cout << "r2_xixi: " << r2_xixi << std::endl;

  if (fabs(FADUTILS::CastToDouble(L)) < COLINEARTOL)
      dserror("Linearization of point to line projection is zero, choose tighter search boxes!");

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1;j++)
    {
      B(0,j)+= - r2_xi(i)*N1(i,j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim2;j++)
    {
      B(0,j+dim1)+= -delta_r(i)*N2_xi(i,j) + r2_xi(i)*N2(i,j);
    }
  }

  // finally the linearizations / directional derivatives
  for (int i=0;i<dim1+dim2;i++)
  {
    delta_eta(i) = 1.0/L*B(0,i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Lin. of contact point coordinate eta with fixed xi
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Lin. of contact point coordinate xi with fixed eta        meier 12/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinXiFixEta( LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi,
                                                                          const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi)
{

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // matrices to compute Lin_Xi and Lin_Eta
  TYPE L=0.0;
  LINALG::TMatrix<TYPE,1,dim1+dim2> B(true);

  // compute L elementwise
  L=FADUTILS::ScalarProduct(r1_xi, r1_xi) + FADUTILS::ScalarProduct(delta_r, r1_xixi);

  if (fabs(FADUTILS::CastToDouble(L)) < COLINEARTOL)
      dserror("Linearization of point to line projection is zero, choose tighter search boxes!");

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1;j++)
    {
      B(0,j)+= -delta_r(i)*N1_xi(i,j) - r1_xi(i)*N1(i,j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim2;j++)
    {
      B(0,j+dim1)+= r1_xi(i)*N2(i,j);
    }
  }

  // finally the linearizations / directional derivatives
  for (int i=0;i<dim1+dim2;i++)
  {
    delta_xi(i) = 1.0/L*B(0,i);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Lin. of contact point coordinate xi with fixed eta
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of integration interval bounds      meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinXiBound(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi_bound,
                                                                        TYPE& eta1_bound,
                                                                        TYPE eta2)
{

    // vectors for shape functions and their derivatives
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
    LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta

    // coords and derivatives of the two contacting points
    LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
    LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
    LINALG::TMatrix<TYPE, 3, 1> r1_xi(true);                            // = r1,xi
    LINALG::TMatrix<TYPE, 3, 1> r2_xi(true);                            // = r2,eta
    LINALG::TMatrix<TYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
    LINALG::TMatrix<TYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
    LINALG::TMatrix<TYPE, 3, 1> delta_r(true);                          // = r1-r2

    // update shape functions and their derivatives
    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1_bound, eta2);
    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    delta_r=FADUTILS::DiffVector(r1,r2);

    const int dim1 = 3*numnodes*numnodalvalues;
    const int dim2 = 3*numnodes*numnodalvalues;

    // matrices to compute Lin_Xi and Lin_Eta
    TYPE a_11(0.0);
    LINALG::TMatrix<TYPE,2,dim1+dim2> B(true);

    a_11=FADUTILS::ScalarProduct(r1_xi, r1_xi) + FADUTILS::ScalarProduct(delta_r, r1_xixi);

    #ifdef CHANGEENDPOINTPROJECTION
    TYPE a_21(0.0);
    a_21=FADUTILS::ScalarProduct(r1_xi, r2_xi);
    #endif


    for (int i=0;i<3;i++)
    {
      for (int j=0;j<dim1;j++)
      {
        B(0,j)+= -delta_r(i)*N1_xi(i,j) - r1_xi(i)*N1(i,j);
        B(1,j)+= - r2_xi(i)*N1(i,j);
      }
    }

    for (int i=0;i<3;i++)
    {
      for (int j=0;j<dim2;j++)
      {
        B(0,j+dim1)+= r1_xi(i)*N2(i,j);
        B(1,j+dim1)+= -delta_r(i)*N2_xi(i,j) + r2_xi(i)*N2(i,j);
      }
    }

    #ifndef CHANGEENDPOINTPROJECTION
      // finally the linearizations / directional derivatives in case the orthogonality condition is fulfilled on beam1
      for (int i=0;i<dim1+dim2;i++)
      {
        delta_xi_bound(i) = B(0,i)/a_11;
      }
    #else
      // finally the linearizations / directional derivatives in case the orthogonality condition is fulfilled on beam2
      for (int i=0;i<dim1+dim2;i++)
      {
        delta_xi_bound(i) = B(1,i)/a_21;
      }
    #endif

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute linearization of integration interval bounds
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of gap                              meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinGap( LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_gap,
                                                                     const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                     const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                     const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                     const TYPE& norm_delta_r,
                                                                     const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                     const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                     const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                     const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // delta g := delta_r/||delta_r||*auxiliary_matri1 delta d, with auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  LINALG::TMatrix<TYPE,3,dim1+dim2>  auxiliary_matrix1(true);

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      auxiliary_matrix1(i,j)+=r1_xi(i)*delta_xi(j)-r2_xi(i)*delta_eta(j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1;j++)
    {
      auxiliary_matrix1(i,j)+= N1(i,j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim2;j++)
    {
      auxiliary_matrix1(i,j+dim1)+= -N2(i,j);
    }
  }

  // compute linearization of gap
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      delta_gap(j) +=  delta_r(i) * auxiliary_matrix1(i,j)/norm_delta_r;
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of gap
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of cosine of contact angle          meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinCosContactAngle( LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_coscontactangle,
                                                                                 LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi,
                                                                                 LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                                 const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                                 const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                                 const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                                 const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                                 const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                                 const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  TYPE norm_r1xi = FADUTILS::VectorNorm<3>(r1_xi);
  TYPE norm_r2xi = FADUTILS::VectorNorm<3>(r2_xi);
  LINALG::TMatrix<TYPE, 3, 1> r1_xi_unit(r1_xi);
  LINALG::TMatrix<TYPE, 3, 1> r2_xi_unit(r2_xi);
  r1_xi_unit.Scale(1.0/norm_r1xi);
  r2_xi_unit.Scale(1.0/norm_r2xi);
  TYPE r1xi_unit_r2xi_unit = FADUTILS::ScalarProduct(r1_xi_unit,r2_xi_unit);

  //Pre-factor representing the modulus, since s=|r1xi_unit_r2xi_unit|
  double modulus_factor = 1.0;

  if(r1xi_unit_r2xi_unit<0.0)
    modulus_factor=-1.0;

  LINALG::TMatrix<TYPE, 3, 1> v1(r2_xi_unit);
  LINALG::TMatrix<TYPE, 3, 1> v2(r1_xi_unit);
  v1.Update(-r1xi_unit_r2xi_unit,r1_xi_unit,1.0);
  v2.Update(-r1xi_unit_r2xi_unit,r2_xi_unit,1.0);
  v1.Scale(1.0/norm_r1xi);
  v2.Scale(1.0/norm_r2xi);

  LINALG::TMatrix<TYPE, 3, dim1+dim2> delta_r1_xi(true);
  LINALG::TMatrix<TYPE, 3, dim1+dim2> delta_r2_xi(true);

  for(int i=0;i<3;i++)
  {
    for(int j=0;j<dim1+dim2;j++)
    {
      delta_r1_xi(i,j)=r1_xixi(i)*delta_xi(j);

      if(j<dim1)
        delta_r1_xi(i,j)+=N1_xi(i,j);
    }
  }

  for(int i=0;i<3;i++)
  {
    for(int j=0;j<dim1+dim2;j++)
    {
      delta_r2_xi(i,j)=r2_xixi(i)*delta_eta(j);

      if(j>=dim1)
        delta_r2_xi(i,j)+=N2_xi(i,j-dim1);
    }
  }

  LINALG::TMatrix<TYPE, 1, dim1+dim2> v1_delta_r1_xi(true);
  LINALG::TMatrix<TYPE, 1, dim1+dim2> v2_delta_r2_xi(true);
  v1_delta_r1_xi.MultiplyTN(v1,delta_r1_xi);
  v2_delta_r2_xi.MultiplyTN(v2,delta_r2_xi);

  for (int j=0;j<dim1+dim2;j++)
  {
    delta_coscontactangle(j)= modulus_factor*(v1_delta_r1_xi(j)+v2_delta_r2_xi(j));
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of cosine of contact angle
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of normal vector                    meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeLinNormal( LINALG::TMatrix<TYPE, 3, 2*3*numnodes*numnodalvalues>& delta_normal,
                                                                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                        const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  //delta n := auxiliary_matri2*auxiliary_matrix1* delta d, with auxiliary_matri2 = (I-nxn)/||r1-r2||
  //and auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  TYPE norm_delta_r = FADUTILS::VectorNorm<3>(delta_r);
  LINALG::TMatrix<TYPE, 3, 1> normal(delta_r);
  normal.Scale(1.0/norm_delta_r);

  LINALG::TMatrix<TYPE,3,dim1+dim2>  auxiliary_matrix1(true);
  LINALG::TMatrix<TYPE,3,3>  auxiliary_matrix2(true);

  //compute auxiliary_matrix1
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      auxiliary_matrix1(i,j)+=r1_xi(i)*delta_xi(j)-r2_xi(i)*delta_eta(j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1;j++)
    {
      auxiliary_matrix1(i,j)+= N1(i,j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim2;j++)
    {
      auxiliary_matrix1(i,j+dim1)+= -N2(i,j);
    }
  }

  //compute auxiliary_matrix2
  for (int i=0;i<3;i++)
  {
    auxiliary_matrix2(i,i)+= 1.0/norm_delta_r;
    for (int j=0;j<3;j++)
    {
      auxiliary_matrix2(i,j)+= -normal(i)*normal(j)/norm_delta_r;
    }
  }

  // compute linearization of normal vector
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<dim1+dim2;k++)
        delta_normal(i,k) +=  auxiliary_matrix2(i,j) * auxiliary_matrix1(j,k);

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of normal vector
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetShapeFunctions( LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                            LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                            LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                            LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                                                                            LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
                                                                            LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi,
                                                                            const TYPE& eta1,
                                                                            const TYPE& eta2)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = element1_->Shape();
  const DRT::Element::DiscretizationType distype2 = element2_->Shape();

  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N1_i(true);
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N1_i_xi(true);
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N1_i_xixi(true);
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N2_i(true);
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N2_i_xi(true);
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N2_i_xixi(true);

  if (numnodalvalues==1)
  {
    // get values and derivatives of shape functions
    DRT::UTILS::shape_function_1D(N1_i, eta1, distype1);
    DRT::UTILS::shape_function_1D(N2_i, eta2, distype2);
    DRT::UTILS::shape_function_1D_deriv1(N1_i_xi, eta1, distype1);
    DRT::UTILS::shape_function_1D_deriv1(N2_i_xi, eta2, distype2);
    DRT::UTILS::shape_function_1D_deriv2(N1_i_xixi, eta1, distype1);
    DRT::UTILS::shape_function_1D_deriv2(N2_i_xixi, eta2, distype2);
  }
  else if (numnodalvalues==2)
  {

    if ( element1_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    if ( element2_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    double length1 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();
    double length2 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element2_))->jacobi();

    // get values and derivatives of shape functions
    DRT::UTILS::shape_function_hermite_1D(N1_i, eta1, length1, distype1);
    DRT::UTILS::shape_function_hermite_1D(N2_i, eta2, length2, distype2);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi, eta1, length1, distype1);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N2_i_xi, eta2, length2, distype2);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N1_i_xixi, eta1, length1, distype1);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N2_i_xixi, eta2, length2, distype2);

  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  //Assemble the individual shape functions in matrices, such that: r1=N1*d1, r1_xi=N1_xi*d1, r1_xixi=N1_xixi*d1, r2=N2*d2, r2_xi=N2_xi*d2, r2_xixi=N2_xixi*d2
  AssembleShapefunctions(N1_i, N1_i_xi, N1_i_xixi, N1, N1_xi, N1_xixi);
  AssembleShapefunctions(N2_i, N2_i_xi, N2_i_xixi, N2, N2_xi, N2_xixi);

  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::GetShapeFunctions( LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N,
                                                                         const TYPE& eta,
                                                                         int deriv,
                                                                         DRT::Element* ele)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype = ele->Shape();
  LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues> N_i(true);

  if (numnodalvalues==1)
  {
    // get values and derivatives of shape functions
    switch (deriv)
    {
      case 0:
      {
        DRT::UTILS::shape_function_1D(N_i, eta, distype);
        break;
      }
      case 1:
      {
        DRT::UTILS::shape_function_1D_deriv1(N_i, eta, distype);
        break;
      }
      case 2:
      {
        DRT::UTILS::shape_function_1D_deriv2(N_i, eta, distype);
        break;
      }
    }
  }
  else if (numnodalvalues==2)
  {

    if ( ele->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    double length = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(ele))->jacobi();

    // get values and derivatives of shape functions
    switch (deriv)
    {
      case 0:
      {
        DRT::UTILS::shape_function_hermite_1D(N_i, eta, length, distype);
        break;
      }
      case 1:
      {
        DRT::UTILS::shape_function_hermite_1D_deriv1(N_i, eta, length, distype);
        break;
      }
      case 2:
      {
        DRT::UTILS::shape_function_hermite_1D_deriv2(N_i, eta, length, distype);
        break;
      }
    }
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  //Assemble the individual shape functions in matrices, such that: r1=N1*d1, r1_xi=N1_xi*d1, r1_xixi=N1_xixi*d1, r2=N2*d2, r2_xi=N2_xi*d2, r2_xixi=N2_xixi*d2
  AssembleShapefunctions(N_i, N);

  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble one shape function matrix                                                             meier 10/14|
 *-----------------------------------------------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::AssembleShapefunctions(const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i,
                                                                             LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N)
{
  //assembly_N is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N
  int assembly_N[3][3*numnodes*numnodalvalues];

  //Initialize to zero
  for (int i=0;i<3*numnodes*numnodalvalues;i++)
  {
    for (int j=0;j<3; j++)
    {
      assembly_N[j][i]=0.0;
    }
  }

  /*
  Set number of shape functions for each 3*3 block:
  e.g. second order Reissner beam (numnodes=3, numnodalvalues=1)
  int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
                           {0,1,0,0,2,0,0,3,0},
                           {0,0,1,0,0,2,0,0,3}};

  e.g. Kirchhoff beam (numnodes=2, numnodalvalues=2)
  int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
                           {0,1,0,0,2,0,0,3,0,0,4,0},
                           {0,0,1,0,0,2,0,0,3,0,0,4}};
  */

  for (int i=0;i<numnodes*numnodalvalues;i++)
  {
    assembly_N[0][3*i]=i+1;
    assembly_N[1][3*i+1]=i+1;
    assembly_N[2][3*i+2]=i+1;
  }

  //Assemble the matrices of the shape functions
  for (int i=0; i<3*numnodes*numnodalvalues; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N(j,i)=0;
      }
      else
      {
        N(j,i)=N_i(assembly_N[j][i]-1);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble all shape functions                                                                  meier 01/14|
 *-----------------------------------------------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::AssembleShapefunctions(const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i,
                                                                                const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i_xi,
                                                                                const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i_xixi,
                                                                                LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N,
                                                                                LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N_xi,
                                                                                LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N_xixi)
{
  //assembly_N is just an array to help assemble the matrices of the shape functions
  //it determines, which shape function is used in which column of N
  int assembly_N[3][3*numnodes*numnodalvalues];

  //Initialize to zero
  for (int i=0;i<3*numnodes*numnodalvalues;i++)
  {
    for (int j=0;j<3; j++)
    {
      assembly_N[j][i]=0.0;
    }
  }

  /*
  Set number of shape functions for each 3*3 block:
  e.g. second order Reissner beam (numnodes=3, numnodalvalues=1)
  int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
                           {0,1,0,0,2,0,0,3,0},
                           {0,0,1,0,0,2,0,0,3}};

  e.g. Kirchhoff beam (numnodes=2, numnodalvalues=2)
  int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
                           {0,1,0,0,2,0,0,3,0,0,4,0},
                           {0,0,1,0,0,2,0,0,3,0,0,4}};
  */

  for (int i=0;i<numnodes*numnodalvalues;i++)
  {
    assembly_N[0][3*i]=i+1;
    assembly_N[1][3*i+1]=i+1;
    assembly_N[2][3*i+2]=i+1;
  }

  //Assemble the matrices of the shape functions
  for (int i=0; i<3*numnodes*numnodalvalues; i++)
  {
    for (int j=0; j<3; j++)
    {
      if(assembly_N[j][i]==0)
      {
        N(j,i)=0;
        N_xi(j,i)=0;
        N_xixi(j,i)=0;
      }
      else
      {
        N(j,i)=N_i(assembly_N[j][i]-1);
        N_xi(j,i)=N_i_xi(assembly_N[j][i]-1);
        N_xixi(j,i)=N_i_xixi(assembly_N[j][i]-1);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | compute position at given curve point                  meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
LINALG::TMatrix<TYPE,3,1> CONTACT::Beam3contact<numnodes, numnodalvalues>::r(const TYPE& eta, DRT::Element* ele)
{

  LINALG::TMatrix<TYPE,3,1> r(true);
  LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues> N(true);
  GetShapeFunctions(N,eta,0,ele);

  if(ele->Id()==element1_->Id())
  {
    // compute output variable
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnodes*numnodalvalues;j++)
      {
        r(i)+=N(i,j)*ele1pos_(j);
      }
    }
  }
  else if(ele->Id()==element2_->Id())
  {
    // compute output variable
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnodes*numnodalvalues;j++)
      {
        r(i)+=N(i,j)*ele2pos_(j);
      }
    }
  }
  else
    dserror("This method can only applied to element1_ and element2_!");

  return r;
}
/*----------------------------------------------------------------------*
 | end: compute position at given curve point              meier 02/14|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute tangent at given curve point                  meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
LINALG::TMatrix<TYPE,3,1> CONTACT::Beam3contact<numnodes, numnodalvalues>::r_xi(const TYPE& eta, DRT::Element* ele)
{

  LINALG::TMatrix<TYPE,3,1> r_xi(true);
  LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues> N_xi(true);
  GetShapeFunctions(N_xi,eta,1,ele);

  if(ele->Id()==element1_->Id())
  {
    // compute output variable
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnodes*numnodalvalues;j++)
      {
        r_xi(i)+=N_xi(i,j)*ele1pos_(j);
      }
    }
  }
  else if(ele->Id()==element2_->Id())
  {
    // compute output variable
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnodes*numnodalvalues;j++)
      {
        r_xi(i)+=N_xi(i,j)*ele2pos_(j);
      }
    }
  }
  else
    dserror("This method can only applied to element1_ and element2_!");

  return r_xi;
}
/*----------------------------------------------------------------------*
 | end: compute tangent at given curve point              meier 02/14|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute contact point coordinates and their derivatives   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeCoordsAndDerivs(LINALG::TMatrix<TYPE,3,1>& r1,
                                                                                LINALG::TMatrix<TYPE,3,1>& r2,
                                                                                LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                                LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                                LINALG::TMatrix<TYPE,3,1>& r1_xixi,
                                                                                LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1_xi,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2_xi,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1_xixi,
                                                                                const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2_xixi)
{
  r1.Clear();
  r2.Clear();
  r1_xi.Clear();
  r2_xi.Clear();
  r1_xixi.Clear();
  r2_xixi.Clear();

#ifdef AUTOMATICDIFF
  BEAMCONTACT::SetFADDispDofs<numnodes, numnodalvalues>(ele1pos_,ele2pos_);
#endif

  // compute output variable
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<3*numnodes*numnodalvalues;j++)
    {
      r1(i)+=N1(i,j)*ele1pos_(j);
      r2(i)+=N2(i,j)*ele2pos_(j);
      r1_xi(i)+=N1_xi(i,j)*ele1pos_(j);
      r2_xi(i)+=N2_xi(i,j)*ele2pos_(j);
      r1_xixi(i)+=N1_xixi(i,j)*ele1pos_(j);
      r2_xixi(i)+=N2_xixi(i,j)*ele2pos_(j);
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | end: compute contact point coordinates and their derivatives         |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate function f in CPP                               meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateOrthogonalityCondition(LINALG::TMatrix<TYPE,2,1>& f,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                                                                        const double norm_delta_r,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& t1,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& t2)
{
  // reset f
  f.Clear();

  int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(bcparams_,"BEAMS_SMOOTHING");
  // evaluate f
  // see Wriggers, Computational Contact Mechanics, equation (12.5)
  if (smoothing == INPAR::BEAMCONTACT::bsm_none) //non-smoothed
  {
    for (int i=0;i<3;i++)
    {
      f(0) += delta_r(i)*r1_xi(i) / norm_delta_r;
      f(1) += -delta_r(i)*r2_xi(i) / norm_delta_r;
    }
  }
  else //smoothed
  {
    dserror("The smoothing procedure is not consistent linearized so far! Thereto, the quantities lin_xi and "
        "lin_eta have to be calculated consistent to the smoothed orthogonality condition below!");
    for (int i=0;i<3;i++)
    {
      f(0) += delta_r(i)*t1(i) / norm_delta_r;
      f(1) += -delta_r(i)*t2(i) / norm_delta_r;
    }
  }


  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate function f in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df in CPP                              meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateLinOrthogonalityCondition( LINALG::TMatrix<TYPE,2,2>& df,
                                                                                            LINALG::TMatrix<TYPE,2,2>& dfinv,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                                                                            const double norm_delta_r,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r1_xixi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& t1,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& t2,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& t1_xi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& t2_xi,
                                                                                            bool& elementscolinear)

{
  // reset df and dfinv
  df.Clear();
  dfinv.Clear();

  int smoothing = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(bcparams_,"BEAMS_SMOOTHING");

  // evaluate df
  // see Wriggers, Computational Contact Mechanics, equation (12.7)
  if (smoothing == INPAR::BEAMCONTACT::bsm_none) //non-smoothed
  {
    for(int i=0;i<3;i++)
    {
      df(0,0) += (r1_xi(i)*r1_xi(i) + delta_r(i)*r1_xixi(i)) / norm_delta_r;
      df(0,1) += -r1_xi(i)*r2_xi(i) / norm_delta_r;
      df(1,0) += -r2_xi(i)*r1_xi(i) / norm_delta_r;
      df(1,1) += (r2_xi(i)*r2_xi(i) - delta_r(i)*r2_xixi(i)) / norm_delta_r;
    }
  }
  else //smoothed
  {
    for(int i=0;i<3;i++)
    {
      df(0,0) += (r1_xi(i)*t1(i) + delta_r(i)*t1_xi(i)) / norm_delta_r;
      df(0,1) += -t1(i)*r2_xi(i) / norm_delta_r;
      df(1,0) += -t2(i)*t1_xi(i) / norm_delta_r;
      df(1,1) += (r2_xi(i)*t2(i) - delta_r(i)*t2_xi(i)) / norm_delta_r;
    }
  }

  // Inverting (2x2) matrix df by hard coded formula, so that it is
  // possible to handle colinear vectors, because they lead to det(df) =0
  TYPE det_df = df(0,0)*df(1,1)-df(1,0)*df(0,1);

  //********************************************************************
  // ASSUMPTION:
  // If det_df=0 we assume, that the two elements have an identical
  // neutral axis. These contact objects will be rejected. The outcome
  // of this physically rare phenomenon is that handling of line contact
  // is not possible with this approach.
  //********************************************************************

  // singular df
  if (FADUTILS::CastToDouble(FADUTILS::Norm(det_df))<COLINEARTOL)
  {
    // sort out
    elementscolinear = true;
  }
  // regular df (inversion possible)
  else
  {
    // do not sort out
    elementscolinear = false;

    // invert df
    dfinv(0,0)=df(1,1)/det_df;
    dfinv(0,1)=-df(0,1)/det_df;
    dfinv(1,0)=-df(1,0)/det_df;
    dfinv(1,1)=df(0,0)/det_df;
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate orthogonality cond. of point to line projeciton  meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluatePTLOrthogonalityCondition(TYPE& f,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                                                                        const double norm_delta_r,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                                        const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                                        bool orthogonalprojection)
{
  // reset f
  f=0.0;

  // evaluate f
  if(orthogonalprojection==false)//standard case
  {
    for (int i=0;i<3;i++)
    {
      f += -delta_r(i)*r2_xi(i) / norm_delta_r;
    }
  }
  else
  {
    for (int i=0;i<3;i++)
    {
      f += -delta_r(i)*r1_xi(i) / norm_delta_r;
    }
  }


  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate orthogonality cond. of point to line projeciton
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df of PTLOrthogonalityCondition        meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateLinPTLOrthogonalityCondition( TYPE& df,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                                                                            const double norm_delta_r,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                                                                            const LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                                                                                            bool orthogonalprojection)

{
  // reset df
  df=0.0;

  // evaluate df
  if(orthogonalprojection==false)//standard case
  {
    for(int i=0;i<3;i++)
    {
      df += (r2_xi(i)*r2_xi(i) - delta_r(i)*r2_xixi(i)) / norm_delta_r;
    }
  }
  else
  {
    for(int i=0;i<3;i++)
    {
      df += r1_xi(i)*r2_xi(i) / norm_delta_r;
    }
  }

  //check, if df=0: This can happen e.g. when the master beam 2 describes a circle geometry and the projectiong slave point
  //coincides with the cetern of the circle

  if(fabs(FADUTILS::CastToDouble(df))<COLINEARTOL)
    return false;
  else
    return true;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute normal vector in contact point                   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeNormal(LINALG::TMatrix<TYPE, 3, 1>& r1,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r2,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                    Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
                                                                    int contacttype)
{

  // compute non-unit normal
  LINALG::TMatrix<TYPE, 3, 1> delta_r = FADUTILS::DiffVector(r1,r2);

  // compute length of normal
  TYPE norm_delta_r = FADUTILS::VectorNorm<3>(delta_r);

  if (FADUTILS::CastToDouble(norm_delta_r) < NORMTOL)
    dserror("ERROR: Normal of length zero! --> change time step!");

  // unit normal
  LINALG::TMatrix<TYPE, 3, 1> normal(true);
  normal.Update(1.0/norm_delta_r,delta_r,0.0);

  TYPE gap=norm_delta_r - R1_ - R2_;

  //TODO
  if (FADUTILS::CastToDouble(gap)<-MAXPENETRATIONSAFETYFAC*(R1_+R2_) and numstep_>0)
  {
    std::cout << "element1_->Id(): " << element1_->Id() << std::endl;
    std::cout << "element2_->Id(): " << element2_->Id() << std::endl;
    std::cout << "gap: " << FADUTILS::CastToDouble(gap) << std::endl;
    std::cout << "xi: " << variables->GetCP().first << std::endl;
    std::cout << "eta: " << variables->GetCP().second << std::endl;
    //std::cout << "angle: " << FADUTILS::CastToDouble(BEAMCONTACT::CalcAngle(r1_xi,r2_xi) << std::endl;
    std::cout << "contacttype: " << contacttype << std::endl;
    dserror("Gap too small, danger of penetration. Choose smaller time step or higher penalty!");
  }

  variables->SetGap(gap);
  variables->SetNormal(normal);
  variables->SetAngle(BEAMCONTACT::CalcAngle(FADUTILS::CastToDouble<TYPE,3,1>(r1_xi),FADUTILS::CastToDouble<TYPE,3,1>(r2_xi)));

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute normal vector in contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if conact is active or inactive                    meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::CheckContactStatus(const double& gap)
{

  //First parameter for contact force regularization
  double g0 = bcparams_.get<double>("BEAMS_PENREGPARAM_G0",-1.0);
  bool contactflag = false;

  int penaltylaw = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW");

  if(penaltylaw==INPAR::BEAMCONTACT::pl_lp)
  {
    //linear penalty force law
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if(penaltylaw==INPAR::BEAMCONTACT::pl_qp)
  {
    //quadratic penalty force law
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if(penaltylaw==INPAR::BEAMCONTACT::pl_lpqp)
  {
    //penalty laws with regularization for positive gaps
    if(g0==-1.0)
      dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

    //Parameter to shift penalty law
    double gbar=bcparams_.get<double>("BEAMS_GAPSHIFTPARAM",0.0);
    TYPE g=gap+gbar;

    if (g < g0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if(penaltylaw==INPAR::BEAMCONTACT::pl_lpcp or penaltylaw==INPAR::BEAMCONTACT::pl_lpdqp or penaltylaw==INPAR::BEAMCONTACT::pl_lpep)
  {
    //penalty laws with regularization for positive gaps
    if(g0==-1.0)
      dserror("Invalid value of regularization parameter BEAMS_PENREGPARAM_G0!");

    if (gap < g0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }
  else if(penaltylaw==INPAR::BEAMCONTACT::pl_lnqp)
  {
    //penalty law with quadratic regularization for negative gaps
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }

  return contactflag;
}
/*----------------------------------------------------------------------*
 |  end: Check if conact is active or inactive
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if damping force is active or inactive             meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3contact<numnodes, numnodalvalues>::CheckDampingStatus(const double& gap)
{

  bool dampingcontactflag = false;

  if (DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_DAMPING")!=INPAR::BEAMCONTACT::bd_no)
  {
    //First parameter for contact force regularization
    double gd1 = bcparams_.get<double>("BEAMS_DAMPREGPARAM1",-1000.0);
    if(gd1==-1000.0)
      dserror("Damping parameter BEAMS_DAMPINGPARAM, BEAMS_DAMPREGPARAM1 and BEAMS_DAMPREGPARAM2 have to be chosen!");

    if (gap < gd1)
    {
      dampingcontactflag = true;
    }
    else
    {
      dampingcontactflag = false;
    }
  }

  return dampingcontactflag;
}
/*----------------------------------------------------------------------*
 |  end: Check if damping force is active or inactive
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get global dofs of a node                                 meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
std::vector<int> CONTACT::Beam3contact<numnodes, numnodalvalues>::GetGlobalDofs(const DRT::Node* node)
{
  // get dofs in beam contact discretization
  std::vector<int> cdofs = ContactDiscret().Dof(node);

  // get dofs in problem discretization via offset
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k=0;k<(int)(cdofs.size());++k)
  {
    pdofs[k]=(dofoffsetmap_.find(cdofs[k]))->second;
  }

  return pdofs;
}
/*----------------------------------------------------------------------*
 |  end: Get global dofs of a node
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Reset Uzawa-based Lagrange mutliplier                  meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::Resetlmuzawa()
{

  dserror("Uzawa algorithm does not work for beam3contact elements!");

  return;
}
/*----------------------------------------------------------------------*
 |  end: Reset Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Set all class variables                                meier 08/2014|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::SetClassVariables(Teuchos::ParameterList& timeintparams)
{
  iter_ = timeintparams.get<int>("iter",-10);
  numstep_ = timeintparams.get<int>("numstep",-10);

  if(iter_==-10.0 or numstep_==-10.0)
    dserror("Invalid time integration parameter!");

  cpvariables_.clear();
  gpvariables_.clear();
  epvariables_.clear();

  double kappa_max = timeintparams.get<double>("kappa_max",-1.0);
  if(kappa_max<0)
    dserror("Maximal curvature should be a positive value!");

  //TODO
//  //Check, if maximal curvature bound is exceeded:
//  double crosssection_to_curvature_ratio = max(R1_,R2_)*kappa_max;
//  if(crosssection_to_curvature_ratio>MAXCROSSSECTIONTOCURVATURE)
//    dserror("Curvature too large. Choose larger value MAXCROSSSECTIONTOCURVATURE and adapt shifting angles!");
}
/*----------------------------------------------------------------------*
 |  end: Set all class variables
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update Uzawa-based Lagrange mutliplier                 meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::Updatelmuzawa(const double& currentpp)
{
  dserror("Uzawa algorithm does not work for beam3contact so far!");

  return;
}
/*----------------------------------------------------------------------*
 |  end: Update Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)                        meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                                                                      Epetra_SerialDenseMatrix& newele2pos)
{
  for (int i=0;i<3*numnodalvalues;i++)
  {
    for (int j=0;j<numnodes;j++)
    {
      ele1pos_(3*numnodalvalues*j+i)=newele1pos(i,j);
      ele2pos_(3*numnodalvalues*j+i)=newele2pos(i,j);
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Update nodal coordinates (public)
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update nodal tangents for tangent smoothing (public)      meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::UpdateEleSmoothTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{

  //Tangent smoothing is only possible for Reissner beam elements --> dserror() otherwise
  if (numnodalvalues>1)
    dserror("Tangent smoothing only possible for Reissner beam elements (numnodalvalues=1)!!!");

  LINALG::Matrix<3*numnodes,1> elepos_aux(true);
  //Tangent smoothing only possible with data type double (not with Sacado FAD)
  for (int i=0;i<3*numnodes;i++)
    elepos_aux(i)=FADUTILS::CastToDouble(ele1pos_(i));

  nodaltangentssmooth1_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element1_,neighbors1_);

  elepos_aux.Clear();
  //Tangent smoothing only possible with data type double (not with Sacado FAD)
  for (int i=0;i<3*numnodes;i++)
    elepos_aux(i)=FADUTILS::CastToDouble(ele2pos_(i));

  nodaltangentssmooth2_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element2_,neighbors2_);

}
/*----------------------------------------------------------------------*
 |  end: Update nodal coordinates (public)
 *----------------------------------------------------------------------*/

template<const int numnodes , const int numnodalvalues>
double CONTACT::Beam3contact<numnodes, numnodalvalues>::GetJacobi(DRT::Element* element1)
{
  double jacobi = 1.0;
  const DRT::ElementType & eot1 = element1->ElementType();

  //The jacobi factor is only needed in order to scale the CPP condition. Therefore, we only use
  //the jacobi_ factor corresponding to the first gauss point of the beam element
  if (eot1 == DRT::ELEMENTS::Beam3ebType::Instance())
  {
    jacobi = (static_cast<DRT::ELEMENTS::Beam3eb*>(element1))->GetJacobi();
  }
  else if (eot1 == DRT::ELEMENTS::Beam3Type::Instance())
  {
    jacobi = (static_cast<DRT::ELEMENTS::Beam3*>(element1))->GetJacobi();
  }
  else if (eot1 == DRT::ELEMENTS::Beam3rType::Instance())
  {
    jacobi = (static_cast<DRT::ELEMENTS::Beam3r*>(element1))->GetJacobi();
  }
  else
  {
    std::cout << "      Warning: No valid jacobi weight in CPP supported by applied beam element!!!" << std::endl;
  }

  return jacobi;
}

#ifdef FADCHECKS
  /*----------------------------------------------------------------------*
   |  FAD-Check for Linearizations of contact point            meier 02/14|
   *----------------------------------------------------------------------*/
  template<const int numnodes , const int numnodalvalues>
  void CONTACT::Beam3contact<numnodes, numnodalvalues>::FADCheckLinXiAndLinEta(const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                                  const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                                  const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                                  const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                                                                  const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                                                                  const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                                  const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                                  const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                                  const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
  {
    LINALG::TMatrix<TYPE,2,1>f(true);
    LINALG::TMatrix<TYPE,3,1>t1_dummy(true);
    LINALG::TMatrix<TYPE,3,1>t2_dummy(true);

    // compute norm of difference vector to scale the equations
    // (this yields better conditioning)
    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
    // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
    double norm_delta_r = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(delta_r));

    EvaluateOrthogonalityCondition(f,delta_r,norm_delta_r,r1_xi,r2_xi,t1_dummy,t2_dummy);

    //**********************************************************************
    // we have to solve the following system of equations:
    //  _              _       _      _       _              _      _       _
    // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
    // |                | *  |         | =  |                 | * |        |
    // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
    //
    // this can be done easily because it is a linear 2x2-system.
    // we obtain the solution by inverting matrix L:
    //
    // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
    //
    //**********************************************************************

    const int dim1 = 3*numnodes*numnodalvalues;
    const int dim2 = 3*numnodes*numnodalvalues;

    // matrices to compute Lin_Xi and Lin_Eta
    LINALG::TMatrix<TYPE,2,2> L(true);
    LINALG::TMatrix<TYPE,2,2> L_inv(true);
    LINALG::TMatrix<TYPE,2,dim1+dim2> B(true);
    LINALG::TMatrix<TYPE,2,dim1+dim2> D(true);

    // compute L elementwise
    L(0,0)= f(0).dx(2*3*numnodes*numnodalvalues);
    L(0,1)= f(0).dx(2*3*numnodes*numnodalvalues+1);
    L(1,0)= f(1).dx(2*3*numnodes*numnodalvalues);
    L(1,1)= f(1).dx(2*3*numnodes*numnodalvalues+1);

    // invert L by hand
    TYPE det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
    if (FADUTILS::CastToDouble(FADUTILS::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
    L_inv(0,0) =  L(1,1) / det_L;
    L_inv(0,1) = -L(0,1) / det_L;
    L_inv(1,0) = -L(1,0) / det_L;
    L_inv(1,1) =  L(0,0) / det_L;

    for (int j=0;j<dim1+dim2;j++)
    {
      B(0,j)= -f(0).dx(j);
      B(1,j)= -f(1).dx(j);
    }

    // compute D = L^-1 * B
    D.Multiply(L_inv, B);

    std::cout << "linxi and lineta: " << std::endl;

    std::cout << D << std::endl;

    return;
  }
  /*----------------------------------------------------------------------*
   |  End: FAD-Check for Linearizations of contact point
   *----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*
   |  FAD-Check for Linearizations of CCP                      meier 02/14|
   *----------------------------------------------------------------------*/
  template<const int numnodes , const int numnodalvalues>
  void CONTACT::Beam3contact<numnodes, numnodalvalues>::FADCheckLinOrthogonalityCondition(const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                                                                          const double& norm_delta_r,
                                                                                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                                          const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                                          const LINALG::TMatrix<TYPE, 3, 1>& t1,
                                                                                          const LINALG::TMatrix<TYPE, 3, 1>& t2)
  {
    LINALG::TMatrix<TYPE,2,1>f(true);

    EvaluateOrthogonalityCondition(f, delta_r,norm_delta_r,r1_xi,r2_xi,t1,t2);

    LINALG::TMatrix<TYPE,2,2>df(true);

    for (int i=0;i<2;i++)
    {
      for (int j=0;j<2;j++)
      {
        df(i,j)=f(i).dx(2*3*numnodes*numnodalvalues+j);
      }
    }

    std::cout << "df_FAD: " << std::endl;

    std::cout << df << std::endl;

    return;
  }
  /*----------------------------------------------------------------------*
   |  End: FAD-Check for Linearizations of CPP
   *----------------------------------------------------------------------*/
#endif //#ifdef FADCHECKS

/*----------------------------------------------------------------------*
|  FD-Check of stiffness matrix                              meier 11/14|
*-----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::FDCheck(LINALG::SparseMatrix& stiffmatrix,
                                                              Epetra_Vector& fint,
                                                              const double& pp,
                                                              std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap,
                                                              Teuchos::ParameterList& timeintparams,
                                                              bool fdcheck)
{
  //This FD-Check is very general, since it applies the complete method "Evaluate" recursively.
  //Therefore, all changes within this class are automatically considered and have not to be
  //adapted in this finite difference check!
  if(fint.GlobalLength()>2*3*numnodes*numnodalvalues)
    dserror("So far, this FDCheck only works for simualtions with two elements!!!");

  Epetra_Vector fint1(fint);
  fint1.PutScalar(0.0);
  Epetra_Vector fint2(fint);
  fint2.PutScalar(0.0);

  LINALG::SparseMatrix stiffmatrix_analyt(stiffmatrix);
  stiffmatrix_analyt.PutScalar(0.0);

  LINALG::SparseMatrix stiffmatrix_dummy(stiffmatrix);
  stiffmatrix_dummy.PutScalar(0.0);

  LINALG::Matrix<2*3*numnodes*numnodalvalues,2*3*numnodes*numnodalvalues> stiffmatrix_fd(true);

  LINALG::Matrix<2*3*numnodes*numnodalvalues,2*3*numnodes*numnodalvalues> fint2_mat(true);

  double delta = 1.0e-10;

  std::cout << "undisturbed configuration: " << std::endl;

  this->Evaluate(stiffmatrix_analyt,fint1,pp,contactpairmap,timeintparams,true);

//  std::cout << std::setprecision(25) << "fint1: " << std::endl;

//  fint1.Print(std::cout);

  std::vector<double> xi1((int)gpvariables_.size(),0.0);
  std::vector<double> eta1((int)gpvariables_.size(),0.0);

  std::vector<LINALG::Matrix<2*3*numnodes*numnodalvalues,1> > xi2((int)gpvariables_.size());
  std::vector<LINALG::Matrix<2*3*numnodes*numnodalvalues,1> > eta2((int)gpvariables_.size());


  for (int i=0;i<(int)gpvariables_.size();i++)
  {
    xi1[i]=FADUTILS::CastToDouble(gpvariables_[i]->GetCP().first);
    eta1[i]=FADUTILS::CastToDouble(gpvariables_[i]->GetCP().second);
  }

  for(int dof=0;dof<2*3*numnodes*numnodalvalues;dof++)
  {
    std::cout << "disturbed configuration: " << std::endl;

    //Add delta
    if(dof<3*numnodes*numnodalvalues)
      ele1pos_(dof)+=delta;
    else
      ele2pos_(dof-3*numnodes*numnodalvalues)+=delta;

    fint2.PutScalar(0.0);
    stiffmatrix_dummy.PutScalar(0.0);

    this->Evaluate(stiffmatrix_dummy,fint2,pp,contactpairmap,timeintparams,true);

//    std::cout << std::setprecision(25) << "fint2: " << std::endl;
//
//    fint2.Print(std::cout);

    for(int i=0;i<(int)gpvariables_.size();i++)
    {
      xi2[i](dof)=FADUTILS::CastToDouble(gpvariables_[i]->GetCP().first);
      eta2[i](dof)=FADUTILS::CastToDouble(gpvariables_[i]->GetCP().second);
    }

    for(int i=0;i<2*3*numnodalvalues*numnodes;i++)
    {
      fint2_mat(i,dof)=fint2[i];

      if(fabs(fint2[i])<1.0e-10 and fabs(fint1[i])<1.0e-10)
        stiffmatrix_fd(i,dof) = 999999999999;
      else
        stiffmatrix_fd(i,dof) = -(fint2[i]-fint1[i])/delta;
    }

    //restore original displacements
    if(dof<3*numnodes*numnodalvalues)
      ele1pos_(dof)-=delta;
    else
      ele2pos_(dof-3*numnodes*numnodalvalues)-=delta;
  }

  std::cout << "FD_LIN: " << std::endl;

  for(int i=0;i<2*3*numnodes*numnodalvalues;i++)
  {
    for(int j=0;j<2*3*numnodes*numnodalvalues;j++)
    {
      std::cout << "row: " << i << "   " << "col: " << j << "   " << stiffmatrix_fd(i,j) << "   fint2: " << fint2_mat(i,j) << "   fint1: " << fint1[i] << std::endl;
    }
  }

  std::cout << "ANALYT_LIN: " << std::endl;

  std::cout << (*(stiffmatrix_analyt.EpetraMatrix())) << std::endl;
}
/*----------------------------------------------------------------------*
|  end: FD-Check of stiffness matrix
*-----------------------------------------------------------------------*/

//Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3contact<2,1>;
template class CONTACT::Beam3contact<3,1>;
template class CONTACT::Beam3contact<4,1>;
template class CONTACT::Beam3contact<5,1>;
template class CONTACT::Beam3contact<2,2>;
