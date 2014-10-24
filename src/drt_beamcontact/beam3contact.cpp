/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact.cpp
\brief One beam contact pair (two beam elements) consisting of several contact segments

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

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
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_inpar/inpar_statmech.H"

#include "Teuchos_TimeMonitor.hpp"

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
R1_(BEAMCONTACT::CalcEleRadius(element1)),
R2_(BEAMCONTACT::CalcEleRadius(element2)),
maxactivegap_(GetMaxActiveDist()),
maxsegdist1_(0.0),
maxsegdist2_(0.0),
numseg1_(1),
numseg2_(1)
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
    if(eot1 != DRT::ELEMENTS::Beam3Type::Instance() and eot1 != DRT::ELEMENTS::Beam3iiType::Instance())
      dserror("Tangent smoothing only implemented for beams of type beam3 and beam3ii!");

    //For both elements the 2 direct neighbor elements are determined and saved in the B3CNeighbor-Class
    //variables neighbors1_ and neighbors2_.
    {
      neighbors1_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element1);
      neighbors2_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element2);
    }
  }

  //Calculate initial length of beam elements (approximation for initially curved elements!)
  LINALG::TMatrix<TYPE,3,1> lvec1(true);
  LINALG::TMatrix<TYPE,3,1> lvec2(true);
  for(int i=0;i<3;i++)
  {
    lvec1(i)=(element1_->Nodes())[0]->X()[i]-(element1_->Nodes())[1]->X()[i];
    lvec2(i)=(element2_->Nodes())[0]->X()[i]-(element2_->Nodes())[1]->X()[i];
  }

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

  if(DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Strategy>(beamcontactparams,"BEAMS_STRATEGY") == INPAR::BEAMCONTACT::bstr_uzawa)
    dserror("Uzawa is not implemented for beam3contact elements so far!");

  if (DRT::INPUT::IntegralValue<int>(bcparams_,"BEAMS_DAMPING")!=INPAR::BEAMCONTACT::bd_no)
    dserror("Damping is not implemented for beam3contact elements so far!");

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
                                                                Teuchos::ParameterList& timeintparams)
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

  //All updates that have to be done in every iteration have do be done here,
  //since most of the elements leave directly after the closest point projection!
  SetClassVariables(timeintparams);

  //Subdevide the two elements in segments with linear approximation
  std::vector<LINALG::TMatrix<double,3,1> > endpoints1(0);
  std::vector<LINALG::TMatrix<double,3,1> > endpoints2(0);
  maxsegdist1_=CreateSegments(element1_, endpoints1, numseg1_);
  maxsegdist2_=CreateSegments(element2_, endpoints2, numseg2_);

//  std::cout << "endpoints1.size(): " << endpoints1.size() << std::endl;
//  std::cout << "endpoints2.size(): " << endpoints2.size() << std::endl;

  //Make pairs of close segments: Most of the pairs are already sorted out
  //at this point and don't have to be considered further in the following CPP
  //Additionally, we store the relative orientation of the pairs

  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> > closelargeanglesegments;
  std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> > closesmallanglesegments;
  GetCloseSegments(endpoints1,endpoints2,closesmallanglesegments,closelargeanglesegments,maxactivegap_);

  //std::cout << "closelargeanglesegments.size(): " << closelargeanglesegments.size() << std::endl;

  //**********************************************************************
  // (1) Closest Point Projection for all close large angle segments(CPP)
  //**********************************************************************
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

    //The method ClosestPointProjection() only delivers a valid solution (validpairfound=true), if eta1,eta2 \in [-1,1] and gap<maxactivegap_
    validpairfound=ClosestPointProjection(eta_left1,eta_left2,l1,l2,segmentdata,closestpoint);

    if(validpairfound)
    {
      cpvariables_.push_back(Teuchos::rcp (new CONTACT::Beam3contactvariables<numnodes,numnodalvalues>(closestpoint,leftpoint_ids,pp)));
      //std::cout << "eta_left1: " << eta_left1 << std::endl;
      //std::cout << "eta_left2: " << eta_left2 << std::endl;
      //std::cout << "closestpoint: " << closestpoint.first << " / " << closestpoint.second << std::endl;
    }
  }

  //std::cout << "cpvariables_.size(): " << cpvariables_.size() << std::endl;

  //Evaluate contact contribution for all closest points found before
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

    //std::cout << "eta1: " << eta1 << std::endl;
    //std::cout << "eta2: " << eta2 << std::endl;

    // update shape functions and their derivatives
    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);

    // call function to compute scaled normal and gap of contact point
    ComputeNormal(r1, r2, r1_xi, r2_xi, cpvariables_[numcp]);

    // call function to compute penalty force
    CalcPenaltyLaw(cpvariables_[numcp]);

    //std::cout << "cpvariables_[numcp]->GetNormal(): " << cpvariables_[numcp]->GetNormal() << std::endl;

//    std::cout << "numcp: " << numcp << std::endl;
//    std::cout << "xi: " << cpvariables_[numcp]->GetCP().first << std::endl;
//    std::cout << "eta: " << cpvariables_[numcp]->GetCP().second << std::endl;
//    std::cout << "gap: " << cpvariables_[numcp]->GetGap() << std::endl;
//    std::cout << "angle: " << BEAMCONTACT::CalcAngle(r1_xi,r2_xi)/M_PI*180 << std::endl;
//    std::cout << "r1_xi: " << r1_xi << std::endl;
//    std::cout << "r2_xi: " << r2_xi << std::endl;
//    std::cout << "|r1_xi|: " << r1_xi.Norm2() << std::endl;
//    std::cout << "|r2_xi|: " << r2_xi.Norm2() << std::endl;
//    std::cout << "r1_xi*r2_xi: " << BEAMCONTACT::ScalarProduct(r1_xi,r2_xi) << std::endl;
    //std::cout << "cpvariables_[numcp]->Getfp(): " << cpvariables_[numcp]->Getfp() << std::endl;

    // call function to compute contact contribution to residual vector
    EvaluateFcContact(&fint, N1, N2,cpvariables_[numcp]);

    //std::cout << "fint: " << std::endl;
    //fint.Print(std::cout);

    //std::cout << "stiffmatrix1: " << std::endl;
    //std::cout << *(stiffmatrix.EpetraMatrix()) << std::endl;

    // call function to compute contact contribution to stiffness matrix
    EvaluateStiffcContact(stiffmatrix, r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi,cpvariables_[numcp]);

    //std::cout << "stiffmatrix2: " << std::endl;
    //std::cout << *(stiffmatrix.EpetraMatrix()) << std::endl;


  }

  return (true);
}
/*----------------------------------------------------------------------*
 |  end: Evaluate the element
 *---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Calculate scalar contact force                           meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::CalcPenaltyLaw(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables)
{
    //First parameter for contact force regularization
    double g0 = bcparams_.get<double>("BEAMS_PENREGPARAM_G0",-1.0);
    TYPE fp=0.0;
    TYPE dfp=0.0;
    double pp=cpvariables->GetPP();
    TYPE gap=cpvariables->GetGap();

    switch (DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW"))
    {
      case INPAR::BEAMCONTACT::pl_lp:               //linear penalty force law
      {
        fp= - pp*gap;
        dfp=-pp;

        break;
      }
      case INPAR::BEAMCONTACT::pl_qp:               //quadratic penalty force law
      {
        fp=pp*gap*gap;
        dfp=2*pp*gap;

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

        double f0=g0*pp/2.0;
        double factor_a=pp/(g0) -f0/(g0*g0);
        double factor_b=-pp;
        double factor_c=f0;
        if(gap>0)
        {
          //std::cout << "Regularized Penalty!" << std::endl;
          fp=factor_a*gap*gap+factor_b*gap+factor_c;
          dfp=2*factor_a*gap+factor_b;
        }
        else
        {
          //std::cout << "Linear Penalty!" << std::endl;
          fp=f0 - pp*gap;
          dfp=-pp;
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

    cpvariables->Setfp(fp);
    cpvariables->Setdfp(dfp);

  return;
}
/*----------------------------------------------------------------------*
 |  end: Calculate scalar contact force
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Subdivide elements into segments for CPP                 meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
double CONTACT::Beam3contact<numnodes, numnodalvalues>::CreateSegments(DRT::Element* ele,
                                                                       std::vector<LINALG::TMatrix<double,3,1> >& endpoints_final,
                                                                       int& numsegment)
{
  //endpoints of the segments
  std::vector<LINALG::TMatrix<double,3,1> > endpoints((int)MAXNUMSEG+1,LINALG::TMatrix<double,3,1>(true));

  numsegment=1;
  double deltaxi=2.0;
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
      xi1= 0.0;
      xi2= 0.0;
      xi1= -1.0 +i/((double)numsegment)*2.0;
      xi2= -1.0 +(i+1)/((double)numsegment)*2.0;  //The cast to double is necessary here to avoid integer round-off
      r1=BEAMCONTACT::CastToDouble<TYPE,3,1>(r(xi1,ele));
      r2=BEAMCONTACT::CastToDouble<TYPE,3,1>(r(xi2,ele));
      t1=BEAMCONTACT::CastToDouble<TYPE,3,1>(r_xi(xi1,ele));
      t2=BEAMCONTACT::CastToDouble<TYPE,3,1>(r_xi(xi2,ele));
      rm=BEAMCONTACT::CastToDouble<TYPE,3,1>(r((xi1+xi2)/2.0,ele));
      endpoints[i]=r1;
      endpoints[i+1]=r2;
      l=BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::DiffVector(r1,r2));
      double safetyfac2=1.1;
      segdist=safetyfac2*l/2.0*(1.0-cos(SEGANGLE))/sin(SEGANGLE);

      if(segdist>maxsegdist)
        maxsegdist=segdist;
      if(!CheckSegment(r1,t1,r2,t2,rm,segdist))
        moresegments = true;
    }
    deltaxi=deltaxi/2;
    numsegment=numsegment*2;
    if (numsegment>(int)MAXNUMSEG)
      dserror("Not more segments than MAXNUMSEG per element possible! Increase MAXNUMSEG or apply finer discretization!");
  }
  numsegment=numsegment/2;

  //std::cout << "numsegment: " << numsegment << std::endl;

  endpoints_final.resize(numsegment+1);

  for (int i=0;i<numsegment+1;i++)
  {
    endpoints_final[i]=endpoints[i];
  }

  //Update class variable

  if(maxsegdist<R1_ or maxsegdist<R2_)
    std::cout << "      Warning: maxsegist is smaller than the beam radius, smaller number of segments is sufficient. Choose smaller segment angle!" << std::endl;

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
//
//  std::cout << "r1: " << r1 << std::endl;
//  std::cout << "r2: " << r2 << std::endl;
//  std::cout << "t1: " << t1 << std::endl;
//  std::cout << "t2: " << t2 << std::endl;
//  std::cout << "rm: " << rm << std::endl;

  //Calculate tangent and midpint of linear nodal interpolation
  for(int i=0;i<3;i++)
  {
    t_lin(i)=r2(i)-r1(i);
    rm_lin(i)=(r2(i)+r1(i))/2.0;
  }
//
//  std::cout << "t_lin: " << t_lin << std::endl;
//  std::cout << "rm_lin: " << rm_lin << std::endl;

  LINALG::TMatrix<double,3,1> diffvec(true);
  diffvec = BEAMCONTACT::DiffVector(rm_lin, rm);
  dist = (double)BEAMCONTACT::VectorNorm<3>(diffvec);
  angle1 = (double)BEAMCONTACT::CalcAngle(t1,t_lin);
  angle2 = (double)BEAMCONTACT::CalcAngle(t2,t_lin);

//  std::cout << "angle1: " << angle1 << std::endl;
//  std::cout << "angle2: " << angle2 << std::endl;

  if(fabs(angle1)<SEGANGLE and fabs(angle2)<SEGANGLE) //segment distribution is fine enough
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
                                                                        double maxactivedist)
{
  LINALG::TMatrix<double,3,1> t1(true);
  LINALG::TMatrix<double,3,1> t2(true);
  LINALG::TMatrix<double,3,1> r1_a(true);
  LINALG::TMatrix<double,3,1> r1_b(true);
  LINALG::TMatrix<double,3,1> r2_a(true);
  LINALG::TMatrix<double,3,1> r2_b(true);
  double angle(0.0);

  //Safety factor for determination of close segments
  double safetyfac = 1.1;
  //Distance at which intersection happens
  double distancelimit=safetyfac*(maxsegdist1_+maxsegdist2_+maxactivedist+R1_+R2_);

//  std::cout << "maxactivedist: " << maxactivedist << std::endl;
//  std::cout << "R1_: " << R1_ << std::endl;
//  std::cout << "R2_: " << R2_ << std::endl;
//  std::cout << "maxsegdist1_: " << maxsegdist1_ << std::endl;
//  std::cout << "maxsegdist2_: " << maxsegdist2_ << std::endl;
//  std::cout << "distancelimit: " << distancelimit << std::endl;
//
//  for (int i=0;i<2;i++)
//  {
//    std::cout << "endpoints1[i]" << endpoints1[i] << std::endl;
//    std::cout << "endpoints2[i]" << endpoints2[i] << std::endl;
//  }

  //TODO: This check is implemented in a brute force way. However, this should be efficient enough as long
  //as the number of segments per element remains small!
  for (int i=0;i<(int)endpoints1.size()-1;i++)
  {
    r1_a=endpoints1[i];
    r1_b=endpoints1[i+1];
    t1=BEAMCONTACT::DiffVector(r1_b,r1_a);
    for (int j=0;j<(int)endpoints2.size()-1;j++)
    {
      r2_a=endpoints2[j];
      r2_b=endpoints2[j+1];
      t2=BEAMCONTACT::DiffVector(r2_b,r2_a);
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

          //Check, if we have sensible parameters
          if(DELTASMALLANGLE < ANGLETOL or DELTALARGEANGLE < ANGLETOL)
            dserror("Tolerance ANGLETOL has to be chosen smaller!");

          //Add new small angle pair since we always have: DELTASMALLANGLE > ANGLETOL and DELTALARGEANGLE > ANGLETOL
          closesmallanglesegments[std::make_pair(i,j)]=segmentdata;
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
          if(fabs(angle)<DELTASMALLANGLE)
            closesmallanglesegments[std::make_pair(i,j)]=segmentdata;
          if(fabs(angle)>DELTALARGEANGLE)
            closelargeanglesegments[std::make_pair(i,j)]=segmentdata;
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
                                                                              std::pair<TYPE,TYPE>& solutionpoints)
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
      delta_r=BEAMCONTACT::DiffVector(r1,r2);

      // compute norm of difference vector to scale the equations
      // (this yields better conditioning)
      // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
      // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
      double norm_delta_r = BEAMCONTACT::CastToDouble(BEAMCONTACT::VectorNorm<3>(delta_r));
      gap=norm_delta_r-R1_-R2_;

      // The closer the beams get, the smaller is norm_delta_r, but
      // norm_delta_r is not allowed to be too small, else numerical problems occur.
      // It can happen quite often that the centerlines of two beam elements of the same physical beam
      // cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they will be sorted out later anyways.
      //std::cout << "norm_delta_r: " << norm_delta_r << std::endl;
      if (norm_delta_r < NORMTOL)
      {
        // this exludes pairs with IDs i and i+2, i.e. contact with the next but one element
        if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(eta1)) + BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(eta2)) < NEIGHBORTOL)
        {
          dserror("Beam axis identical, choose smaller time step!");
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
      //residual = f(0)*f(0)/(jacobi1*jacobi1) + f(1)*f(1)/(jacobi2*jacobi2);
      residual = sqrt(BEAMCONTACT::CastToDouble(f(0)*f(0)/(jacobi1*jacobi1) + f(1)*f(1)/(jacobi2*jacobi2)));

//      std::cout << "iter: " << iter << std::endl;
//      std::cout << "residual: " << residual << std::endl;

      if(iter==1)
        residual0=residual;

      // check if Newton iteration has converged
      if (BEAMCONTACT::CastToDouble(residual) < BEAMCONTACTTOL)
      { converged=true;
        break;
      }

      // evaluate Jacobian of f at current eta1, eta2
      // Note: Parallel elements can not be handled with this beam contact formulation;
      EvaluateLinOrthogonalityCondition(df, dfinv, delta_r, norm_delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, t1, t2, t1_xi, t2_xi, elementscolinear);

      #ifdef FADCHECKS
        std::cout << "f: " << f << std::endl;
        std::cout << "df: " << df << std::endl;
        BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
        FADCheckLinOrthogonalityCondition(delta_r, norm_delta_r, r1_xi,r2_xi, t1, t2);
      #endif

      if (elementscolinear) break;

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

//    std::cout << "numstartpoint: " << numstartpoint << std::endl;
//    std::cout << "iter: " << iter << std::endl;
//    std::cout << "eta1: " << eta1 << std::endl;
//    std::cout << "eta2: " << eta2 << std::endl;

    // Newton iteration unconverged after BEAMCONTACTMAXITER
    if (!converged)
    {
      //Check for cases where the residual has already been decreased considerably, but the tolerance has not been reached:
      if(residual/residual0 < 1.0e-08 and fabs(eta1)< 1.0+XIETATOL and fabs(eta2)< 1.0+XIETATOL)
      {
        std::cout << "iter: " << iter << std::endl;
        std::cout << "residual0: " << residual0 << std::endl;
        std::cout << "lastresidual: " << lastresidual << std::endl;
        std::cout << "residual: " << residual << std::endl;
        std::cout << "eta1: " << eta1 << std::endl;
        std::cout << "eta2: " << eta2 << std::endl;
        dserror("Relative CPP residual norm is smaller than 1.0e-08 but Newton is not converged. Adapt BEAMCONTACTTOL or the maximal number BEAMCONTACTMAXITER of iterations!");
      }

      eta1 = 1e+12;
      eta2 = 1e+12;
    }
    else
    {
      //if we have already found a converged solution with valid closest points eta1,eta2 \in [-1.0;1.0], we can finish here and don't have to apply more starting points
      if( eta_left1<=eta1 and eta1<= eta_right1 and eta_left2<=eta2 and eta2<=eta_right2 and (CheckContactStatus(gap) or CheckDampingStatus(gap)))
      {
        solutionpoints.first=BEAMCONTACT::CastToDouble(eta1);
        solutionpoints.second=BEAMCONTACT::CastToDouble(eta2);
        validpairfound=true;
        break;
      }
    }
  }//for (int numstartpoint=0;numstartpoint<startingpoints.size();numstartpoint++)

  // Set eta1 and eta2 as primary variables for automatic differentiation
  // The dependence between the infinitesimal changes delta eta1 and delta eta2 and the
  // the increments of the primary displacement variables delta disp have to be given explicitly, since
  // no explicit relation between the finite quantities eta1, eta2 and disp exists.
  // The latter would have been necessary if the full linearization had to be computed directly with Sacado!!!
  #ifdef AUTOMATICDIFF
    BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(solutionpoints.first, solutionpoints.second);
  #endif

  return validpairfound;

}
/*----------------------------------------------------------------------*
|  end: Closest point projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   meier 10/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::EvaluateFcContact(Epetra_Vector* fint,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables,
                                                                        LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD,
                                                                        LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD)
{
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

  TYPE gap = cpvariables->GetGap();
  LINALG::TMatrix<TYPE,3,1> normal = cpvariables->GetNormal();
  TYPE fp = cpvariables->Getfp();

  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (CheckContactStatus(BEAMCONTACT::CastToDouble(gap)))
  {
    DoNotAssemble=false;
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    for (int i=0;i<dim1;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc1(i) +=  N1(j,i)*normal(j)*fp;
      }
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    for (int i=0;i<dim2;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc2(i) +=  -N2(j,i)*normal(j)*fp;
      }
    }
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
        fcontact1[i] = BEAMCONTACT::CastToDouble(fc1(i));
    }
    for (int i=0;i<dim2;++i)
    {
        fcontact2[i] = BEAMCONTACT::CastToDouble(fc2(i));
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
                                                                             Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables)
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
  TYPE gap = cpvariables->GetGap();

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
  if (CheckContactStatus(BEAMCONTACT::CastToDouble(gap)))
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

    LINALG::TMatrix<TYPE, 3, 1> delta_r = BEAMCONTACT::DiffVector(r1,r2);
    TYPE norm_delta_r = BEAMCONTACT::VectorNorm<3>(delta_r);
    LINALG::TMatrix<TYPE, 3, 1> normal = cpvariables->GetNormal();
    TYPE fp = cpvariables->Getfp();
    TYPE dfp = cpvariables->Getdfp();

    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************
    // linearization of contact point
    ComputeLinXiAndLinEta(delta_xi,delta_eta,delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);

    // linearization of gap function which is equal to delta d
    ComputeLinGap(delta_gap,delta_xi,delta_eta,delta_r,norm_delta_r,r1_xi,r2_xi,N1,N2);

    // linearization of normal vector
    ComputeLinNormal(delta_n,delta_xi,delta_eta,delta_r,r1_xi,r2_xi,N1,N2);

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
        stiffc1(i,j) += basicstiffweightfac* dfp * N1T_normal(i) * delta_gap(j);
      }
    }

    //The geoemtric part is only applied for gap_ < -BEAMS_BASICSTIFFGAP
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
              stiffc1(j,k) += fp*N1(i,j)*delta_n(i,k);
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
          stiffc1(i,j) += fp*N1xiT_normal(i)*delta_xi(j);
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
        stiffc2(i,j) += -basicstiffweightfac* dfp * N2T_normal(i) * delta_gap(j);
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
              stiffc2(j,k) += -fp*N2(i,j)*delta_n(i,k);
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
          stiffc2(i,j) += -fp*N2xiT_normal(i)*delta_eta(j);
        }
      }
    }
    //*************End of standard linearization of penalty contact forces****************

    // automatic differentiation for debugging
    #ifdef AUTOMATICDIFF
      LINALG::TMatrix<TYPE, dim1, 1> fc1_FAD(true);
      LINALG::TMatrix<TYPE, dim2, 1> fc2_FAD(true);
      EvaluateFcContact(NULL, N1, N2, cpvariables, &fc1_FAD, &fc2_FAD);
      for (int j=0;j<dim1+dim2;j++)
      {
        for (int i=0;i<dim1;i++)
          stiffc1_FAD(i,j) = (fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc1_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
        for (int i=0;i<dim2;i++)
          stiffc2_FAD(i,j) = (fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc2_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
      }

      std::cout << "Pair: " << element1_->Id() << " / " << element2_->Id() << std::endl;

      std::cout << "stiffc1: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          std::cout << stiffc1(i,j).val() << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << "stiffc1_FAD: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          std::cout << stiffc1_FAD(i,j).val() << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << "stiffc2: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          std::cout << stiffc2(i,j).val() << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << "stiffc2_FAD: " << std::endl;
      for (int i=0;i<dim1;i++)
      {
        for (int j=0;j<dim1+dim2;j++)
        {
          std::cout << stiffc2_FAD(i,j).val() << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
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
          stiffcontact1(i,j) = -BEAMCONTACT::CastToDouble(stiffc1(i,j));
        for (int i=0;i<dim2;i++)
          stiffcontact2(i,j) = -BEAMCONTACT::CastToDouble(stiffc2(i,j));
      }
    #else
      for (int j=0;j<dim1+dim2;j++)
      {
        for (int i=0;i<dim1;i++)
          stiffcontact1(i,j) = -BEAMCONTACT::CastToDouble(stiffc1_FAD(i,j));
        for (int i=0;i<dim2;i++)
          stiffcontact2(i,j) = -BEAMCONTACT::CastToDouble(stiffc2_FAD(i,j));
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
  L(0,0)=::BEAMCONTACT::ScalarProduct(r1_xi, r1_xi) + ::BEAMCONTACT::ScalarProduct(delta_r, r1_xixi);
  L(1,1)=-::BEAMCONTACT::ScalarProduct(r2_xi, r2_xi) + ::BEAMCONTACT::ScalarProduct(delta_r, r2_xixi);
  L(0,1)=-::BEAMCONTACT::ScalarProduct(r2_xi, r1_xi);
  L(1,0)=-L(0,1);

  // invert L by hand
  TYPE det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
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

  TYPE norm_delta_r = BEAMCONTACT::VectorNorm<3>(delta_r);
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
 | compute coordinate at given curve point                  meier 10/14|
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
 | end: compute coordinate at given curve point              meier 02/14|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute coordinate at given curve point                  meier 10/14|
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
 | end: compute coordinate at given curve point              meier 02/14|
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

//  // store coordinates of contact point into class variables
//  for(int i=0;i<3;i++)
//  {
//    r1_(i)=r1(i);
//    r2_(i)=r2(i);
//    r1_xi_(i)=r1_xi(i);
//    r2_xi_(i)=r2_xi(i);
//  }

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
  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_df))<COLINEARTOL)
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
 |  Compute normal vector in contact point                   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3contact<numnodes, numnodalvalues>::ComputeNormal(LINALG::TMatrix<TYPE, 3, 1>& r1,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r2,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                                                    LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                                                    Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables)
{

  // compute non-unit normal
  LINALG::TMatrix<TYPE, 3, 1> delta_r = BEAMCONTACT::DiffVector(r1,r2);

  // compute length of normal
  TYPE norm_delta_r = BEAMCONTACT::VectorNorm<3>(delta_r);

  if (BEAMCONTACT::CastToDouble(norm_delta_r) < NORMTOL)
    dserror("ERROR: Normal of length zero! --> change time step!");

  // unit normal
  LINALG::TMatrix<TYPE, 3, 1> normal(true);
  normal.Update(1.0/norm_delta_r,delta_r,0.0);

  TYPE gap=norm_delta_r - R1_ - R2_;

  if (BEAMCONTACT::CastToDouble(gap)<-0.9*(R1_+R2_))
    dserror("Gap to small, danger of penetration. Choose smaller time step!");

  cpvariables->SetGap(gap);
  cpvariables->SetNormal(normal);
  cpvariables->SetAngle(BEAMCONTACT::CalcAngle(r1_xi,r2_xi));

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

  if(penaltylaw==INPAR::BEAMCONTACT::pl_qp)
  {
    //quadratic penalty force law
    if (gap < 0)
    {
      contactflag = true;
    }
    else
      contactflag = false;
  }

  if(penaltylaw==INPAR::BEAMCONTACT::pl_lpqp or penaltylaw==INPAR::BEAMCONTACT::pl_lpcp or penaltylaw==INPAR::BEAMCONTACT::pl_lpdqp or penaltylaw==INPAR::BEAMCONTACT::pl_lpep)
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

  if(penaltylaw==INPAR::BEAMCONTACT::pl_lnqp)
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

  if(iter_==-10.0)
    dserror("Invalid time integration parameter!");

  cpvariables_.clear();
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
    elepos_aux(i)=BEAMCONTACT::CastToDouble(ele1pos_(i));

  nodaltangentssmooth1_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element1_,neighbors1_);

  elepos_aux.Clear();
  //Tangent smoothing only possible with data type double (not with Sacado FAD)
  for (int i=0;i<3*numnodes;i++)
    elepos_aux(i)=BEAMCONTACT::CastToDouble(ele2pos_(i));

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
  else if (eot1 == DRT::ELEMENTS::Beam3iiType::Instance())
  {
    jacobi = (static_cast<DRT::ELEMENTS::Beam3ii*>(element1))->GetJacobi();
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
    double norm_delta_r = BEAMCONTACT::CastToDouble(BEAMCONTACT::VectorNorm<3>(delta_r));

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
    if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
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

//Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3contact<2,1>;
template class CONTACT::Beam3contact<3,1>;
template class CONTACT::Beam3contact<4,1>;
template class CONTACT::Beam3contact<5,1>;
template class CONTACT::Beam3contact<2,2>;
