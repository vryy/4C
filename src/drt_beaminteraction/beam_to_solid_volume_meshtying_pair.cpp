/*----------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_pair.cpp

\brief meshtying element for meshtying between a 3D beam and a 3D solid element

\level 3

\maintainer Alexander Popp
*/
/*----------------------------------------------------------------------------*/


#include "beam_to_solid_volume_meshtying_pair.H"

#include "beam_contact_pair.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "../linalg/linalg_utils.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"

#include "beam_contact_evaluation_data.H"
#include "beam_to_solid_volume_meshtying_evaluation_data.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::BeamToSolidVolumeMeshtyingPair():
BeamContactPair(),
atleastone_(false),
n_gp_(0),
ele1gid_(0),
contactflag_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Setup()
{

  CheckInit();

  // call setup of base class first
  BeamContactPair::Setup();

  ele1gid_ = Element1()->Id();

  // Initialize the vector that keeps track of which Gauss points have a valid projection in this pair
  // get number of Gaus points of the beam element
  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GAUSSRULE);
  n_gp_ = gaussPoints.nquad;

  // declare an entry in the Gauss Point Evaluate Tracker for every BTSVOLMT pair
  std::map<int, std::vector<bool> >* GaussPointEvaluateTracker = EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointEvaluateTracker();
  std::map<int, std::vector<bool> >* GaussPointProjectionTracker = EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointProjectionTracker();


  std::vector<bool> temp;
  temp.resize(n_gp_, false);

  (*GaussPointEvaluateTracker)[ele1gid_] = temp;
  (*GaussPointProjectionTracker)[ele1gid_] = temp;


  /*************************************************************************************************
   *************************************************************************************************/
  // Initialize reference nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < 3*numnodes*numnodalvalues; i++)
    ele1posref_(i) = 0.0;

  // Set reference nodal positions (and tangents) for beam element
  for (unsigned int n=0;n<numnodes;++n)
  {
    const DRT::Node* node = Element1()->Nodes()[n];
    for (int d=0;d<3;++d)
      ele1posref_(3*numnodalvalues*n + d) = node->X()[d];

    // tangents
    if (numnodalvalues==2)
    {
      LINALG::Matrix<3,1> tan;
      const DRT::ElementType& eot = Element1()->ElementType();
      if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        dserror("ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3 element");
      }
      else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(Element1());
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror("ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3r element w/o Hermite CL");

      }
      else if (eot == DRT::ELEMENTS::Beam3kType::Instance())
      {
        const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(Element1());
        tan = ele->Tref()[n];
      }
      else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
      {
        const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(Element1());
        tan = ele->Tref()[n];
      }
      else
      {
        dserror("ERROR: Beam3tosolidmeshtying: Invalid beam element type");
      }

      for (int d=0;d<3;++d)
        ele1posref_(3*numnodalvalues*n+d+3) = tan(d,0);
    }
  }

  // Initialize reference nodal positions for solid surface element
  for (unsigned int i = 0; i < 3*numnodessol; i++)
    ele2posref_(i) = 0.0;

  // Set reference nodal positions for solid surface element
  for (unsigned int n=0;n<numnodessol;++n)
  {
    const DRT::Node* node = Element2()->Nodes()[n];
    for (int d=0;d<3;++d)
      ele2posref_(3*n+d) = node->X()[d];
  }

  // Initialize current nodal positions (and tangents) for beam element
  for (unsigned int i = 0; i < 3*numnodes*numnodalvalues; i++)
    ele1pos_(i) = 0.0;

  // Initialize current nodal positions for solid surface element
  for (unsigned int i = 0; i < 3*numnodessol; i++)
    ele2pos_(i) = 0.0;


  // Init the const-values vectors --> They have size nquad
  // initialize
  beamjacobi_.resize(n_gp_);
  xi1_.resize(n_gp_);
  xi2_.resize(n_gp_);
  xi3_.resize(n_gp_);
  eta_.resize(n_gp_);
  wgp_.resize(n_gp_);
  localGaussPointTracker_.resize(n_gp_, false);

  // Set atleastone flag to false
  atleastone_ = false;
  /*************************************************************************************************
   *************************************************************************************************/


  issetup_ = true;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::PreEvaluate()
{

  bool output = false;

  /******************************************************************************************
   Initialize some data
   ******************************************************************************************/
  bool proj;

  TYPEBTS xi1;
  TYPEBTS xi2;
  TYPEBTS xi3;

  std::map<int, std::vector<bool> >* GaussPointProjectionTracker =
      EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointProjectionTracker();

  // get number of Gaus points of the beam element
  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GAUSSRULE);
  /******************************************************************************************
   End: Initialize some data
   ******************************************************************************************/

  // check if there is already an entry for that element pair in the Gauss point tracker map
  if (GaussPointProjectionTracker->count(ele1gid_) == 0)
  {
    std::vector<bool> temp;
    temp.resize(n_gp_, false);

    (*GaussPointProjectionTracker)[ele1gid_] = temp;
  }

  // iterate over all Gauss points of the beam element
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    // check if the Gauss point already has a valid projection
    if ((*GaussPointProjectionTracker)[ele1gid_][i_gp] == false)
    {
      proj = false;
      // fill beam jacobi and the gauss point and weight vector
      double x_gp = gaussPoints.qxg[i_gp][0];
      GetBeamJacobi(x_gp, beamjacobi_[i_gp]);
      eta_[i_gp] = x_gp;
      wgp_[i_gp] = gaussPoints.qwgt[i_gp];

      Projection(xi1, xi2, xi3, x_gp, proj);

      if (proj)
      {

        if (output)
        {
        std::cout << "Projection was valid!" << std::endl;

        std::cout << "Gauss Point " << i_gp << " with eta = " << x_gp << " was projected at xi1 = " << xi1 << " xi2 = " << xi2 << " xi3 = " << xi3 << std::endl;
        }

        // set class variables xi_i
        xi1_[i_gp] = xi1;
        xi2_[i_gp] = xi2;
        xi3_[i_gp] = xi3;

        // set this Gauss point in the Projection map to true
        (*GaussPointProjectionTracker)[ele1gid_][i_gp] = true;

        // set this bool in the local Projection map to true aswell
        localGaussPointTracker_[i_gp] = true;

        // at least one GP projects
        if (!atleastone_)
        {
          atleastone_ = true;
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Evaluate(
    LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{


  /******************************************************************************************
   Initialize some data
   ******************************************************************************************/
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // resize and initialize variables to zero
  if (forcevec1 != NULL)
    forcevec1->Size(dim1);
  if (forcevec2 != NULL)
    forcevec2->Size(dim2);

  if (stiffmat11 != NULL)
    stiffmat11->Shape(dim1,dim1);
  if (stiffmat12 != NULL)
    stiffmat12->Shape(dim1,dim2);
  if (stiffmat21 != NULL)
    stiffmat21->Shape(dim2,dim1);
  if (stiffmat22 != NULL)
    stiffmat22->Shape(dim2,dim2);


  // in case NULL-pointers are passed for stiffmatij, local correlates have to be defined
  // because the computation of forcevec depends on stiffmat...
  LINALG::SerialDenseMatrix local_stiffmat11;
  LINALG::SerialDenseMatrix local_stiffmat12;
  LINALG::SerialDenseMatrix local_stiffmat21;
  LINALG::SerialDenseMatrix local_stiffmat22;

  double pp = Params()->BeamToSolidVolumeMeshtyingParams()->BeamToSolidVolumeMeshtyingPenaltyParam();

  std::map<int, std::vector<bool> >* GaussPointProjectionTracker =
      EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointProjectionTracker();

  std::map<int, std::vector<bool> >* GaussPointEvaluateTracker =
      EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointEvaluateTracker();

  TYPEBTS xi1;
  TYPEBTS xi2;
  TYPEBTS xi3;
  TYPEBTS eta;
  TYPEBTS jac;

  // Vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2(true);                     // = N2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1(true);                 // = N2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2(true);                 // = N2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi3(true);                 // = N2,xi3

  // Coords and derivatives of beam and surface element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 1> x2(true);                                 // = x2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1(true);                             // = x2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2(true);                             // = x2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi3(true);                             // = x2,xi3
  /******************************************************************************************
   End: Initialize some data
   ******************************************************************************************/


  // count number of valid projections
  int numberOfValidProjectionsPair = 0;
  int numberOfValidProjectionsBeamElement = 0;
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    if (localGaussPointTracker_[i_gp] == 1)
    {
      numberOfValidProjectionsPair++;
    }
    if ((*GaussPointProjectionTracker)[ele1gid_][i_gp] == true)
    {
      numberOfValidProjectionsBeamElement++;
    }
  }

  // Check if a boundary segmentation is necessary
  if ( numberOfValidProjectionsBeamElement < n_gp_ && numberOfValidProjectionsPair > 0)
  {
    EvaluateBoundarySegmentation();
  }


  if (atleastone_)
  {
  // does at least one of the Gauss points of the beam elements project??
  // if none projects, we do nothing because this beam element does not contribute
  // to the beam to solid meshtying forces

    for (int i_gp = 0; i_gp < n_gp_; i_gp++)
    {

      if (localGaussPointTracker_[i_gp] == true
          && (*GaussPointEvaluateTracker)[ele1gid_][i_gp] == false)
      {

        // re-initialize local stiffness matrices
        local_stiffmat11.Shape(dim1,dim1);
        local_stiffmat12.Shape(dim1,dim2);
        local_stiffmat21.Shape(dim2,dim1);
        local_stiffmat22.Shape(dim2,dim2);

        // Get constant values from projection
        const double w_gp = wgp_[i_gp];
        eta = eta_[i_gp];
        jac = beamjacobi_[i_gp];
        xi1 = xi1_[i_gp];
        xi2 = xi2_[i_gp];
        xi3 = xi3_[i_gp];


        // Evaluate the part of the stiffness matrix meshtying vector
        GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
        GetSolShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi3, xi1, xi2, xi3);


        EvaluateStiffmtMeshtying(local_stiffmat11, local_stiffmat12, local_stiffmat21,
            local_stiffmat22, w_gp, N1, N2, jac, pp);

        if (stiffmat11 != NULL && stiffmat12 != NULL && stiffmat21 != NULL && stiffmat22 != NULL)
        {
          FillElementStiffmtMeshtying(local_stiffmat11, local_stiffmat12, local_stiffmat21,
              local_stiffmat22, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
        }


        if (forcevec1 != NULL && forcevec2 != NULL)
        {
          EvaluateFmtMeshtying(*forcevec1, *forcevec2, local_stiffmat11, local_stiffmat12,
              local_stiffmat21, local_stiffmat22);
        }


        // Update GaussPointEvaluateTracker
        (*GaussPointEvaluateTracker)[ele1gid_][i_gp] = true;



      }
    }
  }


  if ( ( forcevec1 != NULL and forcevec1->NormInf() > 0 ) or
      ( forcevec2 != NULL and forcevec2->NormInf() > 0 ) or
      ( stiffmat11 != NULL and stiffmat11->NormInf() > 0 ) or
      ( stiffmat12 != NULL and stiffmat12->NormInf() > 0 ) or
      ( stiffmat21 != NULL and stiffmat21->NormInf() > 0 ) or
      ( stiffmat22 != NULL and stiffmat22->NormInf() > 0 ) )
  {
    // TODO: Maybe think of a better place to set this bool?
    contactflag_ = true;
    return true;
  }
  else
  {
    return false;
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::EvaluateStiffmtMeshtying(
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22,
    const double& w_gp,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2,
    TYPEBTS& jacobi,
    const double& pp)
{

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;


  // Evaluate meshtying stiffness for beam element N_1^T * N_1
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1; j++)
      for (int k = 0; k < dim1; k++)
        stiffmat11(j,k) += pp * jacobi * w_gp * N1(i,j) * N1(i,k);

  // Evaluate meshtying stiffness for beam element "mixed" N_1^T * N_2
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1; j++)
      for (int k = 0; k < dim2; k++)
        stiffmat12(j,k) += pp * jacobi * w_gp * N1(i,j) * ( - N2(i,k) );

  // Evaluate meshtying stiffness for surface element "mixed" N_2^T * N_1
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1; j++)
      for (int k = 0; k < dim2; k++)
        stiffmat21(k,j) += pp * jacobi * w_gp * N2(i,k) * ( - N1(i,j) );

  // Evaluate meshtying stiffness for surface element N_2^T * N_2
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim2; j++)
      for (int k = 0; k < dim2; k++)
        stiffmat22(k,j) += pp * jacobi * w_gp * N2(i,k) * N2(i,j);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::EvaluateFmtMeshtying(
    LINALG::SerialDenseVector& forcevec1,
    LINALG::SerialDenseVector& forcevec2,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;


  // Evaluate meshtying forces for beam element
  for (int i = 0; i < dim1; i++)
  {
    for (int j = 0; j < dim1; j++)
    {
      forcevec1(i) -= stiffmat11(i,j) * ele1pos_(j);
    }
  }

  for (int i = 0; i < dim1; i++)
  {
    for (int j = 0; j < dim2; j++)
    {
      forcevec1(i) -= stiffmat12(i,j) * ele2pos_(j);
    }
  }

  // Evaluate meshtying forces for surface element
  for (int i = 0; i < dim2; i++)
  {
    for (int j = 0; j < dim1; j++)
    {
      forcevec2(i) -= stiffmat21(i,j) * ele1pos_(j);
    }
  }

  for (int i = 0; i < dim2; i++)
  {
    for (int j = 0; j < dim2; j++)
    {
      forcevec2(i) -= stiffmat22(i,j) * ele2pos_(j);
    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::FillElementStiffmtMeshtying(
    LINALG::SerialDenseMatrix& local_stiffmat11,
    LINALG::SerialDenseMatrix& local_stiffmat12,
    LINALG::SerialDenseMatrix& local_stiffmat21,
    LINALG::SerialDenseMatrix& local_stiffmat22,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22)
{

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim1; j++)
      stiffmat11(i,j) += local_stiffmat11(i,j);

  for (int i = 0; i < dim1; i++)
  {
    for (int j = 0; j < dim2; j++)
    {
      stiffmat12(i,j) += local_stiffmat12(i,j);
      stiffmat21(j,i) += local_stiffmat21(j,i);
    }
  }

  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < dim2; j++)
      stiffmat22(i,j) += local_stiffmat22(i,j);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::EvaluateBoundarySegmentation()
{


  /******************************************************************************************
   Initialize some data
   ******************************************************************************************/
  std::map<int, std::vector<bool> >* GaussPointProjectionTracker =
      EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointProjectionTracker();

  bool proj = false;
  bool startBool = (*GaussPointProjectionTracker)[ele1gid_][0];
  int change = -1;
  TYPEBTS xi1;
  TYPEBTS xi2;
  TYPEBTS xi3;
  TYPEBTS eta;
  const bool output = false;
  /******************************************************************************************
   End: Initialize some data
   ******************************************************************************************/


  // here, we check if we have a beam element that overlaps with the discretization
  // such that the first few Gauss points project, then some do not project and
  // then some project again, or vice versa --> we do not want such a case --> ERROR
  // one change of the projection is ok --> we can perform boundary segmentation
  for (int i = 0; i < n_gp_; i++)
  {
    if (change < 0 && (*GaussPointProjectionTracker)[ele1gid_][i] != startBool)
    {
      // we have a change in the projection flag
      change = i;
    }
    if (change > 0 && (*GaussPointProjectionTracker)[ele1gid_][i] == startBool)
    {
      // we have a second change in the projection flag --> ERROR
      std::cout << "Slave element " << ele1gid_ << " could not be projected" << std::endl;
      dserror("EXIT here");
    }
  }

  // later we will need the first Gauss point that does not project no matter in which
  // direction they are ordered

  double compareEta;
  // if the last few Gauss points do not project, i.e. if startBool == true
  // the first Gauss point that does not project is the one at position [change] in eta_
  if (startBool)
  {
    compareEta = eta_[change];

    if (output)
    {
      std::cout << "BOUNDARY SEGMENTATION for GP no. " << change << ", eta = " << compareEta
          << " of Beam element " << ele1gid_ << " and solid element " << Element2()->Id() << std::endl;
    }
  }
  // if the first few Gauss points do not project, i.e. if startBool == false
  // the first Gauss point that does not project is the one at position [change - 1] in eta_
  else
  {
    compareEta = eta_[change-1];

    if (output)
    {
      std::cout << "BOUNDARY SEGMENTATION for GP no. " << change-1 << ", eta = " << compareEta
          << " of Beam element " << ele1gid_ << std::endl;
    }
  }


  double curreta = 1000.0;

  // now we project all edges of the solid element in order to
  // find the correct beam coordinate
  bool onefound = false;

  if (numnodessol == 8 || numnodessol == 20 || numnodessol == 27)
  {
    for (int j = 0; j < 3; j++)
    // j = 0 --> xi1 fixed, j = 1 --> xi2 fixed, j = 2 --> xi3 fixed
    {
      // project edge xi1, xi2 or xi3 = 1.0
      SurfProjection(xi1, xi2, xi3, eta, j, 1.0, proj);
      if (proj)
      {
        onefound = true;
        // if the projection is valid, check if it is closer to the relevant Gauss
        // point compareEta
        if (fabs(eta - compareEta) < fabs(curreta - compareEta))
        {
          curreta = eta;
          if (output)
          {
            std::cout << "FOUND new point for beam element " << ele1gid_
                << " projecting on solid element " << Element2()->Id() << std::endl;
            std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
                << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
          }
        }
      }
      // project edge xi1, xi2 or xi3 = -1.0
      SurfProjection(xi1, xi2, xi3, eta, j, -1.0, proj);
      if (proj)
      {
        onefound = true;
        // if the projection is valid, check if it is closer to the relevant Gauss
        // point compareEta
        if (fabs(eta - compareEta) < fabs(curreta - compareEta))
        {
          curreta = eta;
          if (output)
          {
            std::cout << "FOUND new point for beam element " << ele1gid_
                << " projecting on solid element " << Element2()->Id() << std::endl;
            std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
                << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
          }
        }
      }
    }
  }
  else if (numnodessol == 4 || numnodessol == 10)
  {
    // project xi1 fixed at 0.0
    SurfProjection(xi1, xi2, xi3, eta, 0, 0.0, proj);
    if (proj)
    {
      onefound = true;
      // if the projection is valid, check if it is closer to the relevant Gauss
      // point compareEta
      if (fabs(eta - compareEta) < fabs(curreta - compareEta))
      {
        curreta = eta;
        if (output)
        {
          std::cout << "FOUND new point for beam element " << ele1gid_
              << " projecting on solid element " << Element2()->Id() << std::endl;
          std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
              << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
        }
      }
    }
    // project xi2 fixed at 0.0
    SurfProjection(xi1, xi2, xi3, eta, 1, 0.0, proj);
    if (proj)
    {
      onefound = true;
      // if the projection is valid, check if it is closer to the relevant Gauss
      // point compareEta
      if (fabs(eta - compareEta) < fabs(curreta - compareEta))
      {
        curreta = eta;
        if (output)
        {
          std::cout << "FOUND new point for beam element " << ele1gid_
              << " projecting on solid element " << Element2()->Id() << std::endl;
          std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
              << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
        }
      }
    }
    // project xi3 fixed at 0.0
    SurfProjection(xi1, xi2, xi3, eta, 2, 0.0, proj);
    if (proj)
    {
      onefound = true;
      // if the projection is valid, check if it is closer to the relevant Gauss
      // point compareEta
      if (fabs(eta - compareEta) < fabs(curreta - compareEta))
      {
        curreta = eta;
        if (output)
        {
          std::cout << "FOUND new point for beam element " << ele1gid_
              << " projecting on solid element " << Element2()->Id() << std::endl;
          std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
              << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
        }
      }
    }
    // project third surface of tetahedron xi1 + xi2 + xi3 = 1
    SurfProjection(xi1, xi2, xi3, eta, 3, 0.0, proj);
    if (proj)
    {
      onefound = true;
      // if the projection is valid, check if it is closer to the relevant Gauss
      // point compareEta
      if (fabs(eta - compareEta) < fabs(curreta - compareEta))
      {
        curreta = eta;
        if (output)
        {
          std::cout << "FOUND new point for beam element " << ele1gid_
              << " projecting on solid element " << Element2()->Id() << std::endl;
          std::cout << "Gauss point " << compareEta << ", Found eta = " << eta
              << ", xi1 = " << xi1 << ", xi2 = " << xi2 << ", xi3 = " << xi3 << std::endl;
        }
      }
    }
  }
  else
    dserror ("Only numnodessol = 8, 20, 27, 4 or 10 is valid for Boundary Segmentation so far");


  if (!onefound)
  {
    std::cout << "Boundary Segmentation could not be performed, no Edge Projection found"
        " for Beam Element " << ele1gid_ << std::endl;
    dserror ("EXIT HERE");
  }

  // compute the borders of the segment
  double etaA;
  double etaB;
  if (startBool)
  {
    // segment goes from [-1.0; curreta] because the first original Gauss point(s) could
    // be projected
    etaA = -1.0;
    etaB = curreta;
  }
  else
  {
    // segment goes from [curreta; 1.0] because the first original Gauss point could
    // not to be protected but the last one(s) could
    etaA = curreta;
    etaB = 1.0;
  }


  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GAUSSRULE);

  // get determinant for coordinate transformation
  double determinant = (etaB - etaA)/2.0;
  double dummy;
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    proj = false;
    // compute the coordinate trafo, weights and the jacobi*determinant
    eta_[i_gp] = (etaA + etaB)/2.0 + gaussPoints.qxg[i_gp][0] * (etaB - etaA)/2.0;
    wgp_[i_gp] = gaussPoints.qwgt[i_gp];
    GetBeamJacobi(eta_[i_gp], dummy);
    beamjacobi_[i_gp] = dummy*determinant;
    // compute the solid coords by projecting the new Gauss points

    Projection(xi1, xi2, xi3, eta_[i_gp], proj);
    if (proj)
    {
    // fill vectors
      xi1_[i_gp] = xi1;
      xi2_[i_gp] = xi2;
      xi3_[i_gp] = xi3;
      localGaussPointTracker_[i_gp] = true;

      if (output)
      {
        std::cout << "BOUNDARY SEGMENTATION FINISHED for GP no. " << i_gp
            << " of beam element " << ele1gid_ << std::endl;
        if (localGaussPointTracker_[i_gp] == true)
        {
          std::cout << "Gauss point " << eta_[i_gp] << " projects on Solid element "
              << Element2()->Id() << ", xi1 = " << xi1_[i_gp] << ", xi2 = " << xi2_[i_gp]
              << ", xi3 = " << xi3_[i_gp] << std::endl;
        }
      }

    }
    if (!proj)
    {
      std::cout << "In boundary segmentation: Gauss point " << eta_[i_gp]
          << " could not be projected" << std::endl;
      dserror("EXIT HERE");
    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::SurfProjection(
    TYPEBTS& xi1,
    TYPEBTS& xi2,
    TYPEBTS& xi3,
    TYPEBTS& eta,
    const int fixedPar,
    double fixedAt,
    bool& proj_allowed)
{
  proj_allowed = true;
  bool parallel = false;

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + XIETATOL;

  const bool output = false;

  // Vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2(true);                     // = N2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1(true);                 // = N2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2(true);                 // = N2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi3(true);                 // = N2,xi3

  // Coords and derivatives of beam and surface element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 1> x2(true);                                 // = x2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1(true);                             // = x2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2(true);                             // = x2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi3(true);                             // = x2,xi3


  // reset iteration variables
  eta = 0.0;
  if (numnodessol == 8 || numnodessol == 20 || numnodessol == 27)
  {
    if (fixedPar == 0) // xi1 fixed
    {
      xi1 = fixedAt;
      xi2 = 0.0;
      xi3 = 0.0;
    }
    else if (fixedPar == 1) // xi2 fixed
    {
      xi1 = 0.0;
      xi2 = fixedAt;
      xi3 = 0.0;
    }
    else if (fixedPar == 2) // xi3 fixed
    {
      xi1 = 0.0;
      xi2 = 0.0;
      xi3 = fixedAt;
    }
    else
      dserror ("only fixedPar = 0, fixedPar = 1 or fixedPar = 2 possible for hex elements");
  }
  else if (numnodessol == 4 || numnodessol == 10)
  {

    if (fixedPar == 0) // xi1 fixed
    {
      xi1 = fixedAt;
      xi2 = 0.25;
      xi3 = 0.25;
    }
    else if (fixedPar == 1) // xi2 fixed
    {
      xi1 = 0.25;
      xi2 = fixedAt;
      xi3 = 0.25;
    }
    else if (fixedPar == 2) // xi3 fixed
    {
      xi1 = 0.25;
      xi2 = 0.25;
      xi3 = fixedAt;
    }
    else if (fixedPar == 3) // fourth surface xi1 + xi2 + xi3 = 1 fixed
    {
      xi1 = 0.25;
      xi2 = 0.25;
      xi3 = 0.5;
    }
    else
      dserror ("only fixedPar = 0, fixedPar = 1, fixedPar = 2 or fixedPar = 3 possible for tet elements");
  }
  else
    dserror ("Only numnodessol = 8, 20, 27, 4 or 10 is valid for Boundary Segmentation so far");

  if (output)
  {
    std::cout << "Projection output:" << std::endl;
    std::cout << "Start parameters xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << std::endl;
  }

  // Initialize function f and Jacobian J for Newton iteration
  LINALG::TMatrix<TYPEBTS,3,1> f(true);
  LINALG::TMatrix<TYPEBTS,3,3> J(true);
  LINALG::TMatrix<TYPEBTS,3,3> Jinv(true);

  // Initial scalar residual (L2-norm of f)
  TYPEBTS residual = 0.0;

  // Local newton iteration
  // -----------------------------------------------------------------

  int iter;
  for (iter = 0; iter < BEAMCONTACTMAXITER; iter++)
  {
    // Update shape functions and their derivatives for beam and surface element
    GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
    GetSolShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi3, xi1, xi2, xi3);

    // Update coordinates and derivatives for beam and surface element
    ComputeBeamCoordsAndDerivsRef(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);
    ComputeSolCoordsAndDerivsRef(x2, x2_xi1, x2_xi2, x2_xi3, N2, N2_xi1, N2_xi2, N2_xi3);

    //std::cout << "n2: " << n2 << std::endl;
    // Evaluate f at current xi1, xi2, alpha
    f.Clear();
    for (int i = 0; i < 3; i++)
    {
      f(i) = x2(i) - r1(i);
    }
    // Compute scalar residuum
    residual = FADUTILS::VectorNorm<3>(f);

    // Reset matrices
    J.Clear();

    if (fixedPar == 0) // xi1 fixed
    {
      // df/dxi_2
      J(0,0) = x2_xi2(0);
      J(1,0) = x2_xi2(1);
      J(2,0) = x2_xi2(2);

      // df/dxi_3
      J(0,1) = x2_xi3(0);
      J(1,1) = x2_xi3(1);
      J(2,1) = x2_xi3(2);
    }
    else if (fixedPar == 1) // xi2 fixed
    {
      // df/dxi_1
      J(0,0) = x2_xi1(0);
      J(1,0) = x2_xi1(1);
      J(2,0) = x2_xi1(2);

      // df/dxi_3
      J(0,1) = x2_xi3(0);
      J(1,1) = x2_xi3(1);
      J(2,1) = x2_xi3(2);
    }
    else if (fixedPar == 2) // xi3 fixed
    {
      // df/dxi1
      J(0,0) = x2_xi1(0);
      J(1,0) = x2_xi1(1);
      J(2,0) = x2_xi1(2);

      // df/dxi_2
      J(0,1) = x2_xi2(0);
      J(1,1) = x2_xi2(1);
      J(2,1) = x2_xi2(2);
    }
    else if (fixedPar == 3) // xi3 fixed at xi3 = 1.0 - xi1 - xi2
    {
      // df/dxi1 = df/dxi1 - df/dxi3
      J(0,0) = x2_xi1(0) - x2_xi3(0);
      J(1,0) = x2_xi1(1) - x2_xi3(1);
      J(2,0) = x2_xi1(2) - x2_xi3(2);

      // df/dxi_2 = df/dxi2 - df/dxi3
      J(0,1) = x2_xi2(0) - x2_xi3(0);
      J(1,1) = x2_xi2(1) - x2_xi3(1);
      J(2,1) = x2_xi2(2) - x2_xi3(2);
    }

    // df/deta
    J(0,2) = -r1_eta(0);
    J(1,2) = -r1_eta(1);
    J(2,2) = -r1_eta(2);


    double jacdet = J.Determinant();

    // If det_J = 0 we assume, that the beam centerline and the surface edge are parallel.
    // These projection is not needed due the fact that the contact interval can also be
    // identified by two contact interval borders found with the GetContactLines method
    parallel = fabs(jacdet) < COLINEARTOL;
    if (!parallel)
      jacdet = J.Invert();

    // Check if the local Newton iteration has converged
    // If the start point fulfills the orthogonalty conditions (residual < BEAMCONTACTTOL), we also check if
    // the beam centerline and the surface edge are parallel. This is done by calculating det_J before checking
    // if the local Newton iteration has converged by fulfilling the condition residual < BEAMCONTACTTOL
    if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL && !parallel)
    {
      if (output)
      {
        std::cout << "Local Newton iteration converged after " << iter << " iterations" << std::endl;
        std::cout << "Found point at xi1: " << xi1 << ", xi2: " << xi2 << ", xi3: " << xi3
            << ", eta: " << eta << " with residual: " << residual << std::endl;
        std::cout << "r1:\n" << r1 << ", x2:\n" << x2 << std::endl;
      }
      // Local Newton iteration converged
      break;
    }
    else if (output && iter > 0)
    {
      std::cout << "New point at xi1: " << xi1 << ", xi2: " << xi2 << ", xi3: " << xi3
          << ", eta: " << eta << " with residual: " << residual <<  std::endl;
    }

    // Singular J
    if (parallel)
    {
      // Sort out
      if (output)
      {
        std::cout << "elementscolinear: det_J = " << jacdet << std::endl;
      }
      break;
    }
    // Regular J (inversion possible)

    if (fixedPar == 0)
    {
      xi2 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
      xi3 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    }
    else if (fixedPar == 1)
    {
      xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
      xi3 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    }
    else if (fixedPar == 2)
    {
      xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
      xi2 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    }
    else if (fixedPar == 3)
    {
      xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
      xi2 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
      xi3 = 1.0 - xi1 - xi2;
    }

    eta += -J(2, 0) * f(0) - J(2, 1) * f(1) - J(2, 2) * f(2);

  }
  // -----------------------------------------------------------------
  // End: Local Newton iteration


  // Local Newton iteration unconverged after BEAMCONTACTMAXITER
  if (residual > BEAMCONTACTTOL || parallel)
  {

    xi1 = 1e+12;
    xi2 = 1e+12;
    xi3 = 1e+12;
    eta = 1e+12;

    if (output)
      std::cout << "Local Newton iteration unconverged (!) after " << iter + 1 << " iterations" << std::endl;
  }
  // hex8, hex27, hex20 --> xi1, xi2, xi3 has to be in [-1; 1] for valid projection
  if (numnodessol == 8 || numnodessol == 20 || numnodessol == 27)
  {
    if (fabs(xi1) > limit || fabs(xi2) > limit || fabs(xi3) > limit || fabs(eta) > limit)
    {
      proj_allowed = false;
      if (output)
        std::cout << "Projection not allowed" << std::endl;
    }
  }
  // tet4, tet10 --> xi1 > 0; xi2 > 0; xi3 > 0 and xi1 + xi2 + xi3 < 1 for valid projection
  else if (numnodessol == 4 || numnodessol == 10)
  {
    if (xi1 < -XIETATOL || xi2 < -XIETATOL || xi3 < -XIETATOL
        || xi1 + xi2 + xi3 > 1.0 + XIETATOL || fabs(eta) > limit)
    {
      proj_allowed = false;
      if (output)
        std::cout << "Projection not allowed" << std::endl;
    }
  }
  else
    dserror("numnodessol = 8, 20, 27, 4 or 10 is valid for Boundary Segmentation so far");

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Projection(
    TYPEBTS& xi1,
    TYPEBTS& xi2,
    TYPEBTS& xi3,
    TYPEBTS& eta,
    bool& proj_allowed)
{
  proj_allowed = true;
  bool parallel = false;

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + XIETATOL;

  const bool output = false;

  // Vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2(true);                     // = N2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1(true);                 // = N2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2(true);                 // = N2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi3(true);                 // = N2,xi3

  // Coords and derivatives of beam and surface element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 1> x2(true);                                 // = x2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1(true);                             // = x2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2(true);                             // = x2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi3(true);                             // = x2,xi3

  if (numnodessol == 8 || numnodessol == 27 || numnodessol == 20)
  {
    xi1 = 0.0;
    xi2 = 0.0;
    xi3 = 0.0;
  }
  else if (numnodessol == 4 || numnodessol == 10)
  {
    xi1 = 0.25;
    xi2 = 0.25;
    xi3 = 0.25;
  }
  else
    dserror("Only numnodessol = 8, 20, 27, 4, 10 is valid");

  if (output)
  {
    std::cout << "Projection output:" << std::endl;
    std::cout << "Start parameters xi1: " << xi1 << ", xi2: " << xi2 << ", xi3: " << xi3 << std::endl;
  }

  // Initialize function f and Jacobian J for Newton iteration
  LINALG::TMatrix<TYPEBTS,3,1> f(true);
  LINALG::TMatrix<TYPEBTS,3,3> J(true);
  LINALG::TMatrix<TYPEBTS,3,3> Jinv(true);

  // Initial scalar residual (L2-norm of f)
  TYPEBTS residual = 0.0;

  // Local newton iteration
  // -----------------------------------------------------------------

  int iter;

  for (iter = 0; iter < BEAMCONTACTMAXITER; iter++)
  {

    // Update shape functions and their derivatives for beam and solid element
    GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
    GetSolShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi3, xi1, xi2, xi3);

    // Update coordinates and derivatives for beam and solid element
    ComputeBeamCoordsAndDerivsRef(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);
    ComputeSolCoordsAndDerivsRef(x2, x2_xi1, x2_xi2, x2_xi3, N2, N2_xi1, N2_xi2, N2_xi3);

    // Evaluate f at current xi1, xi2, alpha
    f.Clear();
    for (int i = 0; i < 3; i++)
    {
      f(i) = x2(i) - r1(i);
    }

    // Compute scalar residuum
    residual = FADUTILS::VectorNorm<3>(f);

    // Reset matrices
    J.Clear();

    // df/dxi_1
    J(0,0) = x2_xi1(0);
    J(1,0) = x2_xi1(1);
    J(2,0) = x2_xi1(2);

    // df/dxi_2
    J(0,1) = x2_xi2(0);
    J(1,1) = x2_xi2(1);
    J(2,1) = x2_xi2(2);

    // df/dxi_3
    J(0,2) = x2_xi3(0);
    J(1,2) = x2_xi3(1);
    J(2,2) = x2_xi3(2);


    double jacdet = J.Determinant();

    // If det_J = 0 we assume, that the beam centerline and the surface edge are parallel.
    // These projection is not needed due the fact that the contact interval can also be
    // identified by two contact interval borders found with the GetContactLines method
    parallel = fabs(jacdet) < COLINEARTOL;
    if (!parallel)
      jacdet = J.Invert();

    // Check if the local Newton iteration has converged
    // If the start point fulfills the orthogonalty conditions (residual < BEAMCONTACTTOL), we also check if
    // the beam centerline and the surface edge are parallel. This is done by calculating det_J before checking
    // if the local Newton iteration has converged by fulfilling the condition residual < BEAMCONTACTTOL
    if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL && !parallel)
    {
      if (output)
      {
        std::cout << "Local Newton iteration converged after " << iter << " iterations" << std::endl;
        std::cout << "Found point at xi1: " << xi1 << ", xi2: " << xi2 << ", xi3: " << xi3
            << " with residual: " << residual << std::endl;
        std::cout << "r1:\n" << r1 << ", x2:\n" << x2 << std::endl;
      }
      // Local Newton iteration converged
      break;
    }
    else if (output && iter > 0)
    {
      std::cout << "New point at xi1: " << xi1 << ", xi2: " << xi2 << ", xi3: " << xi3
          << " with residual: " << residual <<  std::endl;
    }

    // Singular J
    if (parallel)
    {
      // Sort out
      if (output)
      {
        std::cout << "elementscolinear: det_J = " << jacdet << std::endl;
      }
      break;
    }
    // Regular J (inversion possible)

    xi1 += -J(0, 0) * f(0) - J(0, 1) * f(1) - J(0, 2) * f(2);
    xi2 += -J(1, 0) * f(0) - J(1, 1) * f(1) - J(1, 2) * f(2);
    xi3 += -J(2, 0) * f(0) - J(2, 1) * f(1) - J(2, 2) * f(2);

  }
  // -----------------------------------------------------------------
  // End: Local Newton iteration


  // Local Newton iteration unconverged after BEAMCONTACTMAXITER
  if (residual > BEAMCONTACTTOL || parallel)
  {

    xi1 = 1e+12;
    xi2 = 1e+12;
    xi3 = 1e+12;

    if (output)
      std::cout << "Local Newton iteration unconverged (!) after " << iter + 1 << " iterations" << std::endl;
  }
  // HEX8, HEX27, HEX20 --> xi1, xi2, xi3 has to be in [-1; 1] for valid projection
  if (numnodessol == 8 || numnodessol == 27 || numnodessol == 20)
  {
    if (fabs(xi1) > limit || fabs(xi2) > limit || fabs(xi3) > limit)
    {
      proj_allowed = false;
      if (output)
        std::cout << "Projection not allowed" << std::endl;
    }
    else
    {
      if (output)
        std::cout << "Projection allowed" << std::endl;
    }
  }
  // TET4, TET10 --> xi1 > 0; xi2 > 0; xi3 > 0 and xi1 + xi2 + xi3 < 1 for valid projection
  else if (numnodessol == 4 || numnodessol == 10)
  {
    if (xi1 < -XIETATOL || xi2 < -XIETATOL || xi3 < -XIETATOL || xi1 + xi2 + xi3 > 1.0 + XIETATOL)
    {
      proj_allowed = false;
      if (output)
        std::cout << "Projection not allowed" << std::endl;
    }
    else
    {
      if (output)
        std::cout << "Projection allowed" << std::endl;
    }
  }
  else
    dserror("Only numnodessol = 8, 20, 27, 4 or 10 is valid");

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::GetSolShapeFunctions(
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi3,
    const TYPEBTS& xi1,
    const TYPEBTS& xi2,
    const TYPEBTS& xi3)
{
  // Clear shape functions and derivatives
  N.Clear();
  N_xi1.Clear();
  N_xi2.Clear();
  N_xi3.Clear();


  LINALG::TMatrix<TYPEBTS, 1, numnodessol> N_i(true);
  LINALG::TMatrix<TYPEBTS, 3, numnodessol> N_i_xi(true);

  DRT::UTILS::shape_function_3D(N_i, xi1, xi2, xi3, Element2()->Shape());
  DRT::UTILS::shape_function_3D_deriv1(N_i_xi, xi1, xi2, xi3, Element2()->Shape());

  // Assemble the individual shape functions in matrices, such that: r = N * d, r_xi = N_xi * d, ...
  AssembleSolShapefunctions(N_i, N_i_xi, N, N_xi1, N_xi2, N_xi3);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::AssembleSolShapefunctions(
    const LINALG::TMatrix<TYPEBTS,1,numnodessol>& N_i,
    const LINALG::TMatrix<TYPEBTS,3,numnodessol>& N_i_xi,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi3)
{
  // assembly_N is just an array to help assemble the matrices of the shape functions
  // it determines, which shape function is used in which column of N
  int assembly_N[3][3*numnodessol];

  // Initialize to zero
  for (unsigned int i = 0; i < 3*numnodessol; i++)
    for (int j = 0; j < 3; j++)
      assembly_N[j][i] = 0.0;


  // Set number of shape functions for each 3x3 block:
  // e.g. quad4 surface element (numnodesol = 4)
  // int assembly_N[3][12] = {{1,0,0,2,0,0,3,0,0,4,0,0},
  //                          {0,1,0,0,2,0,0,3,0,0,4,0},
  //                          {0,0,1,0,0,2,0,0,3,0,0,4}};

  for (unsigned int i = 0; i < numnodessol; i++)
  {
    assembly_N[0][3*i] = i + 1;
    assembly_N[1][3*i + 1] = i + 1;
    assembly_N[2][3*i + 2] = i + 1;
  }

  // Assemble the matrices of the shape functions
  for (unsigned int i = 0; i < 3*numnodessol; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (assembly_N[j][i] == 0)
      {
        N(j,i) = 0;
        N_xi1(j,i) = 0;
        N_xi2(j,i) = 0;
        N_xi3(j,i) = 0;
      }
      else
      {
        int k = assembly_N[j][i] - 1;
        N(j,i) = N_i(k);
        N_xi1(j,i) = N_i_xi(0,k);
        N_xi2(j,i) = N_i_xi(1,k);
        N_xi3(j,i) = N_i_xi(2,k);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::ComputeSolCoordsAndDerivsRef(
    LINALG::TMatrix<TYPEBTS,3,1>& r,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi1,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi3,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi3)
{
  r.Clear();
  r_xi1.Clear();
  r_xi2.Clear();
  r_xi3.Clear();

  // Compute output variable
  for (int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3*numnodessol; j++)
    {
      r(i) += N(i,j) * ele2posref_(j);
      r_xi1(i) += N_xi1(i,j) * ele2posref_(j);
      r_xi2(i) += N_xi2(i,j) * ele2posref_(j);
      r_xi3(i) += N_xi3(i,j) * ele2posref_(j);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::ComputeSolCoordsAndDerivsMom(
    LINALG::TMatrix<TYPEBTS,3,1>& r,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi1,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi3,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi3)
{
  r.Clear();
  r_xi1.Clear();
  r_xi2.Clear();
  r_xi3.Clear();

  // Compute output variable
  for (int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3*numnodessol; j++)
    {
      r(i) += N(i,j) * ele2pos_(j);
      r_xi1(i) += N_xi1(i,j) * ele2pos_(j);
      r_xi2(i) += N_xi2(i,j) * ele2pos_(j);
      r_xi3(i) += N_xi3(i,j) * ele2pos_(j);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::GetBeamJacobi(
    const TYPEBTS& eta,
    TYPEBTS& jac)
{
  // Vectors for beam shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  // Coords and derivatives of beam element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  // Get Shape Functions type
  GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);

  // Update coordinates and derivatives for beam element
  ComputeBeamCoordsAndDerivsRef(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);

  // Calculate the Jacobian
  jac = FADUTILS::VectorNorm<3>(r1_eta);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::GetBeamShapeFunctions(
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N_eta,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N_etaeta,
    const TYPEBTS& eta)
{
  // Clear shape functions and derivatives
  N.Clear();
  N_eta.Clear();
  N_etaeta.Clear();

  // Get discretization type
  const DRT::Element::DiscretizationType distype = Element1()->Shape();
  //std::cout << "distype " << distype << std::endl;

  LINALG::TMatrix<TYPEBTS, 1, numnodes*numnodalvalues> N_i(true);
  LINALG::TMatrix<TYPEBTS, 1, numnodes*numnodalvalues> N_i_eta(true);
  LINALG::TMatrix<TYPEBTS, 1, numnodes*numnodalvalues> N_i_etaeta(true);

  if (numnodalvalues == 1)
  {
    // Get values and derivatives of shape functions
    DRT::UTILS::shape_function_1D(N_i, eta, distype);
    DRT::UTILS::shape_function_1D_deriv1(N_i_eta, eta, distype);
    DRT::UTILS::shape_function_1D_deriv2(N_i_etaeta, eta, distype);
  }
  else if (numnodalvalues == 2)
  {

    double length = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1()))->RefLength();

    const DRT::Element::DiscretizationType distype1herm = DRT::Element::line2;

    // Get values and derivatives of shape functions
    DRT::UTILS::shape_function_hermite_1D(N_i, eta, length, distype1herm);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_eta, eta, length, distype1herm);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_etaeta, eta, length, distype1herm);
  }
  else
    dserror("Only beam elements with one (nodal positions) or two "
        "(nodal positions + nodal tangents) values are valid!");

  // Assemble the individual shape functions in matrices, such that: r = N * d,
  // r_eta = N_eta * d, r_etaeta = N_etaeta * d
  AssembleBeamShapefunctions(N_i, N_i_eta, N_i_etaeta, N, N_eta, N_etaeta);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::AssembleBeamShapefunctions(
    const LINALG::TMatrix<TYPEBTS,1,numnodes*numnodalvalues>& N_i,
    const LINALG::TMatrix<TYPEBTS,1,numnodes*numnodalvalues>& N_i_eta,
    const LINALG::TMatrix<TYPEBTS,1,numnodes*numnodalvalues>& N_i_etaeta,
    LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_eta,
    LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_etaeta)
{
  // assembly_N is just an array to help assemble the matrices of the shape functions
  // it determines, which shape function is used in which column of N
  int assembly_N[3][3*numnodes*numnodalvalues];

  // Initialize to zero
  for (unsigned int i = 0;i < 3*numnodes*numnodalvalues; i++)
    for (int j = 0; j < 3; j++)
      assembly_N[j][i] = 0.0;

  // Set number of shape functions for each 3x3 block:
  // e.g. second order Reissner beam (numnodes = 3, numnodalvalues = 1)
  // int assembly_N[3][9] = {{1,0,0,2,0,0,3,0,0},
  //                         {0,1,0,0,2,0,0,3,0},
  //                         {0,0,1,0,0,2,0,0,3}};
  //
  // e.g. Kirchhoff beam (numnodes = 2, numnodalvalues = 2)
  // int assembly_N[3][12] = {{1,0,0,2,0,0,3,0,0,4,0,0},
  //                          {0,1,0,0,2,0,0,3,0,0,4,0},
  //                          {0,0,1,0,0,2,0,0,3,0,0,4}};

  for (unsigned int i = 0; i < numnodes*numnodalvalues; i++)
  {
    assembly_N[0][3*i] = i + 1;
    assembly_N[1][3*i + 1] = i + 1;
    assembly_N[2][3*i + 2] = i + 1;
  }

  // Assemble the matrices of the shape functions
  for (unsigned int i = 0; i < 3*numnodes*numnodalvalues; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (assembly_N[j][i]==0)
      {
        N(j,i) = 0;
        N_eta(j,i) = 0;
        N_etaeta(j,i) = 0;
      }
      else
      {
        int k = assembly_N[j][i]-1;
        N(j,i) = N_i(k);
        N_eta(j,i) = N_i_eta(k);
        N_etaeta(j,i) = N_i_etaeta(k);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::ComputeBeamCoordsAndDerivsRef(
    LINALG::TMatrix<TYPEBTS,3,1>& r,
    LINALG::TMatrix<TYPEBTS,3,1>& r_eta,
    LINALG::TMatrix<TYPEBTS,3,1>& r_etaeta,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_eta,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_etaeta)
{
  r.Clear();
  r_eta.Clear();
  r_etaeta.Clear();

  // Compute output variable
  for (int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3*numnodes*numnodalvalues; j++)
    {
      r(i) += N(i,j) * ele1posref_(j);
      r_eta(i) += N_eta(i,j) * ele1posref_(j);
      r_etaeta(i) += N_etaeta(i,j) * ele1posref_(j);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::ComputeBeamCoordsAndDerivsMom(
    LINALG::TMatrix<TYPEBTS,3,1>& r,
    LINALG::TMatrix<TYPEBTS,3,1>& r_eta,
    LINALG::TMatrix<TYPEBTS,3,1>& r_etaeta,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_eta,
    const LINALG::TMatrix<TYPEBTS,3,3*numnodes*numnodalvalues>& N_etaeta)
{
  r.Clear();
  r_eta.Clear();
  r_etaeta.Clear();

  // Compute output variable
  for (int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3*numnodes*numnodalvalues; j++)
    {
      r(i) += N(i,j) * ele1pos_(j);
      r_eta(i) += N_eta(i,j) * ele1pos_(j);
      r_etaeta(i) += N_etaeta(i,j) * ele1pos_(j);
    }
  }

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // beam element
  for (unsigned int i = 0; i < 3*numnodalvalues; i++)
  {
    for (unsigned int j = 0; j < numnodes; j++)
    {
      ele1pos_(3*numnodalvalues*j + i) = beam_centerline_dofvec[3*numnodalvalues*j + i];
    }
  }

  // solid element
  for (int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < numnodessol; j++)
    {
      ele2pos_(3*j + i) = solid_nodal_dofvec[3*j + i];
    }
  }

  // reset the GaussPointEvaluateTracker for the next iteration
  std::map<int, std::vector<bool> >* GaussPointEvaluateTracker
    = EvaluationData()->BeamToSolidVolumeMeshtyingEvaluationData()->BTSVOLMTGaussPointEvaluateTracker();

  std::vector<bool> temp;
  temp.resize(n_gp_, false);

  (*GaussPointEvaluateTracker)[ele1gid_] = temp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::Print(
    std::ostream& out) const
{

  CheckInitSetup();

  // Print some general information: Element IDs and dofvecs
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToSolidVolumeMeshtyingPair"
      << "\nBeam EleGID:  " << ele1gid_
      << "\nSolid EleGID: " << Element2()->Id();

  out << "\n\nele1 dofvec: " << ele1pos_;
  out << "\nele2 dofvec: " << ele2pos_;
  out << "\n";

  // Print status of each gauss point for that pair
  out << "Gauss Point Status:";
  int numberOfValidProjections = 0;
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    out << "\nGauss-Point " << i_gp;
    if (localGaussPointTracker_[i_gp] == 1)
    {
      out << " active";
      numberOfValidProjections++;
    }
    else
    {
      out << " inactive";
    }
  }
  out << "\n";

  out << "------------------------------------------------------------------------\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodessol, unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<numnodessol, numnodes, numnodalvalues>::PrintSummaryOneLinePerActiveSegmentPair(
    std::ostream& out) const
{

  CheckInitSetup();

  // count number of valid projections
  int numberOfValidProjections = 0;
  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    if (localGaussPointTracker_[i_gp] == 1)
    {
      numberOfValidProjections++;
    }
  }

  for (int i_gp = 0; i_gp < n_gp_; i_gp++)
  {
    if ( eta_[i_gp] != 0 )
    {
      out << "\nGP "<< i_gp << "  "<<  ele1gid_ << " (beam)"
          << "  |  " << Element2()->Id() << " (solid) "
          << "  |  " << "VMT"
          << " | xi = (" << xi1_[i_gp] << ", " << xi2_[i_gp] << ", " << xi3_[i_gp] << ")"
          << " | eta = " << eta_[i_gp];
    }
  }

}


//explicit template instantiations
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<8,2,2>;     // Hermite beam element, hex8 solid element
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<20,2,2>;     // Hermite beam element, hex20 solid element
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<27,2,2>;     // Hermite beam element, hex27 solid element
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<4,2,2>;     // Hermite beam element, tet4 solid element
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair<10,2,2>;     // Hermite beam element, tet10 solid element
