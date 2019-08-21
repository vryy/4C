/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief class to handle contact between a 3D beam element and a rigid sphere

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_to_sphere_contact_pair.H"

#include "beam3contact_utils.H"
#include "beam3contact_defines.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"

// Todo check and get rid of outdated header inclusions
#include "../drt_inpar/inpar_beamcontact.H"

#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "beam_contact_params.H"
#include "beam_to_sphere_contact_params.H"



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::BeamToSphereContactPair()
    : beam_element_(NULL), sphere_element_(NULL), contactflag_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::Setup()
{
  CheckInit();

  // call setup of base class first
  BeamContactPair::Setup();


  ele1pos_.Clear();
  ele2pos_.Clear();

  fc1_.Clear();
  fc2_.Clear();

  gap_ = 0.0;

  // initialize class variables for contact point coordinates
  x1_.Clear();
  x2_.Clear();
  normal_.Clear();


  // cast first element to Beam3Base
  beam_element_ = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1());

  if (beam_element_ == NULL)
    dserror(
        "cast to Beam3Base failed! first element in BeamToSphereContactPair pair"
        " must be a beam element!");

  // get radius and stress-free reference length of beam element
  radius1_ = BEAMCONTACT::CalcEleRadius(beam_element_);
  beamele_reflength_ = beam_element_->RefLength();

  // cast second element to RigidSphere
  sphere_element_ = dynamic_cast<const DRT::ELEMENTS::Rigidsphere*>(Element2());

  if (sphere_element_ == NULL)
    dserror(
        "cast to Rigidsphere failed! second element in BeamToSphereContactPair pair"
        " must be a Rigidsphere element!");

  // get radius of sphere element
  radius2_ = sphere_element_->Radius();


  nodalcontactflag_.assign(2, false);

  issetup_ = true;


  // *********************** DEBUG ************************************************
  //  std::cout << "\nSuccessful Creation&Init&Setup of";
  //  Print(std::cout);
  // ******************* END DEBUG ************************************************
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::PreEvaluate()
{
  // do nothing
  return;
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
bool BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  unsigned int dim1 = 3 * numnodes * numnodalvalues;
  unsigned int dim2 = 3;

  // resize and initialize variables to zero
  if (forcevec1 != NULL) forcevec1->Size(dim1);
  if (forcevec2 != NULL) forcevec2->Size(dim2);

  if (stiffmat11 != NULL) stiffmat11->Shape(dim1, dim1);
  if (stiffmat12 != NULL) stiffmat12->Shape(dim1, dim2);
  if (stiffmat21 != NULL) stiffmat21->Shape(dim2, dim1);
  if (stiffmat22 != NULL) stiffmat22->Shape(dim2, dim2);

  //**********************************************************************
  // Evaluation of contact forces and stiffness
  //**********************************************************************
  // (1) Closest Point Projection
  // (2) Compute some auxiliary quantities
  //     -> normal vector, gap, shape functions, contact flag,
  //     -> linearizations of all geometric quantities
  // (3) Compute contact forces and stiffness
  //     -> stiffness terms are directly assembled to global matrix
  //     -> contact forces are only returned as global vector
  // (4) Perform some finite difference checks
  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
  // set class variable for status of gapfunction

  // Note: element1 is always the beam, element2 is the sphere

  ClosestPointProjection();

  // vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i;
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i_xi;
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i_xixi;

  // coords and derivatives of the two contact points
  LINALG::TMatrix<TYPE, 3, 1> x1;    // = x1
  LINALG::TMatrix<TYPE, 3, 1> x2;    // = x2
  LINALG::TMatrix<TYPE, 3, 1> dx1;   // = x1,xi
  LINALG::TMatrix<TYPE, 3, 1> ddx1;  // = x1,xixi

  // initialize
  LINALG::TMatrix<TYPE, 3, 1> normal;
  TYPE gap = 0.0;
  TYPE norm = 0.0;

  bool validclosestpointprojection = true;

  if (std::abs(FADUTILS::CastToDouble(xicontact_)) <
      (1.0 + XIETATOL))  // ToDo when to reset nodalcontactflag_?
  {
    nodalcontactflag_[0] = false;
    nodalcontactflag_[1] = false;
  }
  else
  {
    contactflag_ = false;
    validclosestpointprojection = false;
  }

  if (validclosestpointprojection)
  {
    //**********************************************************************
    // (1) Compute some auxiliary quantities
    //**********************************************************************

    // call function to fill variables for shape functions and their derivatives
    GetShapeFunctions(N1_i, N1_i_xi, N1_i_xixi, xicontact_);

    // call function to fill variables with coords and derivs of the contact point
    ComputeCoordsAndDerivs(x1, x2, dx1, ddx1, N1_i, N1_i_xi, N1_i_xixi);

    // call function to compute scaled normal and gap in possible contact point
    ComputeNormal(normal, gap, norm, x1, x2);

    // call function to evaluate contact status
    CheckAndSetContactStatus();

    //**********************************************************************
    // (2) Compute contact forces and stiffness
    //**********************************************************************

    // call function to evaluate contact forces
    if (forcevec1 != NULL and forcevec2 != NULL)
    {
      EvaluateFcContact(*forcevec1, *forcevec2,
          Params()->BeamToSphereContactParams()->BeamToSpherePenaltyParam(), gap, normal, N1_i,
          contactflag_);
    }

    // call function to evaluate contact stiffness
    if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
    {
      EvaluateStiffcContact(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22,
          Params()->BeamToSphereContactParams()->BeamToSpherePenaltyParam(), gap, normal, norm, x1,
          x2, dx1, ddx1, N1_i, N1_i_xi, N1_i_xixi, contactflag_);
    }
  }
  else  // Fixme do this only for the nodes (physical) end points of beams in case of C1-continuous
        // beam centerline representation
  {
    if (numnodes != 2)
      dserror(
          "BeamToSphereContactPair: Check for nodal contact only implemented for 2-noded beam "
          "elements!");

    //************************************************************************
    // NEW: NODAL-BEAM-TO-SPHERE CONTACT
    //************************************************************************
    // distance of beam nodes and center of sphere is computed
    // -> penalty contribution if distance is smaller than sum of radii
    // (point-to-point contact: xi=+-1, hence no contribution from linearization of xi)
    // this also has consequences for AUTOMATICDIFF: do not use xicontact_ here!
    TYPE XiContact = 0.0;
    for (unsigned int inode = 0; inode < 2; ++inode)
    {
      /* if inode=0 is active contact (gap<0), we do NOT want to run the loop for inode=1
       * because the case validclosestpointprojection=false and nodalcontactflag_[0]=true and
       * nodalcontactflag_[1]=true is highly unlikely and class variables like gap_ are overwritten
       * if we run this loop again this leads to wrong output for gap_ of active BTSPH pairs in
       * Beam3contactmanager */
      if (nodalcontactflag_[0] != true)
      {
        // Set XiContact to +-1.0
        switch (inode)
        {
          case 0:
            XiContact = -1.0;  // first node
            break;
          case 1:
            XiContact = 1.0;  // second node
            break;
          default:
            dserror("This must not happen!");
            break;
        }

        // Now do exactly the same as for "normal" contact based on closest point projection
        // (except for contribution from linearization of xi)

        // TODO: make more efficient: use knowledge that contact point is a node -> matrix of shape
        // functions, etc. simplifies

        //**********************************************************************
        // (1) Compute some auxiliary quantities
        //**********************************************************************

        // call function to fill variables for shape functions and their derivatives
        GetShapeFunctions(N1_i, N1_i_xi, N1_i_xixi, XiContact);

        // call function to fill variables with coords and derivs of the contact point
        ComputeCoordsAndDerivs(x1, x2, dx1, ddx1, N1_i, N1_i_xi, N1_i_xixi);

        // call function to compute scaled normal and gap in possible contact point
        ComputeNormal(normal, gap, norm, x1, x2);

        // evaluate nodal contact status
        if (gap_ < 0)
        {
          nodalcontactflag_[inode] = true;
        }
        else
          nodalcontactflag_[inode] = false;

        //**********************************************************************
        // (2) Compute contact forces and stiffness
        //**********************************************************************

        // call function to evaluate contact forces
        if (forcevec1 != NULL and forcevec2 != NULL)
        {
          EvaluateFcContact(*forcevec1, *forcevec2,
              Params()->BeamToSphereContactParams()->BeamToSpherePenaltyParam(), gap, normal, N1_i,
              nodalcontactflag_[inode]);
        }

        // call function to evaluate contact stiffness
        if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
        {
          EvaluateStiffcContact(*stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22,
              Params()->BeamToSphereContactParams()->BeamToSpherePenaltyParam(), gap, normal, norm,
              x1, x2, dx1, ddx1, N1_i, N1_i_xi, N1_i_xixi, nodalcontactflag_[inode], false);
        }
      }
    }
  }

  return GetContactFlag();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ClosestPointProjection()
{
  // local variable for beam element coordinate
  TYPE eta = 0.0;

  // vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i;
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i_xi;
  LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues> N1_i_xixi;

  // coords and derivatives of the two contact points
  LINALG::TMatrix<TYPE, 3, 1> x1;       // = x1
  LINALG::TMatrix<TYPE, 3, 1> x2;       // = x2
  LINALG::TMatrix<TYPE, 3, 1> dx1;      // = x1,xi
  LINALG::TMatrix<TYPE, 3, 1> ddx1;     // = x1,xixi
  LINALG::TMatrix<TYPE, 3, 1> delta_x;  // = x1 - x2

  // initialize function f and Jacobian df for Newton iteration
  TYPE f;
  TYPE df;

  // initial scalar residual (L2-norm of f)
  double residual = 0.0;

  int iter = 0;

  //**********************************************************************
  // local Newton iteration
  //**********************************************************************
  for (int i = 0; i < BEAMCONTACTMAXITER; ++i)
  {
    iter++;

    // reset shape function variables to zero
    N1_i.Clear();
    N1_i_xi.Clear();
    N1_i_xixi.Clear();

    // update shape functions and their derivatives
    GetShapeFunctions(N1_i, N1_i_xi, N1_i_xixi, eta);

    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(x1, x2, dx1, ddx1, N1_i, N1_i_xi, N1_i_xixi);

    // compute delta_x = x1-x2
    for (unsigned int j = 0; j < 3; ++j) delta_x(j) = x1(j) - x2(j);

    // compute norm of difference vector to scale the equations
    // (this yields better conditioning)
    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type
    // double since this factor is needed for a pure scaling of the nonlinear CCP and has not to be
    // linearized!
    double norm_delta_x = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(delta_x));

    // the closer the beams get, the smaller is norm
    // norm is not allowed to be too small, else numerical problems occur
    if (norm_delta_x < NORMTOL)
    {
      dserror("Contact points x1 and x2 are identical. Choose smaller time step!");
    }

    // evaluate f at current eta
    EvaluateOrthogonalityCondition(f, delta_x, norm_delta_x, dx1);

    // compute the scalar residuum
    residual = abs(FADUTILS::CastToDouble(f));

    // check if Newton iteration has converged
    if (residual < BEAMCONTACTTOL) break;

    // evaluate Jacobian of f at current eta
    EvaluateLinOrthogonalityCondition(df, delta_x, norm_delta_x, dx1, ddx1);

    // singular df
    if (abs(df) < COLINEARTOL)
    {
      dserror("No solution for Closest Point Projection!");
    }
    // regular df (inversion possible)
    else
    {
      // update element coordinates of contact point
      eta += -f / df;
    }
  }
  //**********************************************************************

  // Newton iteration unconverged
  if (residual > BEAMCONTACTTOL)
  {
    eta = 1e+12;
  }

  // store final result and return
  xicontact_ = eta;

#ifdef AUTOMATICDIFF
  // Set xi1_ as (additional) primary variable for automatic differentiation
  // The dependence between the infinitesimal changes delta xi1_ and the
  // the increments of the primary displacement variables delta disp have to be given explicitly,
  // since no explicit relation between the finite quantities xi1_ and disp exists. The latter would
  // have been necessary if the full linearization had to be computed directly with Sacado!!!

  // The 3*numnodes*numnodalvalues+3 primary DoFs are the components of the nodal positions /
  // tangents. The additional degree of freedom (+1) represents the dependency on the beam parameter
  // coordinate xi, which is necessary in beam contact.
  xicontact_.diff(
      (3 * numnodes * numnodalvalues + 3 + 1) - 1, 3 * numnodes * numnodalvalues + 3 + 1);
#endif

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes,
    numnodalvalues>::EvaluateOrthogonalityCondition(TYPE& f,
    const LINALG::TMatrix<TYPE, 3, 1>& delta_x, const double norm_delta_x,
    const LINALG::TMatrix<TYPE, 3, 1>& dx1)
{
  // reset f
  f = 0.0;

  // evaluate f
  // see Wriggers, Computational Contact Mechanics, equation (12.5)

  for (unsigned int i = 0; i < 3; i++)
  {
    f += delta_x(i) * dx1(i) / norm_delta_x;
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes,
    numnodalvalues>::EvaluateLinOrthogonalityCondition(TYPE& df,
    LINALG::TMatrix<TYPE, 3, 1>& delta_x, const double norm_delta_x,
    const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1)

{
  // reset df
  df = 0;

  // evaluate df
  // see Wriggers, Computational Contact Mechanics, equation (12.7)
  for (unsigned int i = 0; i < 3; i++)
  {
    df += (dx1(i) * dx1(i) + delta_x(i) * ddx1(i)) / norm_delta_x;
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::CheckAndSetContactStatus()
{
  // check contact condition
  contactflag_ = (gap_ < 0) ? true : false;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::EvaluateFcContact(
    LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseVector& forcevec2, const double& pp,
    const TYPE& gap, const LINALG::TMatrix<TYPE, 3, 1>& normal,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const bool contactactive)
{
  fc1_.Clear();
  fc2_.Clear();

  // get dimensions for vectors fc1 and fc2
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3;

  // flag indicating assembly
  bool DoNotAssemble = false;

  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (contactactive)
  {
    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    for (unsigned int i = 0; i < numnodes * numnodalvalues; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        fc1_(3 * i + j) = (-pp * gap) * normal(j) * N1_i(i);
      }
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    for (unsigned int j = 0; j < 3; ++j)
    {
      fc2_(j) = pp * gap * normal(j);
    }
  }
  //**********************************************************************
  // no forces for inactive pairs
  //**********************************************************************
  else
  {
    // set flag to avoid assembly
    DoNotAssemble = true;
  }

  //**********************************************************************
  // assemble contact forces
  //**********************************************************************
  if (!DoNotAssemble)
  {
    for (unsigned int i = 0; i < dim1; ++i) forcevec1(i) += FADUTILS::CastToDouble(fc1_(i));

    for (unsigned int i = 0; i < dim2; ++i) forcevec2(i) += FADUTILS::CastToDouble(fc2_(i));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::EvaluateStiffcContact(
    LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22, const double& pp,
    const TYPE& gap, const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm,
    const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
    const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, bool activecontact,
    bool linxi)
{
  // get dimensions for vector fc1 and fc2
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3;

  // temporary matrices for stiffness matrices
  LINALG::TMatrix<TYPE, dim1, dim1 + dim2> stiffc1;
  LINALG::TMatrix<TYPE, dim2, dim1 + dim2> stiffc2;

  // flag indicating assembly
  bool DoNotAssemble = false;

  // initialize stiffness to zero
  for (unsigned int i = 0; i < dim1; i++)
    for (unsigned int j = 0; j < (dim1 + dim2); j++) stiffc1(i, j) = 0.0;
  for (unsigned int i = 0; i < dim2; i++)
    for (unsigned int j = 0; j < (dim1 + dim2); j++) stiffc2(i, j) = 0.0;

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  // initialize delta_xi here because we need it outside the if-environment for FAD Check
  LINALG::TMatrix<TYPE, dim1 + dim2, 1> delta_xi;

  if (activecontact)
  {
    // auxiliary stiffmatrix for part III of linearization to avoid tensor notation
    LINALG::TMatrix<TYPE, dim1 + dim2, dim1 + dim2> stiffc_III;

    // initialize storage for linearizations
    LINALG::TMatrix<TYPE, 3, 1> distance;
    TYPE normdist = 0.0;
    LINALG::TMatrix<TYPE, dim1 + dim2, 1> delta_gap;
    LINALG::TMatrix<TYPE, 3, dim1 + dim2> delta_x1_minus_x2;
    LINALG::TMatrix<TYPE, 3, dim1 + dim2> delta_n;

    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************
    // linearization of contact point
    // (not needed in case of check for nodal points in contact)
    if (linxi)
      ComputeLinXi(delta_xi, x1, x2, dx1, ddx1, N1_i, N1_i_xi);
    else
    {
      delta_xi.Clear();
    }

    // evaluation of distance
    ComputeDistance(distance, normdist, normal, norm);

    // linearization of gap function which is equal to delta d
    ComputeLinGap(
        delta_gap, delta_xi, x1, x2, dx1, N1_i, normdist, normal, norm, gap, delta_x1_minus_x2);

    // linearization of normal vector
    ComputeLinNormal(delta_n, delta_xi, normal, norm, dx1, N1_i);

    //********************************************************************
    // evaluate contact stiffness
    // (1) stiffc1 of first element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************
    LINALG::TMatrix<TYPE, dim1, 1> N1T_normal;
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
      {
        N1T_normal(3 * j + i) += normal(i) * N1_i(j);
      }
    }

    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < (dim1 + dim2); j++)
      {
        stiffc1(i, j) = -pp * N1T_normal(i) * delta_gap(j);
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
      {
        for (unsigned int k = 0; k < dim1 + dim2; k++)
        {
          stiffc1(3 * j + i, k) -= pp * gap * N1_i(j) * delta_n(i, k);
        }
      }
    }

    //********************************************************************
    // part III
    //********************************************************************
    LINALG::TMatrix<TYPE, dim1, 1> N1xiT_normal;
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
      {
        N1xiT_normal(3 * j + i) += normal(i) * N1_i_xi(j);
      }
    }

    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        stiffc1(i, j) -= pp * gap * N1xiT_normal(i) * delta_xi(j);
      }
    }

    //********************************************************************
    // evaluate contact stiffness
    // (2) stiffc2 of second element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************

    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        stiffc2(i, j) = pp * normal(i) * delta_gap(j);
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int k = 0; k < dim1 + dim2; k++)
      {
        stiffc2(i, k) -= pp * gap * delta_n(i, k);
      }
    }
    //********************************************************************
    // part III
    //********************************************************************

    // There is no third part for spheres, since no shape functions are applied!
  }
  //**********************************************************************
  // no stiffness for inactive pairs
  //**********************************************************************
  else
  {
    // set flag to avoid assembly
    DoNotAssemble = true;
  }

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  if (not DoNotAssemble)
  {
#ifndef AUTOMATICDIFF
    // change sign of stiffc1 and stiffc2 due to time integration.
    // according to analytical derivation there is no minus sign, but for
    // our time integration methods the negative stiffness must be assembled.

    // Fixme check the minus sign !!!

    for (unsigned int j = 0; j < dim1; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat11(i, j) += -FADUTILS::CastToDouble(stiffc1(i, j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat21(i, j) += -FADUTILS::CastToDouble(stiffc2(i, j));
    }
    for (unsigned int j = 0; j < dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffmat12(i, j) += -FADUTILS::CastToDouble(stiffc1(i, dim1 + j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffmat22(i, j) += -FADUTILS::CastToDouble(stiffc2(i, dim1 + j));
    }
#else
    dserror("check implementation of AUTOMATICDIFF for BeamToSphereContactPair before using it!");

    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      for (unsigned int i = 0; i < dim1; i++)
        stiffc1_copy(i, j) =
            -FADUTILS::CastToDouble(fc1_(i).dx(j) + fc1_(i).dx(dim1 + dim2) * delta_xi(j));
      for (unsigned int i = 0; i < dim2; i++)
        stiffc2_copy(i, j) =
            -FADUTILS::CastToDouble(fc2_(i).dx(j) + fc2_(i).dx(dim1 + dim2) * delta_xi(j));
    }

#ifdef FADCHECKS
    dserror("check implementation of FADCHECKS for BeamToSphereContactPair before using it!");

    std::cout << "BTSPH Contact Pair: " << Element1()->Id() << " / " << Element2()->Id()
              << std::endl;

    std::cout << "stiffc1: " << std::endl;
    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        std::cout << std::setw(14) << -stiffc1(i, j).val();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "stiffc1_FAD: " << std::endl;
    for (unsigned int i = 0; i < dim1; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        std::cout << std::setw(14) << stiffc1_copy(i, j);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "stiffc2: " << std::endl;
    for (unsigned int i = 0; i < dim2; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        std::cout << std::setw(14) << -stiffc2(i, j).val();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "stiffc2_FAD: " << std::endl;
    for (unsigned int i = 0; i < dim2; i++)
    {
      for (unsigned int j = 0; j < dim1 + dim2; j++)
      {
        std::cout << std::setw(14) << stiffc2_copy(i, j);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

#endif  // FADCHECKS

#endif  // AUTOMATICDIFF
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeNormal(
    LINALG::TMatrix<TYPE, 3, 1>& normal, TYPE& gap, TYPE& norm,
    const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2)
{
  // compute non-unit normal
  for (unsigned int i = 0; i < 3; i++) normal(i) = x1(i) - x2(i);

  // compute length of normal
  norm = FADUTILS::VectorNorm<3>(normal);
  if (norm < NORMTOL) dserror("ERROR: Normal of length zero! --> change time step!");

  // compute unit normal and store it in class variable
  for (unsigned int i = 0; i < 3; i++)
  {
    normal(i) /= norm;
    normal_(i) = normal(i);
  }

  // evaluate scalar gap function
  ComputeGap(gap, norm);

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeGap(
    TYPE& gap, const TYPE& norm)
{
  // compute gap to be returned
  gap = norm - radius1_ - radius2_;
  gap_ = gap;

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeCoordsAndDerivs(
    LINALG::TMatrix<TYPE, 3, 1>& x1, LINALG::TMatrix<TYPE, 3, 1>& x2,
    LINALG::TMatrix<TYPE, 3, 1>& dx1, LINALG::TMatrix<TYPE, 3, 1>& ddx1,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi)
{
  // reset input variables
  x1.Clear();
  x2.Clear();
  dx1.Clear();
  ddx1.Clear();

#ifdef AUTOMATICDIFF
  // The 3*numnodes*numnodalvalues+3 primary DoFs are the components of the nodal positions /
  // tangents. The additional degree of freedom (+1) represents the dependency on the parameter
  // coordinates xi, which is necessary in beam contact.
  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; i++)
    ele1pos_(i).diff(i, 3 * numnodes * numnodalvalues + 3 + 1);

  for (unsigned int i = 0; i < 3; i++)
    ele2pos_(i).diff(3 * numnodes * numnodalvalues + i, 3 * numnodes * numnodalvalues + 3 + 1);
#endif  // AUTOMATICDIFF

  // compute output variables
  for (unsigned int i = 0; i < 3; i++) x2(i) = ele2pos_(i);

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
    {
      x1(i) += N1_i(j) * ele1pos_(3 * j + i);         // x1 = N1 * x~1
      dx1(i) += N1_i_xi(j) * ele1pos_(3 * j + i);     // dx1 = N1,xi * x~1
      ddx1(i) += N1_i_xixi(j) * ele1pos_(3 * j + i);  // ddx1 = N1,xixi * x~1
    }
  }

  // store coordinates of contact point into class variables
  for (unsigned int i = 0; i < 3; i++)
  {
    x1_(i) = x1(i);
    x2_(i) = x2(i);
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::GetShapeFunctions(
    LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
    LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
    LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, const TYPE& eta)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = BeamElement()->Shape();

  if (numnodalvalues == 1)
  {
    // get values and derivatives of shape functions
    DRT::UTILS::shape_function_1D(N1_i, eta, distype1);
    DRT::UTILS::shape_function_1D_deriv1(N1_i_xi, eta, distype1);
    DRT::UTILS::shape_function_1D_deriv2(N1_i_xixi, eta, distype1);
  }
  else if (numnodalvalues == 2)
  {
    /* TODO hard set distype to line2 in case of numnodalvalues_=2 because
     *  only 3rd order Hermite interpolation is used (always 2 nodes) */
    const DRT::Element::DiscretizationType distype1herm = DRT::Element::line2;

    // get values and derivatives of shape functions
    DRT::UTILS::shape_function_hermite_1D(N1_i, eta, beamele_reflength_, distype1herm);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi, eta, beamele_reflength_, distype1herm);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N1_i_xixi, eta, beamele_reflength_, distype1herm);
  }
  else
    dserror(
        "Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) "
        "values are valid!");

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeLinXi(
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
    const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
    const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3;

  // matrices to compute Lin_Xi (in case of beam-sphere contact: scalar and vector)
  TYPE L = 0.0;
  LINALG::TMatrix<TYPE, 1, dim1 + dim2> B;

  LINALG::TMatrix<TYPE, 3, 1> delta_x;

  // compute L
  for (unsigned int i = 0; i < 3; i++)
  {
    delta_x(i) = x1(i) - x2(i);
    L += dx1(i) * dx1(i) + delta_x(i) * ddx1(i);
  }

  if (L == 0) dserror("ERROR: L = 0");


  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < numnodes * numnodalvalues; ++j)
    {
      B(3 * j + i) += -delta_x(i) * N1_i_xi(j) - dx1(i) * N1_i(j);
    }
  }

  for (unsigned int i = 0; i < dim2; ++i) B(dim1 + i) += dx1(i);

  // finally the linearizations / directional derivatives
  for (unsigned int i = 0; i < dim1 + dim2; i++)
  {
    delta_xi(i) = B(i) / L;
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeDistance(
    LINALG::TMatrix<TYPE, 3, 1>& distance, TYPE& normdist,
    const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm)
{
  // compute distance vector
  for (unsigned int i = 0; i < 3; i++) distance(i) = normal(i) * norm;

  // compute scalar distance
  normdist = norm;

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeLinGap(
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_gap,
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
    const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
    const LINALG::TMatrix<TYPE, 3, 1>& dx1,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const TYPE& normdist,
    const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm, const TYPE& gap,
    LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_x1_minus_x2)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3;


  // delta g := delta_r/||delta_r||*auxiliary_matri1 delta d, with auxiliary_matri1 =
  // (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  LINALG::TMatrix<TYPE, 3, dim1 + dim2> auxiliary_matrix1(true);

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      auxiliary_matrix1(i, j) += dx1(i) * delta_xi(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
    {
      auxiliary_matrix1(i, 3 * j + i) += N1_i(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    auxiliary_matrix1(i, i + dim1) += -1.0;
  }

  // compute linearization of gap
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      delta_gap(j) += (x1(i) - x2(i)) * auxiliary_matrix1(i, j) / normdist;
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ComputeLinNormal(
    LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_normal,
    const LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
    const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm_delta_x,
    const LINALG::TMatrix<TYPE, 3, 1>& x1_xi,
    const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i)
{
  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3;

  // delta n := auxiliary_matri2*auxiliary_matrix1* delta d, with auxiliary_matri2 =
  // (I-nxn)/||r1-r2|| and auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  LINALG::TMatrix<TYPE, 3, dim1 + dim2> auxiliary_matrix1(true);
  LINALG::TMatrix<TYPE, 3, 3> auxiliary_matrix2(true);

  // compute auxiliary_matrix1
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < dim1 + dim2; j++)
    {
      auxiliary_matrix1(i, j) += x1_xi(i) * delta_xi(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < numnodes * numnodalvalues; j++)
    {
      auxiliary_matrix1(i, 3 * j + i) += N1_i(j);
    }
  }

  for (unsigned int i = 0; i < 3; i++)
  {
    auxiliary_matrix1(i, dim1 + i) += -1.0;
  }

  // compute auxiliary_matrix2
  for (unsigned int i = 0; i < 3; i++)
  {
    auxiliary_matrix2(i, i) += 1.0 / norm_delta_x;
    for (unsigned int j = 0; j < 3; j++)
    {
      auxiliary_matrix2(i, j) += -normal(i) * normal(j) / norm_delta_x;
    }
  }

  // compute linearization of normal vector
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < dim1 + dim2; k++)
        delta_normal(i, k) += auxiliary_matrix2(i, j) * auxiliary_matrix1(j, k);

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::ResetState(
    const std::vector<double>& centerline_dofvec_ele1,
    const std::vector<double>& centerline_dofvec_ele2)
{
  if (centerline_dofvec_ele1.size() != 3 * numnodes * numnodalvalues)
    dserror("size mismatch! expected %d values for centerline_dofvec_ele1, but got %d",
        3 * numnodes * numnodalvalues, centerline_dofvec_ele1.size());

  if (centerline_dofvec_ele2.size() != 3)
    dserror("size mismatch! expected %d values for centerline_dofvec_ele2, but got %d", 3,
        centerline_dofvec_ele1.size());

  for (unsigned int i = 0; i < 3 * numnodalvalues * numnodalvalues; ++i)
    ele1pos_(i) = centerline_dofvec_ele1[i];

  for (unsigned int i = 0; i < 3; ++i) ele2pos_(i) = centerline_dofvec_ele2[i];
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes, numnodalvalues>::Print(
    std::ostream& out) const
{
  // ToDo add further information here
  out << "\nInstance of BeamToSphereContactPair: element GIDs " << Element1()->Id()
      << " (beam) and " << Element2()->Id() << " (sphere)";
  out << "\ncontactflag_: " << contactflag_ << "\tnodalcontactflag: " << nodalcontactflag_[0]
      << nodalcontactflag_[1];
  out << "\ngap_: " << gap_;
  out << "\nxicontact_: " << xicontact_;
  out << "\nnormal_: " << normal_;
  out << "\nx1_: " << x1_;
  out << "\nx2_: " << x2_;
  out << "\nele1pos_: " << ele1pos_;
  out << "\nele2pos_: " << ele2pos_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToSphereContactPair<numnodes,
    numnodalvalues>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  if (GetContactFlag())
  {
    out << "    " << std::setw(5) << std::left << Element1()->Id() << "(" << std::setw(3)
        << std::right << -1 << "/" << std::setw(3) << -1 << ")"
        << " " << std::setw(5) << std::left << Element2()->Id() << "(" << std::setw(3) << std::right
        << -1 << "/" << std::setw(3) << -1 << ")"
        << "  BTSPH ";

    out << std::setw(9) << std::left << std::setprecision(2) << xicontact_ << std::setw(9)
        << std::left << std::setprecision(2) << -1 << std::setw(9) << std::left
        << std::setprecision(3) << -1 << std::setw(12) << std::left << std::scientific << gap_
        << std::setw(12) << std::left << std::scientific
        << FADUTILS::CastToDouble<TYPE, 3, 1>(fc2_).Norm2() << std::setprecision(6)
        << std::resetiosflags(std::ios::scientific) << std::right;

    out << "\n";
  }
}

// explicit template instantiations
template class BEAMINTERACTION::BeamToSphereContactPair<2, 1>;
template class BEAMINTERACTION::BeamToSphereContactPair<3, 1>;
template class BEAMINTERACTION::BeamToSphereContactPair<4, 1>;
template class BEAMINTERACTION::BeamToSphereContactPair<5, 1>;
template class BEAMINTERACTION::BeamToSphereContactPair<2, 2>;
