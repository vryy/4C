/*----------------------------------------------------------------------------*/
/*!
\file beam3tosolidcontact.cpp

\brief One beam and solid contact pair (two elements)

\level 3

\maintainer Alexander Popp
*/
/*----------------------------------------------------------------------------*/

#include "beam3tosolidcontact.H"
#include "../drt_beaminteraction/beam3contact_utils.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_beaminteraction/beam3contact_tangentsmoothing.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

/*----------------------------------------------------------------------*
 | Constructor (public)                                     meier 01/14 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Beam3tosolidcontact(
    const DRT::Discretization& pdiscret,
    const DRT::Discretization& cdiscret,
    const std::map<int,int>& dofoffsetmap,
    DRT::Element* element1,
    DRT::Element* element2,
    Teuchos::ParameterList beamcontactparams):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
sgn_(1.0),
firstcall_(true),
gap_(0.0),
gap_original_(0.0),
contactflag_(false),
elementscolinear_(false),
elementscrossing_(false),
shiftnodalvalues_(false),
xi1_(0.0),
xi2_(0.0),
gmshDebugPoints_()
{
  for (int i = 0; i < 3*numnodes*numnodalvalues; i++)
    ele1pos_(i) = 0.0;

  for (int i = 0; i < 3*numnodessol; i++)
    ele2pos_(i) = 0.0;

  normalsets_.clear();
  normalsets_.resize(0);

  normalsets_old_.clear();
  normalsets_old_.resize(0);

  return;
}
/*----------------------------------------------------------------------*
 | End: constructor                                                     |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | copy-constructor (public)                                            |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Beam3tosolidcontact(const Beam3tosolidcontact& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffsetmap_(old.dofoffsetmap_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}
/*----------------------------------------------------------------------*
 | End: copy-constructor                                                |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate the element (public)                                        |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
bool CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Evaluate(
     LINALG::SparseMatrix& stiffmatrix,
     Epetra_Vector& fint,
     const double& pp)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

#if defined(GMSHDEBUG) || defined(SAVEFORCE) || defined(GMSHFORCE)
  // Clear vector for Gmsh debug
  gmshDebugPoints_.clear();
  gmshDebugPoints_.resize(0);
#endif


  // Find contact interval borders
  // -----------------------------------------------------------------

  // Vector for contact interval border parameter sets (xi1, xi2, eta and index of fixed paramater)
  std::vector<std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int> > parsets;

  // Find contact interval borders
  GetContactIntervalBorders(parsets);

  // -----------------------------------------------------------------
  // End: Find contact interval borders


  // Compute contact forces and stiffness for all contact intervals
  // -----------------------------------------------------------------

  // Temporary vectors for contact forces
  LINALG::TMatrix<TYPEBTS, dim1, 1> fc1(true);
  LINALG::TMatrix<TYPEBTS, dim2, 1> fc2(true);

  // Temporary matrices for contact stiffness
  LINALG::TMatrix<TYPEBTS, dim1, dim1+dim2> stiffc1(true);
  LINALG::TMatrix<TYPEBTS, dim2, dim1+dim2> stiffc2(true);

  // Temporary matrices for contact stiffness calculated via FAD
  // TODO: Declare stiffc1_FAD and stiffc2_FAD only if needed (#ifdef AUTOMATICDIFFBTS)
  LINALG::TMatrix<TYPEBTS, dim1, dim1+dim2> stiffc1_FAD(true);
  LINALG::TMatrix<TYPEBTS, dim2, dim1+dim2> stiffc2_FAD(true);

  // Set total gap for beam to solid contact pair to zero
  // gap_ is used for calculating the contact flag of the beam to solid contact contact pair
  gap_ = 0.0;

  // Clear vector containing pairs of surface unit normal vector nD and beam parameter eta of current time step
  normalsets_.clear();
  normalsets_.resize(0);

  // Initialize flag indicating if the beam to solid contact pair has contact
  // For active contact assemble contact forces and stiffness
  bool doAssemble = false;

  // Loop over all contact intervals, calculate and sum up contact forces (fc1, fc2) and stiffness (stiffc1, stiffc2)
  for (int i = 0; i < (int)parsets.size() - 1; i += 2)
  {
    // Initialize flag indicating if the current contact interval has contact
    bool doAssembleContactInterval = false;

    // Two parameter sets indicate the beginning (a) and the end (b) of the contact interval
    // TODO: Use stiffc1_FAD and stiffc2_FAD only if needed (#ifdef AUTOMATICDIFFBTS)
    EvaluateContactInterval(pp, parsets[i], parsets[i + 1], fc1, fc2, stiffc1, stiffc2,
        doAssembleContactInterval, stiffc1_FAD, stiffc2_FAD);

    // If at least one contact interval has contact, set the doAssemble flag true
    if (!doAssemble && doAssembleContactInterval)
    {
      doAssemble = true;
    }
  }

  // -----------------------------------------------------------------
  // End: Compute contact forces and stiffness for all contact intervals


  // Assemble contact forces and contact stiffness
  // -----------------------------------------------------------------

  // Evaluate contact status for beam to solid contact pair
  contactflag_ = (gap_ < 0.0);

#ifdef FDCHECKSTIFFNESS
  // Calculate contact stiffness with finite difference
  if (doAssemble)
  {
    FDCheckStiffness(pp, fc1, fc2, stiffc1, stiffc2);
  }
#endif

#ifdef FADCHECKSTIFFNESS
  // Compare calculated contact stiffness with FAD
  if (doAssemble)
  {
    FADCheckStiffness(stiffc1, stiffc2, stiffc1_FAD, stiffc2_FAD);
  }
#endif

  // Finally assemble fc1, fc2 and stiffc1, stiffc2 for beam to solid contact pair
  if (doAssemble)
  {
    AssembleFcAndStiffcContact(fc1, fc2, &fint, stiffc1, stiffc2, stiffmatrix);
  }

  // -----------------------------------------------------------------
  // end: Assemble contact forces and contact stiffness

  return true;
}
/*----------------------------------------------------------------------*
 | End: Evaluate the element                                            |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Evaluate contact interval                                            |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateContactInterval(
    const double& pp,
    const std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int>& parset_a,
    const std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int>& parset_b,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 1>& fc1,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 1>& fc2,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2,
    bool& doAssembleContactInterval,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1_FAD,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2_FAD)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // NOTE: Linearization lin() = (),d * lin_d with directional derivative (),d

  // Beam and surface element parameters
  TYPEBTS xi1 = 0.0;
  TYPEBTS xi2 = 0.0;
  TYPEBTS eta = 0.0;

  // Vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2(true);                     // = N2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1(true);                 // = N2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2(true);                 // = N2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1xi1(true);              // = N2,xi1xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2xi2(true);              // = N2,xi2xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1xi2(true);              // = N2,xi1xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2xi1(true);              // = N2,xi2xi1

  // Coords and derivatives of beam and surface element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 1> x2(true);                                 // = x2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1(true);                             // = x2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2(true);                             // = x2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1xi1(true);                          // = x2,xi1xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2xi2(true);                          // = x2,xi2xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1xi2(true);                          // = x2,xi1xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2xi1(true);                          // = x2,xi2xi1

  // Distance vector its norm and distance unit vector
  LINALG::TMatrix<TYPEBTS, 3, 1> rD(true);                                 // = r1 - x2
  TYPEBTS norm_rD = 0.0;                                                   // = ||rD||
  LINALG::TMatrix<TYPEBTS, 3, 1> nD(true);                                 // = rD / norm_rD

  // Surface tangent cross product, its norm and unit surface normal vector
  LINALG::TMatrix<TYPEBTS, 3, 1> a2(true);                                 // = x2_xi1 x x2_xi2
  TYPEBTS norm_a2 = 0.0;                                                   // = ||a||
  LINALG::TMatrix<TYPEBTS, 3, 1> n2(true);                                 // = a / norm_a

  // Sign variable and gap function
  double sgn = 1.0;                                                        // = sign(nD * n2)
  TYPEBTS gap = 0.0;                                                       // = n2 * rD - radius1


  // Calculate radius of beam element via moment of inertia (only valid for beams with circular cross section)
  const DRT::ELEMENTS::Beam3Base* beamele1 = static_cast<const DRT::ELEMENTS::Beam3Base*>(element1_);

  if (beamele1 == NULL)
    dserror("cast to beam base failed!");

  const double radius1 = MANIPULATERADIUS * beamele1->GetCircularCrossSectionRadiusForInteractions();

  // Get contact interval borders eta_a and eta_b
  TYPEBTS eta_a = parset_a.first(2);
  TYPEBTS eta_b = parset_b.first(2);


  // Calculate linearization of contact interval borders
  // TODO: Calculate only if the current contact interval has contact
  // -----------------------------------------------------------------

  // Initialize directional derivatives of contact interval borders
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_a_d(true);
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_b_d(true);

  // Compute directional derivatives eta_a_d (i = 0) and eta_b_d (i = 1)
  // TODO: May this can be done in a more beautiful way
  for (int i = 0; i < 2; i++)
  {
    LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_d(true);
    std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int> parset;
    switch (i)
    {
      case 0:
        parset = parset_a;
        break;
      case 1:
        parset = parset_b;
        break;
    }

    // Get parameters from contact interval border parameter set
    TYPEBTS xi1 = parset.first(0);
    TYPEBTS xi2 = parset.first(1);
    TYPEBTS eta = parset.first(2);
    const int fixed_par = parset.second;

    // Initialize directional derivatives of surface parameters xi1 and xi2. They are used temporary,
    // because we are only interested in the directional derivative of the beam parameter eta
    LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi1_d(true);
    LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi2_d(true);

#ifdef FADCHECKLINCONTACTINTERVALBORDER
    // Set known parameters xi1, xi2, eta and element positions as primary variables
    // for checking the linearization of the contact interval borders with FAD
    BEAMCONTACT::SetFADParCoordDofs<numnodessol, numnodes, numnodalvalues>(xi1, xi2, eta);
    BEAMCONTACT::SetFADDispDofs<numnodessol, numnodes, numnodalvalues>(ele1pos_, ele2pos_, 3);
#endif

    // Update shape functions and their derivatives for beam and surface element
    GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
    GetSurfShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1, xi1, xi2);

    // Update coordinates and derivatives for beam and surface element
    ComputeBeamCoordsAndDerivs(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);
    ComputeSurfCoordsAndDerivs(x2, x2_xi1, x2_xi2, x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1,
        N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1);

    // Compute unit distance vector nD
    ComputeDistanceNormal(r1, x2, rD, norm_rD, nD);

    // Add nD and eta to normalsets_ of the current time step
    std::pair<TYPEBTS, LINALG::TMatrix<TYPEBTS, 3, 1> > normalset;
    normalset.first = eta;
    normalset.second = nD;
    normalsets_.push_back(normalset);

#if defined(GMSHDEBUG) || defined(GMSHFORCE) || defined(SAVEFORCE)
    // Compute some variables at contact interval border
    ComputeSurfaceNormal(x2_xi1, x2_xi2, a2, norm_a2, n2);
    gap = FADUTILS::ScalarProduct(n2, rD) - radius1;
    double fp = 0.0;
    if (gap < 0)
    {
      TYPEBTS dfp = 0.0;
      EvaluatePenaltyForceLaw(pp, gap, fp, dfp);
    }

    // Store variables at contact interval border in gmshDebugPoints_
    gmshDebugPoint gmshDebugPoint;
    gmshDebugPoint.r1 = r1;
    gmshDebugPoint.x2 = x2;
    gmshDebugPoint.n2 = n2;
    gmshDebugPoint.gap = gap;
    gmshDebugPoint.fp = fp;
    gmshDebugPoint.type = -1 + 2 * i; // -1: Left contact interval border, +1: Right contact interval border
    gmshDebugPoints_.push_back(gmshDebugPoint);
#endif

    // If eta is fixed then its linearization is zero, otherwise compute directional derivative
    if (fixed_par == 2)
    {
      // Set directional derivative of contact interval border to zero
      eta_d.Clear();
    }
    else
    {
      xi1_d.Clear();
      xi2_d.Clear();
      eta_d.Clear();

      // Compute directional derivatives of not fixed parameters at contact interval border,
      // only the directional derivative eta_d is needed
      ComputeLinParameter(fixed_par, xi1_d, xi2_d, eta_d, rD, r1_eta, x2_xi1, x2_xi2,
          x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1, N1, N2, N2_xi1, N2_xi2);

#ifdef FADCHECKLINCONTACTINTERVALBORDER
      // Print some information
      std::cout << std::endl << "FAD-Check: Linearization of contact interval borders at "
          << "xi1: " << xi1.val() << ", xi2: " << xi2.val() << ", eta: " << eta.val() << " with fixed_par: " << fixed_par;
      switch (i)
      {
        case 0:
          std::cout << " (eta_a_d): " << std::endl;
          break;
        case 1:
          std::cout << " (eta_b_d): " << std::endl;
          break;
      }

      // Directional derivatives of eta and temporary derivatives of xi1 and xi2
      LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi1_d_FAD(true);
      LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi2_d_FAD(true);
      LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_d_FAD(true);

      // Compute linearization of contact interval borders eta_a (if i == 0) or eta_b (if i == 1) with FAD
      FADCheckLinParameter(fixed_par, rD, x2_xi1, x2_xi2, xi1_d_FAD, xi2_d_FAD, eta_d_FAD, xi1_d, xi2_d, eta_d);
#endif
    }

    switch (i)
    {
      case 0:
        eta_a_d = eta_d;
        break;
      case 1:
        eta_b_d = eta_d;
        break;
    }
  }
  // -----------------------------------------------------------------
  // End: Calculate linearization of contact interval borders


  // Evaluate contact force and stiffness at each Gauss point
  // -----------------------------------------------------------------

  // Initialize directional derivative of eta at Gauss point
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_d(true);

  // Get Gauss points and weights in interval [-1, 1]
  DRT::UTILS::IntegrationPoints1D gaussPoints = DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GAUSSRULE);

  // Loop over all Gauss points in contact interval [eta_a, eta_b]
  for (int i_gp = 0; i_gp < gaussPoints.nquad; i_gp++)
  {
    // Initialize flag indicating if Gauss point has contact
    bool contactflag = false;

    // Initialize flag indicating if projection at Gauss point is allowed
    bool proj_allowed = false;

    // Set surface parameters and gap to zero
    xi1 = 0.0;
    xi2 = 0.0;
    gap = 0.0;

    // Get integration point and weight
    const double x_gp = gaussPoints.qxg[i_gp][0];
    const double w_gp = gaussPoints.qwgt[i_gp];

    // Calculate beam parameter eta at Gauss point using a linear transformation
    eta = 0.5 * (eta_b - eta_a) * x_gp + 0.5 * (eta_a + eta_b);

    // Calculate directional derivative of eta at Gauss point
    eta_d.Clear();
    for (int i = 0; i < dim1+dim2; i++)
    {
      eta_d(i) = 0.5 * (1 - x_gp) * eta_a_d(i) + 0.5 * (1 + x_gp) * eta_b_d(i);
    }

    // Projection of Gauss point on surface with fixed eta to get xi1 and xi2
    Projection(2, xi1, xi2, eta, proj_allowed);

    // Check if the projected point found for this contact pair is really on the
    // considered surface element, i.e. xi1 and xi2 in [-1, 1]
    if (fabs(xi1) > 1.0 + XIETATOL || fabs(xi2) > 1.0 + XIETATOL)
    {
      dserror("Projection on surface is not on the considered surface element.");

      // std::cout << "Projection of Gauss point " << i_gp + 1 << " of " << gaussPoints.nquad
      //     << " (at eta: " << eta << ") on surface is not on the considered surface element." << std::endl;
      // std::cout << "xi1: " << xi1 << std::endl;
      // std::cout << "xi2: " << xi2 << std::endl;
      // continue;
    }

#ifdef FADCHECKSTIFFNESS
    // Set known parameters xi1, xi2, contact interval borders eta_a, eta_b and element
    // positions as primary variables for checking the contact stiffness with FAD
    BEAMCONTACT::SetFADParCoordDofs<numnodessol, numnodes, numnodalvalues>(xi1, xi2, eta_a, eta_b);
    BEAMCONTACT::SetFADDispDofs<numnodessol, numnodes, numnodalvalues>(ele1pos_, ele2pos_, 4);

    // Calculate eta again depending on eta_a and eta_b, because at this point all FAD paramater coordinates are known
    eta = 0.5 * (eta_b - eta_a) * x_gp + 0.5 * (eta_a + eta_b);
#endif

#ifdef FADCHECKLINGAUSSPOINT
    // Set known parameters xi1, xi2, eta and element positions as primary variables for
    // checking the linearization of surface parameters xi1 and xi2, gap and surface unit normal vector n2
    BEAMCONTACT::SetFADParCoordDofs<numnodessol, numnodes, numnodalvalues>(xi1, xi2, eta);
    BEAMCONTACT::SetFADDispDofs<numnodessol, numnodes, numnodalvalues>(ele1pos_, ele2pos_, 3);
#endif

    // Update shape functions and their derivatives for beam and surface element
    GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
    GetSurfShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1, xi1, xi2);

    // Update coordinates and derivatives for beam and surface element
    ComputeBeamCoordsAndDerivs(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);
    ComputeSurfCoordsAndDerivs(x2, x2_xi1, x2_xi2, x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1,
        N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1);

    // Compute distance vector rD, its norm norm_rD and unit distance vector nD
    ComputeDistanceNormal(r1, x2, rD, norm_rD, nD);

    // Compute tangent cross product a2, its norm norm_a2 and surface unit normal vector n2
    ComputeSurfaceNormal(x2_xi1, x2_xi2, a2, norm_a2, n2);

    // Compute sign of projection of distance vector rD on surface unit normal vector n2
    sgn = FADUTILS::CastToDouble(FADUTILS::Signum(FADUTILS::ScalarProduct(nD, n2)));

    // Compute gap at Gauss point
    gap = FADUTILS::ScalarProduct(n2, rD) - radius1; // = sgn * norm_rD - radius1;

    // If the beam element center line lies under the surface element and there is no change of sgn from
    // the last to the current time step, the contact is not active (here called projection not allowed)
    // and contactflag is false (initialization). Otherwise check the contact status at this Gauss point.
    if (proj_allowed)
    {
      CheckContactStatus(pp, gap, contactflag);
    }

#if defined(GMSHDEBUG) || defined(GMSHFORCE) || defined(SAVEFORCE)
    if (sgn > 0 || contactflag)
    {
      // Insert variables at current Gauss point in gmshDebugPoints_ before variables of contact interval end
      double fp = 0.0;
      if (gap < 0)
      {
        TYPEBTS dfp = 0.0;
        EvaluatePenaltyForceLaw(pp, gap, fp, dfp);
      }
      gmshDebugPoint gmshDebugPoint;
      gmshDebugPoint.r1 = r1;
      gmshDebugPoint.x2 = x2;
      gmshDebugPoint.n2 = n2;
      gmshDebugPoint.gap = gap;
      gmshDebugPoint.fp = fp;
      gmshDebugPoint.type = 0; // 0: Gauss point
      gmshDebugPoints_.insert(gmshDebugPoints_.end() - 1, gmshDebugPoint);
    }
#endif

    // Evaluate and sum up contact forces and stiffness only if the current Gauss point has contact
    if (contactflag)
    {
      // Sum up all negative gaps for getting a total gap for the current beam to solid contact pair
      gap_ += gap;

      // Get Jacobi factor of beam element
      const double jacobi = (static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();

      // Evaluate penalty force law with penalty parameter pp and gap to get force fp and derivative dfp = dfp / dgap
      TYPEBTS fp = 0.0;
      TYPEBTS dfp = 0.0;
      EvaluatePenaltyForceLaw(pp, gap, fp, dfp);

      // Evaluate contact forces at current Gauss point
      EvaluateFcContact(fp, fc1, fc2, eta_a, eta_b, w_gp, sgn, nD, n2, N1, N2, jacobi);

      // Evaluate contact stiffness at current Gauss point
      EvaluateStiffcContact(fp, dfp, gap, stiffc1, stiffc2, sgn, eta_a, eta_b, eta_d, eta_a_d, eta_b_d, w_gp,
          rD, norm_rD, nD, a2, norm_a2, n2, r1_eta, x2_xi1, x2_xi2, x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1,
          N1, N1_eta, N2, N2_xi1, N2_xi2, jacobi, stiffc1_FAD, stiffc2_FAD);

      // Set the doAssembleContactInterval flag true, if at least one Gauss point
      // inside the current contact interval has active contact
      if (!doAssembleContactInterval)
      {
        doAssembleContactInterval = true;
      }
    }

    // Save current parameter eta and unit distance vector nD in normalsets_ of the current time step
    // by inserting this pair before the last normalset in normalsets_ (contact interval end)
    std::pair<TYPEBTS, LINALG::TMatrix<TYPEBTS, 3, 1> > normalset;
    normalset.first = eta;
    normalset.second = nD;
    normalsets_.insert(normalsets_.end() - 1, normalset);
  }
  // -----------------------------------------------------------------
  // End: Evaluate contact force and stiffness at each Gauss point

  return;
}
/*----------------------------------------------------------------------*
 | End: Evaluate contact interval                                       |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Evaluate penalty force law                                           |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluatePenaltyForceLaw(
    const double& pp,
    const TYPEBTS& gap,
    TYPEBTS& fp,
    TYPEBTS& dfp)
{
#ifdef LINPENALTY
  // Linear penalty force law
  fp = - pp * gap;
  dfp = -pp;
#endif

#ifdef QUADPENALTY
  // Quadratic penalty force law
  fp = pp * gap * gap;
  dfp = 2 * pp * gap;
#endif

#ifdef ARBITPENALTY
  const double g0 = G0;

  // Quadratic penalty force law
  if (ARBITPENALTY == 1) // Cubic regularization for positive gaps
  {
    double k = 2.0;
    double f0 = pp * g0 / k;
    double factor_a = -pp / (g0 * g0) + 2 * f0 / (g0 * g0 * g0);
    double factor_b = 2 * pp / g0 - 3 * f0 / (g0 * g0);
    double factor_c = -pp;
    double factor_d = f0;
    if (gap > 0.0)
    {
      // std::cout << "Regularized Penalty!" << std::endl;
      fp = factor_a * gap * gap * gap + factor_b * gap * gap + factor_c * gap + factor_d;
      dfp = 3 * factor_a * gap * gap + 2 * factor_b * gap + factor_c;
    }
    else
    {
      // std::cout << "Linear Penalty!" << std::endl;
      fp = f0 - pp * gap;
      dfp= -pp;
    }
  }
  else if (ARBITPENALTY == 2) // Quadratic regularization for negative gaps
  {
    if (gap > -g0)
    {
      // std::cout << "Regularized Penalty!" << std::endl;
      fp = pp / (2.0 * g0) * gap * gap;
      dfp = pp / g0 * gap;
    }
    else
    {
      // std::cout << "Linear Penalty!" << std::endl;
      fp = -pp * (gap + g0 / 2.0);
      dfp = -pp;
    }
  }
  else if (ARBITPENALTY == 3) // Quadratic regularization for positiv gaps
  {
    double f0 = g0 * pp / 2.0;
    double factor_a = pp / g0 - f0 / (g0 * g0);
    double factor_b = -pp;
    double factor_c = f0;

    if (gap > 0)
    {
      // std::cout << "Regularized Penalty!" << std::endl;
      fp = factor_a * gap * gap + factor_b * gap + factor_c;
      dfp = 2 * factor_a * gap + factor_b;
    }
    else
    {
      // std::cout << "Linear Penalty!" << std::endl;
      fp = f0 - pp * gap;
      dfp = -pp;
    }
  }
  else if (ARBITPENALTY == 4) // Double quadratic regularization for positiv gaps
  {

    double f0 = 0.25;
    double g1 = 1.8 * f0 / pp;
    double c_tilde = f0;
    double b_tilde = -pp;
    double a_bar = (2 * f0 - pp * g1) / (2 * g0 * (g0 - g1));
    double b_bar = -2 * g0 * a_bar;
    double c_bar= -g0 * g0 * a_bar -g0 * b_bar;
    double a_tilde = (2 * g1 * a_bar + b_bar - b_tilde) / (2 * g1);

    if (gap > g1)
    {
      // std::cout << "Regularized Penalty: g1 < gap < g0!" << std::endl;
      fp = a_bar * gap * gap + b_bar * gap + c_bar;
      dfp = 2 * a_bar * gap + b_bar;
    }
    else if (gap > 0)
    {
      // std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
      fp = a_tilde * gap * gap + b_tilde * gap + c_tilde;
      dfp = 2 * a_tilde * gap + b_tilde;
    }
    else
    {
      // std::cout << "Linear Penalty!" << std::endl;
      fp = f0 - pp * gap;
      dfp = -pp;
    }
  }
  else if(ARBITPENALTY == 5) // Exponential regularization for positiv gaps
  {
    double f0 = 0.25;

    if (gap > 0)
    {
      // std::cout << "Regularized Penalty: 0 < gap < g1!" << std::endl;
      fp = f0 * exp(-pp * gap / f0);
      dfp = -pp * exp(-pp * gap / f0);
    }
    else
    {
      // std::cout << "Linear Penalty!" << std::endl;
      fp = f0 - pp * gap;
      dfp = -pp;
    }
  }
  else
    dserror("Only the values 1,2,3,4 are possible for ARBITPENALTY!");
#endif
}
/*----------------------------------------------------------------------*
 | End: Evaluate penalty force law                                      |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute contact forces                                               |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateFcContact(
    const TYPEBTS& fp,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 1>& fc1,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 1>& fc2,
    const TYPEBTS& eta_a,
    const TYPEBTS& eta_b,
    const double& w_gp,
    const double& sgn,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& nD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& n2,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2,
    const double& jacobi)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // NOTE: Macauley-bracket <-gap> = -gap if gap < 0


  // Evaluate contact forces acting on beam element
  // -----------------------------------------------------------------

  // Compute fc1 at Gauss point for beam element
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < 3; j++)
      fc1(i) += jacobi * w_gp * fp * N1(j,i) * n2(j) * 0.5 * (eta_b - eta_a);

  // -----------------------------------------------------------------
  // end: Evaluate contact forces acting on beam element


  // Evaluate contact forces acting on solid element
  // -----------------------------------------------------------------

  // Compute fc2 at Gauss point for surface element
  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < 3; j++)
      fc2(i) += -jacobi * w_gp * fp * N2(j,i) * n2(j) * 0.5 * (eta_b - eta_a);

  // -----------------------------------------------------------------
  // end: Evaluate contact forces acting on solid element

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute contact forces                                          |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate contact stiffness                                           |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateStiffcContact(
    const TYPEBTS& fp,
    const TYPEBTS& dfp,
    const TYPEBTS& gap,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2,
    const double& sgn,
    const TYPEBTS& eta_a,
    const TYPEBTS& eta_b,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_a_d,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_b_d,
    const double& w_gp,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const TYPEBTS& norm_rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& nD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& a2,
    const TYPEBTS& norm_a2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& n2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& r1_eta,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1_eta,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi2,
    const double& jacobi,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1_FAD,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2_FAD)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // NOTE: Linearization lin() = (),d * lin_d with directional derivative (),d

  // Initialize storage for directional derivatives
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi1_d(true);
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi2_d(true);
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> gap_d(true);
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> rD_d(true);
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> nD_d(true);
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> n2_d(true);

  // Compute directional derivative of parameters xi1 and xi2 with given eta_d (par_fixed = 2)
  ComputeLinParameter(2, xi1_d, xi2_d, eta_d, rD, r1_eta, x2_xi1, x2_xi2,
      x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1, N1, N2, N2_xi1, N2_xi2);

  // Compute directional derivatives of gap and distance vector rD
  ComputeLinGap(gap_d, xi1_d, xi2_d, eta_d, sgn, rD, nD, n2, norm_rD, r1_eta, x2_xi1, x2_xi2, N1, N2, rD_d);

  // Compute directional derivative of normal vector n2
  ComputeLinNormal(nD_d, nD, norm_rD, n2_d, n2, norm_a2, rD_d, xi1_d, xi2_d, x2_xi1, x2_xi2,
      x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1, N2_xi1, N2_xi2);


#ifdef FADCHECKLINGAUSSPOINT
  // Print information for linearization of surface parameters at Gauss point
  std::cout << std::endl << "FAD-Check: Linearization of surface parameters at Gauss point" << std::endl;

  // Initialize directional derivatives of surface parameters xi1 and xi2 for FAD
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi1_d_FAD(true);
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> xi2_d_FAD(true);

  // At this point the directional derivative of eta is already known and is used for calculating xi1_d and xi2_d
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> eta_d_FAD = eta_d;

  // Compute directional derivatives of surface parameters xi1 and xi2 for fixed eta (par_fixed = 2) with FAD
  FADCheckLinParameter(2, rD, x2_xi1, x2_xi2, xi1_d_FAD, xi2_d_FAD, eta_d_FAD, xi1_d, xi2_d, eta_d);

  // Print information for linearization of gap and distance vector
  std::cout << std::endl << "FAD-Check: Linearization of gap and distance vector" << std::endl;

  // Initialize directional derivative of gap and distance vector rD for FAD
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> gap_d_FAD(true);
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> rD_d_FAD(true);

  // Compute directional derivative of gap and distance vector rD with FAD
  FADCheckLinGapAndDistanceVector(gap, rD, xi1_d, xi2_d, eta_d, gap_d_FAD, rD_d_FAD, gap_d, rD_d);

  // Print information for linearization of normal vector
  std::cout << std::endl << "FAD-Check: Linearization of normal vector" << std::endl;

  // Initialize directional derivatives of unit distance vector and normal vector n2 for FAD
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> nD_d_FAD(true);
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> n2_d_FAD(true);

  // Compute directional derivatives of unit distance vector and normal vector n2 with FAD
  FADCheckLinNormal(nD, n2, xi1_d, xi2_d, eta_d, nD_d_FAD, n2_d_FAD, nD_d, n2_d);
#endif


  // Initialize flag indicating if the complete stiffness is used
  bool completestiff = true;
#ifdef GAP0
  // In order to accelerate convergence, we only apply the basic stiffness part in case of very large negative gaps
  // NOTE: GAP0 > 0
  if (gap < -GAP0)
  {
    completestiff = false;
    std::cout << "gap: " << gap << " < -GAP0: " << -GAP0 << " => completestiff: " << completestiff << std::endl;
  }
#endif


  // NOTE: Macaulay bracket <-gap> = -gap and Heaviside function H(-gap) = 1 if gap < 0,
  // lin_<-gap> = -H(-gap) * lin_gap = -lin_gap <=> <-gap>_d = -gap_d


  // Evaluate contact stiffness of beam element
  // -----------------------------------------------------------------

  // Part I: gap_d
  LINALG::TMatrix<TYPEBTS, dim1, 1> N1T_n2(true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1; j++)
      N1T_n2(j) += N1(i,j) * n2(i);

  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim1+dim2; j++)
      stiffc1(i,j) += jacobi * w_gp * dfp * N1T_n2(i) * gap_d(j) * 0.5 * (eta_b - eta_a);

  if (completestiff)
  {
    // Part II: n2_d
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < dim1;j++)
        for (int k = 0; k < dim1+dim2; k++)
          stiffc1(j,k) += jacobi * w_gp * fp * N1(i,j) * n2_d(i,k) * 0.5 * (eta_b - eta_a);

    // Part III: N1_d
    LINALG::TMatrix<TYPEBTS, dim1, 1> N1_etaT_n2(true);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < dim1; j++)
        N1_etaT_n2(j) += N1_eta(i,j) * n2(i);
    for (int i = 0; i < dim1; i++)
      for (int j = 0; j < dim1+dim2; j++)
        stiffc1(i,j) += jacobi * w_gp * fp * N1_etaT_n2(i) * eta_d(j) * 0.5 * (eta_b - eta_a);

    // Part IV: eta_a_d and eta_b_d (Jacobi factor contact interval)
    for (int i = 0; i < dim1; i++)
      for (int j = 0; j < dim1+dim2; j++)
        stiffc1(i,j) += jacobi * w_gp * fp * N1T_n2(i) * 0.5 * (eta_b_d(j) - eta_a_d(j));
  }

  // -----------------------------------------------------------------
  // end: Evaluate contact stiffness of beam element


  // Evaluate contact stiffness of solid element
  // -----------------------------------------------------------------

  // Part I: gap_d
  LINALG::TMatrix<TYPEBTS, dim2, 1> N2T_n2(true);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim2; j++)
      N2T_n2(j) += N2(i,j) * n2(i);

  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < dim1+dim2; j++)
      stiffc2(i,j) += -jacobi * w_gp * dfp * N2T_n2(i) * gap_d(j) * 0.5 * (eta_b - eta_a);

  if (completestiff)
  {
    // Part II: n2_d
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < dim2;j++)
        for (int k = 0; k < dim1+dim2; k++)
          stiffc2(j,k) += -jacobi * w_gp * fp * N2(i,j) * n2_d(i,k) * 0.5 * (eta_b - eta_a);

    // Part III: N2_d
    LINALG::TMatrix<TYPEBTS, dim2, 1> N2_xi1T_n2(true);
    LINALG::TMatrix<TYPEBTS, dim2, 1> N2_xi2T_n2(true);
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < dim2; j++)
      {
        N2_xi1T_n2(j) += N2_xi1(i,j) * n2(i);
        N2_xi2T_n2(j) += N2_xi2(i,j) * n2(i);
      }
    }
    for (int i = 0; i < dim2; i++)
      for (int j = 0; j < dim1+dim2; j++)
        stiffc2(i,j) += -jacobi * w_gp * fp * (N2_xi1T_n2(i) * xi1_d(j) + N2_xi2T_n2(i) * xi2_d(j)) * 0.5 * (eta_b - eta_a);

    // Part IV: eta_a_d and eta_b_d (Jacobi factor contact interval)
    for (int i = 0; i < dim2; i++)
      for (int j = 0; j < dim1+dim2; j++)
        stiffc2(i,j) += -jacobi * w_gp * fp * N2T_n2(i) * 0.5 * (eta_b_d(j) - eta_a_d(j));
  }

  // -----------------------------------------------------------------
  // end: Evaluate contact stiffness of solid element


#ifdef FADCHECKSTIFFNESS

  LINALG::TMatrix<TYPEBTS, dim1, 1> fc1_FAD(true);
  LINALG::TMatrix<TYPEBTS, dim2, 1> fc2_FAD(true);

  // Evaluate contact forces at current Gauss point (fc1_FAD, fc2_FAD)
  EvaluateFcContact(fp, fc1_FAD, fc2_FAD, eta_a, eta_b, w_gp, sgn, nD, n2, N1, N2, jacobi);

  // Calculate contact stiffness (stiffc1_FAD, stiffc2_FAD) with FAD
  for (int j = 0; j < dim1+dim2; j++)
  {
    // Beam element
    for (int i = 0; i < dim1; i++)
    {
      stiffc1_FAD(i,j) += fc1_FAD(i).dx(j)
          + fc1_FAD(i).dx(dim1+dim2) * xi1_d(j)
          + fc1_FAD(i).dx(dim1+dim2+1) * xi2_d(j)
          + fc1_FAD(i).dx(dim1+dim2+2) * eta_a_d(j)
          + fc1_FAD(i).dx(dim1+dim2+3) * eta_b_d(j);
    }

    // Surface element
    for (int i = 0; i < dim2; i++)
    {
      stiffc2_FAD(i,j) += fc2_FAD(i).dx(j)
          + fc2_FAD(i).dx(dim1+dim2) * xi1_d(j)
          + fc2_FAD(i).dx(dim1+dim2+1) * xi2_d(j)
          + fc2_FAD(i).dx(dim1+dim2+2) * eta_a_d(j)
          + fc2_FAD(i).dx(dim1+dim2+3) * eta_b_d(j);
    }
  }
#endif

  return;
}
/*----------------------------------------------------------------------*
 | End: Evaluate contact stiffness
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute directional derivatives of element parameters                |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinParameter(
    const int& fixed_par,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& r1_eta,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi2)
{
  // Compute linearization of two parameters with one known parameter linearization
  //
  // NOTE: Linearization lin() = (),d * lin_d with directional derivative (),d
  //
  // Solve the following system of equations:
  //                                        _       _
  //  _                               _    | lin_xi1 |    _                _     _      _
  // | Lpar(1,1)  Lpar(1,2)  Lpar(1,3) |   |         |   | B'(1,1)  B'(1,2) |   | lin_d1 |
  // |                                 | * | lin_xi2 | = |                  | * |        |
  // |_Lpar(2,1)  Lpar(2,2)  Lpar(2,3)_|   |         |   |_B'(2,1)  B'(2,2)_|   |_lin_d2_|
  //                                       |_lin_eta_|
  //
  // This is only possible if one parameter linerization is known and that's the case here.
  //
  // Assemble needed entries of Lpar in L:
  //  _              _     _        _     _                _     _      _     _      _
  // | L(1,1)  L(1,2) |   | lin_par1 |   | B'(1,1)  B'(1,2) |   | lin_d1 |   | L(1,3) |
  // |                | * |          | = |                  | * |        | - |        | * lin_par3
  // |_L(2,1)  L(2,2)_|   |_lin_par2_|   |_B'(2,1)  B'(2,2)_|   |_lin_d2_|   |_L(2,3)_|
  //
  // Solution by inverting matrix L:
  //
  // [lin_par1; lin_par2] = L^-1 * (B' * [lin_d1; lin_d2] - [L(1,3); L(2,3)] * lin_par3) = L^-1 * B
  //
  // with B = B' * [lin_d1; lin_d2] - [L(1,3); L(2,3)] * lin_par3
  //
  // [par1_d; par2_d] = D = L^-1 * B

  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;


  // Initialize vector with parameters xi1, xi2, eta
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> par_d[3];
  par_d[0] = xi1_d;
  par_d[1] = xi2_d;
  par_d[2] = eta_d;

  // Get indices of not fixed parameters
  int par_i[2];
  int j = 0;
  for (int i = 0; i < 3; i++)
  {
    if (i != fixed_par)
    {
      par_i[j] = i;
      j++;
    }
  }

  // Clear not fixed directional derivatives
  for (int i = 0; i < 2; i++)
  {
    par_d[par_i[i]].Clear();
  }

  // Initialize matrices
  LINALG::TMatrix<TYPEBTS, 2, 2> L(true);
  LINALG::TMatrix<TYPEBTS, 2, 3> Lpar(true);
  LINALG::TMatrix<TYPEBTS, 2, 2> L_inv(true);
  LINALG::TMatrix<TYPEBTS, 2, dim1+dim2> B(true);
  LINALG::TMatrix<TYPEBTS, 2, dim1+dim2> D(true);

  // Compute Lpar for all parametes
  Lpar(0,0) = -::FADUTILS::ScalarProduct(x2_xi1, x2_xi1) + ::FADUTILS::ScalarProduct(rD, x2_xi1xi1);
  Lpar(1,0) = -::FADUTILS::ScalarProduct(x2_xi1, x2_xi2) + ::FADUTILS::ScalarProduct(rD, x2_xi2xi1);
  Lpar(0,1) = -::FADUTILS::ScalarProduct(x2_xi2, x2_xi1) + ::FADUTILS::ScalarProduct(rD, x2_xi1xi2);
  Lpar(1,1) = -::FADUTILS::ScalarProduct(x2_xi2, x2_xi2) + ::FADUTILS::ScalarProduct(rD, x2_xi2xi2);
  Lpar(0,2) = ::FADUTILS::ScalarProduct(r1_eta, x2_xi1);
  Lpar(1,2) = ::FADUTILS::ScalarProduct(r1_eta, x2_xi2);

  // Assemble needed entries of Lpar in L
  L(0,0) = Lpar(0, par_i[0]);
  L(1,0) = Lpar(1, par_i[0]);
  L(0,1) = Lpar(0, par_i[1]);
  L(1,1) = Lpar(1, par_i[1]);

  // Invert L by hand
  TYPEBTS det_L = L(0,0) * L(1,1) - L(0,1) * L(1,0);
  if (FADUTILS::CastToDouble(FADUTILS::Norm(det_L)) < DETERMINANTTOL)
  {
    dserror("ERROR: Determinant of L = 0");
  }
  L_inv(0,0) =  L(1,1) / det_L;
  L_inv(0,1) = -L(0,1) / det_L;
  L_inv(1,0) = -L(1,0) / det_L;
  L_inv(1,1) =  L(0,0) / det_L;

  // Compute B
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1; j++)
    {
      B(0,j) += -N1(i,j) * x2_xi1(i);
      B(1,j) += -N1(i,j) * x2_xi2(i);
    }
  }
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim2; j++)
    {
      B(0,j+dim1) += x2_xi1(i) * N2(i,j) - rD(i) * N2_xi1(i,j);
      B(1,j+dim1) += x2_xi2(i) * N2(i,j) - rD(i) * N2_xi2(i,j);
    }
  }

  // Evaluation at Gauss point: In this case the directional derivative of the fixed beam parameter eta is already
  // known and is now used to calculate the directional derivatives of the both surface parameters xi1 and xi2
  // NOTE: This is not implemented for fixed xi1 or xi2 yet, because for computing the linearization
  // of the contact interval borders lin_xi1 or lin_xi2 is zero (projection of surface edges)
  if (fixed_par == 2) // eta fixed
  {
    for (int i = 0; i < dim1+dim2; i++)
    {
      B(0,i) += -Lpar(0,2) * eta_d(i);
      B(1,i) += -Lpar(1,2) * eta_d(i);
    }
  }

  // Compute D = L^-1 * B
  D.Multiply(L_inv, B);

  // Finally get directional derivatives
  for (int i = 0; i < dim1+dim2; i++)
  {
    par_d[par_i[0]](i) = D(0,i);
    par_d[par_i[1]](i) = D(1,i);
  }

  xi1_d = par_d[0];
  xi2_d = par_d[1];
  eta_d = par_d[2];

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute directional derivatives of element parameters           |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute directional derivative of gap                                |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinGap(
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& gap_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d,
    const double sgn,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& nD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& n2,
    const TYPEBTS& norm_rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& r1_eta,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues>& N1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& rD_d)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Compute gap_d and rD_d
  //
  // NOTE: Linearization lin() = (),d * lin_d with directional derivative (),d
  //
  // NOTE: lin_<-gap> = -H(-gap) * lin_gap with the Heaviside function H
  //
  // lin_gap = gap_d * lin_d = n2 * rD_d * lin_d
  //
  // with rD_d = r1_eta * eta_d - x2_xi1 * xi1_d - x2_xi2 * xi2_d + [N1, -N2] and gap_d = n2 * rD_d

  // First compute directional derivative of distance vector rD
  rD_d.Clear();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1+dim2; j++)
      rD_d(i,j) += r1_eta(i) * eta_d(j) - x2_xi1(i) * xi1_d(j) - x2_xi2(i) * xi2_d(j);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1; j++)
      rD_d(i,j) += N1(i,j);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim2; j++)
      rD_d(i,j+dim1) += -N2(i,j);


  // TODO: Old/alternative linearization
  // Finally compute directional derivative of gap
  //gap_d.Clear();
  //for (int i = 0; i < 3; i++)
  //  for (int j = 0; j < dim1+dim2; j++)
  //    gap_d(j) += sgn * rD(i) * rD_d(i,j) / norm_rD;
  //std::cout << "gap_d (old linearization): " << gap_d << std::endl;


  // Finally compute directional derivative of gap (new linearization)
  // TODO: Use directly gap_d = n2 * r1_eta * eta_d + n2 * [N1, -N2] instead of
  // gap_d = n2 * rD_d using rD_d, because n2 * x2_xi1 * xi1_d = n2 * x2_xi2 * xi2_d = 0
  gap_d.Clear();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1+dim2; j++)
      gap_d(j) += n2(i) * rD_d(i,j);
  //std::cout << "gap_d (new linearization): " << gap_d << std::endl;

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute directional derivative of gap                           |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute directional derivative of surface unit normal vector         |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinNormal(
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& nD_d,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& nD,
    const TYPEBTS& norm_rD,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& n2_d,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& n2,
    const TYPEBTS& norm_a2,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& rD_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1xi2,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N2_xi2)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Compute nD_d
  //
  // NOTE: Linearization lin() = (),d * lin_d with directional derivative (),d
  //
  // lin_nD = nD_d * lin_d = auxiliary_matrix * a2_d * lin_d
  //
  // with auxiliary_matrix = (I - n2 x n2) / norm_a2

  // Initialize auxiliary_matrix
  LINALG::TMatrix<TYPEBTS, 3, 3> auxiliary_matrix(true);

  // TODO: Old/alternative linearization
  // Compute auxiliary_matrix
  for (int i = 0; i < 3; i++)
  {
    auxiliary_matrix(i,i) += 1.0 / norm_rD;
    for (int j = 0; j < 3; j++)
      auxiliary_matrix(i,j) += -nD(i) * nD(j) / norm_rD;
  }
  // Finally compute derivatives of unit distance vector
  nD_d.Clear();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < dim1+dim2; k++)
        nD_d(i,k) += auxiliary_matrix(i,j) * rD_d(j,k);
  //std::cout << "nD_d (old linearization): " << nD_d << std::endl;


  // Skew-symmetric matrix to replace cross product with matrix multiplication
  // TODO: Use a function for calculating cross products
  LINALG::TMatrix<TYPEBTS, 3, 3> x2_xi1_tilde(true);
  x2_xi1_tilde(2,1) = x2_xi1(0);
  x2_xi1_tilde(2,0) = -x2_xi1(1);
  x2_xi1_tilde(1,0) = x2_xi1(2);
  x2_xi1_tilde(1,2) = -x2_xi1_tilde(2,1);
  x2_xi1_tilde(0,2) = -x2_xi1_tilde(2,0);
  x2_xi1_tilde(0,1) = -x2_xi1_tilde(1,0);
  LINALG::TMatrix<TYPEBTS, 3, 3> x2_xi2_tilde(true);
  x2_xi2_tilde(2,1) = x2_xi2(0);
  x2_xi2_tilde(2,0) = -x2_xi2(1);
  x2_xi2_tilde(1,0) = x2_xi2(2);
  x2_xi2_tilde(1,2) = -x2_xi2_tilde(2,1);
  x2_xi2_tilde(0,2) = -x2_xi2_tilde(2,0);
  x2_xi2_tilde(0,1) = -x2_xi2_tilde(1,0);

  // First compute compute directional derivative of tangent cross product a2
  LINALG::TMatrix<TYPEBTS, 3, dim1+dim2> a2_d(true);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        a2_d(i,j) += (x2_xi1_tilde(i,k) * x2_xi2xi1(k) - x2_xi2_tilde(i,k) * x2_xi1xi1(k)) * xi1_d(j)
            + (x2_xi1_tilde(i,k) * x2_xi2xi2(k) - x2_xi2_tilde(i,k) * x2_xi1xi2(k)) * xi2_d(j);
      }
    }
  }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim2; j++)
      for (int k = 0; k < 3; k++)
        a2_d(i,j+dim1) += x2_xi1_tilde(i,k) * N2_xi2(k,j) - x2_xi2_tilde(i,k) * N2_xi1(k,j);

  // Compute auxiliary_matrix
  auxiliary_matrix.Clear();
  for (int i = 0; i < 3; i++)
  {
    auxiliary_matrix(i,i) += 1.0 / norm_a2;
    for (int j = 0; j < 3; j++)
      auxiliary_matrix(i,j) += -n2(i) * n2(j) / norm_a2;
  }

  // Finally compute directional derivative of surface unit normal vector n2 (new linearization)
  n2_d.Clear();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < dim1+dim2; k++)
        n2_d(i,k) += auxiliary_matrix(i,j) * a2_d(j,k);
  // std::cout << "n2_d (new linearization): " << n2_d << std::endl;

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute directional derivative of surface unit normal vector    |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Assemble contact force and stiffness                                 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::AssembleFcAndStiffcContact(
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 1> fc1,
    const LINALG::TMatrix<TYPEBTS, 3*numnodessol, 1> fc2,
    Epetra_Vector* fint,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol> stiffc1,
    const LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol> stiffc2,
    LINALG::SparseMatrix& stiffmatrix)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Node ids of elements
  const int* node_ids1 = element1_->NodeIds();
  const int* node_ids2 = element2_->NodeIds();

  // Temporary vectors for contact forces, DOF-GIDs and owning procs
  Epetra_SerialDenseVector fcontact1(dim1);
  Epetra_SerialDenseVector fcontact2(dim2);
  std::vector<int> lm1(dim1);
  std::vector<int> lm2(dim2);
  std::vector<int> lmowner1(dim1);
  std::vector<int> lmowner2(dim2);

  // Temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  Epetra_SerialDenseMatrix stiffcontact1(dim1, dim1+dim2);
  Epetra_SerialDenseMatrix stiffcontact2(dim2, dim1+dim2);
  std::vector<int> lmrow1(dim1);
  std::vector<int> lmrow2(dim2);
  std::vector<int> lmrowowner1(dim1);
  std::vector<int> lmrowowner2(dim2);
  std::vector<int> lmcol1(dim1+dim2);
  std::vector<int> lmcol2(dim1+dim2);


  // Assemble contact forces acting on beam and solid
  // -----------------------------------------------------------------

  // Prepare assembly

  // Fill lm1 and lmowner1
  for (int i = 0; i < numnodes; ++i)
  {
    // Get node pointer and dof ids
    DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    for (int j = 0;j < 3*numnodalvalues; ++j)
    {
      lm1[3*numnodalvalues*i + j] = NodeDofGIDs[j];
      lmowner1[3*numnodalvalues*i + j] = node->Owner();
    }
  }

  // Fill lm2 and lmowner2
  for (int i = 0; i < numnodessol; ++i)
  {
    // Get node pointer and dof ids
    DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    for (int j = 0; j < 3; ++j)
    {
      lm2[3*i + j] = NodeDofGIDs[j];
      lmowner2[3*i + j] = node->Owner();
    }
  }

  // Assemble fc1 and fc2 into global contact force vector
  for (int i = 0; i < dim1; i++)
    fcontact1[i] = FADUTILS::CastToDouble(fc1(i));
  for (int i = 0; i < dim2; i++)
    fcontact2[i] = FADUTILS::CastToDouble(fc2(i));

  LINALG::Assemble(*fint, fcontact1, lm1, lmowner1);
  LINALG::Assemble(*fint, fcontact2, lm2, lmowner2);

  // -----------------------------------------------------------------
  // End: Assemble contact forces acting on beam and solid


  // Assemble contact stiffness for beam and solid
  // -----------------------------------------------------------------

  // Prepare assembly

  // Fill lmrow1 and lmrowowner1
  lmrow1 = lm1;
  lmrowowner1 = lmowner1;

  // Fill lmrow2 and lmrowowner2
  lmrow2 = lm2;
  lmrowowner2 = lmowner2;

  // Fill lmcol1 and lmcol2
  for (int i = 0; i < numnodes; ++i)
  {
    // Get pointer and node ids
    DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    for (int j = 0; j < 3*numnodalvalues; ++j)
    {
      lmcol1[3*numnodalvalues*i + j] = NodeDofGIDs[j];
      lmcol2[3*numnodalvalues*i + j] = NodeDofGIDs[j];
    }
  }
  for (int i = 0; i < numnodessol; ++i)
  {
    // Get pointer and node ids
    DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    for (int j = 0; j < 3; ++j)
    {
      lmcol1[dim1 + 3*i + j] = NodeDofGIDs[j];
      lmcol2[dim1 + 3*i + j] = NodeDofGIDs[j];
    }
  }

  // Now finally assemble stiffc1 and stiffc2
  for (int j = 0;j < dim1+dim2; j++)
  {
    // Change sign of stiffc1 and stiffc2 due to time integration. According to analytical derivation there
    // is no minus sign, but for our time integration methods the negative stiffness must be assembled.
    for (int i = 0; i < dim1; i++)
      stiffcontact1(i,j) = -FADUTILS::CastToDouble(stiffc1(i,j));
    for (int i = 0; i < dim2; i++)
      stiffcontact2(i,j) = -FADUTILS::CastToDouble(stiffc2(i,j));
  }

  stiffmatrix.Assemble(0, stiffcontact1, lmrow1, lmrowowner1, lmcol1);
  stiffmatrix.Assemble(0, stiffcontact2, lmrow2, lmrowowner2, lmcol2);

  // ----------------------------------------------------------------
  // End: Assemble contact stiffness for beam and solid

  // Print some debug information
  const bool output = false;
  if (output)
  {
    std::cout << "element1_->Id(): " << element1_->Id() << std::endl;

    std::cout << "lmrow1: " << std::endl;
    for (int i = 0; i < dim1; i++)
      std::cout << lmrow1[i] << std::endl;

    std::cout << "lmcol1: " << std::endl;
    for (int i = 0; i < dim1+dim2; i++)
      std::cout << lmcol1[i] << std::endl;

    std::cout << "element2_->Id(): " << element2_->Id() << std::endl;

    std::cout << "lmrow2: " << std::endl;
      for (int i = 0; i < dim2; i++)
        std::cout << lmrow2[i] << std::endl;

    std::cout << "lmcol2: " << std::endl;
    for (int i = 0; i < dim1+dim2; i++)
      std::cout << lmcol2[i] << std::endl;
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Assemble contact force and stiffness                            |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Get contact interval borders                                         |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetContactIntervalBorders(
    std::vector<std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int> >& parsets)
{
  // Enable or disable console output
  const bool output = false;
  if (output)
  {
    std::cout << "GetContactIntervalBorders output:" << std::endl;
  }

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + XIETATOL;

  // Clear and resize current vector for parameter sets
  parsets.clear();
  parsets.resize(0);

  // Temporary vector containg parsets and flags for allowed or not allowed projections
  // TODO: May this can be done in a more beautiful way without using the temporary vector parsetstmp
  std::vector<std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, LINALG::TMatrix<int, 2, 1> > > parsetstmp;
  std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, LINALG::TMatrix<int, 2, 1> > parsettmp;

  // Initialize flags indicating if the projection of beam start and end point are necessary
  bool projbeamstart = true;
  bool projbeamend = true;


  // Find projection of surface edges xi1 = +-1 and xi2 = +-1
  // Only implemented for quad surface elements such as quad4, quad8
  // TODO: Implementation for tri surface elements
  // -----------------------------------------------------------------

  if (output)
  {
    std::cout << "Projection surface edges:" << std::endl;
  }

  // Fixed parameter xi1 or xi2 (i = 0: xi1 fixed, i = 1: xi2 fixed)
  for (int i = 0; i < 2; i++)
  {
    // Value for xi1 or xi2 (j = 0: xi = -1, j = 1: xi = +1)
    for (int j = 0; j < 2; j++)
    {
      // Set start point for Newton iteration
      TYPEBTS xi1 = 0.0;
      TYPEBTS xi2 = 0.0;
      TYPEBTS eta = 0.0;

      // Set value for fixed parameter xi1 or xi2
      switch(i)
      {
        case 0: xi1 = -1.0 + 2.0 * j; break;
        case 1: xi2 = -1.0 + 2.0 * j; break;
      }

      // Projection of surface edge
      bool proj_allowed = true;
      Projection(i, xi1, xi2, eta, proj_allowed);

      if (output)
      {
        std::cout << "  xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta
            << ", fixed par: " << i << ", proj_allowed: " << proj_allowed << std::endl;
      }

      // Check if the parameter eta found within this projection is already a contact interval
      // border, e.g. if the projected beam element centerline crosses an surface edge point
      bool alreadyfound = false;
      for (int k = 0; k < (int)parsetstmp.size(); k++)
      {
        if (fabs(parsetstmp[k].first(2) - eta) < XIETATOL)
        {
          alreadyfound = true;
          break;
        }
      }

      if (!alreadyfound && fabs(xi1) < limit && fabs(xi2) < limit && fabs(eta) < limit)
      {
        // Insert contact interval border parameter set
        parsettmp.first(0) = xi1;
        parsettmp.first(1) = xi2;
        parsettmp.first(2) = eta;
        parsettmp.second(0) = i;
        parsettmp.second(1) = proj_allowed;
        parsetstmp.push_back(parsettmp);

        // Check if the projection of the beam start and end point is needed
        if (fabs(eta + 1.0) < XIETATOL)
        {
          // Start point projection not needed
          projbeamstart = false;
        }
        if (fabs(eta - 1.0) < XIETATOL)
        {
          // End point projection not needed
          projbeamend = false;
        }
      }
    }
  }

  // Sort vector with parameter sets from smallest to biggest eta
  std::sort(parsetstmp.begin(), parsetstmp.end(), CompareParsets);

  // -----------------------------------------------------------------
  // End: Find projection of surface edges xi1 = +-1 and xi2 = +-1


  // Find projection of beam start and end point eta = +-1
  // -----------------------------------------------------------------

  if (output)
  {
    std::cout << "Projection of beam start and end point:" << std::endl;
  }

  // Check if the projected start point of beam (eta = -1) lies inside the surface
  if (projbeamstart)
  {
    // Set start point for newton iteration
    TYPEBTS xi1 = 0.0;
    TYPEBTS xi2 = 0.0;

    // Set beam parameter
    TYPEBTS eta = -1.0;

    // Projection of beam start point
    bool proj_allowed = true;
    Projection(2, xi1, xi2, eta, proj_allowed);

    if (fabs(xi1) < limit && fabs(xi2) < limit)
    {
      // Insert contact interval border parameter set at first position
      parsettmp.first(0) = xi1;
      parsettmp.first(1) = xi2;
      parsettmp.first(2) = eta;
      parsettmp.second(0) = 2;
      parsettmp.second(1) = proj_allowed;
      parsetstmp.insert(parsetstmp.begin(), parsettmp);
      if (output)
      {
        std::cout << "  xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << ", fixed par: " << 2
            << ", proj_allowed: " << proj_allowed << "(beam start point)" << std::endl;
      }
    }
  }

  // Check if projected end point of beam (eta = 1) lies inside the surface
  if (projbeamend)
  {
    // Set start point for newton iteration
    TYPEBTS xi1 = 0.0;
    TYPEBTS xi2 = 0.0;

    // Set beam parameter
    TYPEBTS eta = 1.0;

    // Projection of beam end point
    bool proj_allowed = true;
    Projection(2, xi1, xi2, eta, proj_allowed);

    if (fabs(xi1) < limit && fabs(xi2) < limit)
    {
      // Insert contact interval border parameter set at last position
      parsettmp.first(0) = xi1;
      parsettmp.first(1) = xi2;
      parsettmp.first(2) = eta;
      parsettmp.second(0) = 2;
      parsettmp.second(1) = proj_allowed;
      parsetstmp.push_back(parsettmp);
      if (output)
      {
        std::cout << "  xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << ", fixed par: " << 2
            << ", proj_allowed: " << proj_allowed << "(beam end point)" << std::endl;
      }
    }
  }

  // -----------------------------------------------------------------
  // End: Find projection of beam start and end point eta = +-1


  // Sort contact interval borders to identify contact intervals
  // -----------------------------------------------------------------

  // Print sorted vector with all parameter sets
  if (output)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << "Number of parameter sets found for this beam to solid contact pair: " << parsetstmp.size() << std::endl;
    if (parsetstmp.size() > 0)
    {
      for (int i = 0; i < (int)parsetstmp.size(); i++)
      {
        std::cout << "  xi1: " << parsetstmp[i].first(0) << ", xi2: " << parsetstmp[i].first(1) << ", eta: " << parsetstmp[i].first(2)
            << ", fixed par: " << parsetstmp[i].second(0) << ", proj_allowed: " << parsetstmp[i].second(1) << std::endl;
      }
    }
  }

  // Check if all found contact intervals are closed
  if (parsetstmp.size() % 2 != 0)
  {
    // TODO: Point contact
    if (parsetstmp.size() == 1)
    {
      parsetstmp.clear();
      parsetstmp.resize(0);
    }
    else
    {
      // std::cout << "Contact interval not closed" << std::endl;
      dserror("Contact interval not closed.");
    }
  }

  // Insert parsetstmp in parsets and sort out contact interval border parameter sets with not allowed projection
  std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int> parset;
  for (int i = 0; i < (int)parsetstmp.size() - 1; i += 2)
  {
    // Check if the contact interval is allowed
    if (parsetstmp[i].second(1) && parsetstmp[i+1].second(1))
    {
      for (int j = 0; j < 2; j++)
      {
        parset.first = parsetstmp[i+j].first;
        parset.second = parsetstmp[i+j].second(0);
        parsets.push_back(parset);
      }
    }
  }

  // Print final sorted vector with parameter sets
  if (output)
  {
    std::cout << "-----------------------------------------------------------------" << std::endl;
    std::cout << "Final number of parameter sets allowed for this beam to solid contact pair: " << parsets.size() << std::endl;
    if (parsets.size() > 0)
    {
      for (int i = 0; i < (int)parsets.size(); i++)
      {
        std::cout << "  xi1: " << parsets[i].first(0) << ", xi2: " << parsets[i].first(1) << ", eta: " << parsets[i].first(2)
            << ", fixed par: " << parsets[i].second << std::endl;
      }
    }
    std::cout << std::endl;
  }

  // -----------------------------------------------------------------
  // End: Sort contact interval borders to identify contact intervals

  return;
}
/*----------------------------------------------------------------------*
 | End: Get contact interval borders                                    |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Projection                                                           |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Projection(
    const int& fixed_par,
    TYPEBTS& xi1,
    TYPEBTS& xi2,
    TYPEBTS& eta,
    bool& proj_allowed)
{
  proj_allowed = true;
  bool parallel = false;

  // Initialize limit for parameter values (interval [-limit, limit])
  const double limit = 1.0 + XIETATOL;

  const bool output = false;
  if (output)
  {
    std::cout << "Projection output:" << std::endl;
    std::cout << "Start parameters xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << " and fixed_par: " << fixed_par << std::endl;
  }

  // Vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1(true);         // = N1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_eta(true);     // = N1,eta
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues> N1_etaeta(true);  // = N1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2(true);                     // = N2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1(true);                 // = N2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2(true);                 // = N2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1xi1(true);              // = N2,xi1xi1
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2xi2(true);              // = N2,xi2xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi1xi2(true);              // = N2,xi1xi2
  LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol> N2_xi2xi1(true);              // = N2,xi2xi1

  // Coords and derivatives of beam and surface element
  LINALG::TMatrix<TYPEBTS, 3, 1> r1(true);                                 // = r1
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_eta(true);                             // = r1,eta
  LINALG::TMatrix<TYPEBTS, 3, 1> r1_etaeta(true);                          // = r1,etaeta

  LINALG::TMatrix<TYPEBTS, 3, 1> x2(true);                                 // = x2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1(true);                             // = x2,xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2(true);                             // = x2,xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1xi1(true);                          // = x2,xi1xi1
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2xi2(true);                          // = x2,xi2xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi1xi2(true);                          // = x2,xi1xi2
  LINALG::TMatrix<TYPEBTS, 3, 1> x2_xi2xi1(true);                          // = x2,xi2xi1

  // Distance vector, its norm and unit distance vector
  LINALG::TMatrix<TYPEBTS, 3, 1> rD(true);                                 // = r1 - x2
  TYPEBTS norm_rD = 0.0;                                                   // = ||rD||
  LINALG::TMatrix<TYPEBTS, 3, 1> nD(true);                                 // = rD / norm_rD

  // Surface tangent cross product, its norm and unit surface normal vector
  LINALG::TMatrix<TYPEBTS, 3, 1> a2(true);                                 // = x2_xi1 x x2_xi2
  TYPEBTS norm_a2 = 0.0;                                                   // = ||a||
  LINALG::TMatrix<TYPEBTS, 3, 1> n2(true);                                 // = a / norm_a


  // Initialize function f and Jacobian J for Newton iteration
  LINALG::TMatrix<TYPEBTS,2,1> f(true);
  LINALG::TMatrix<TYPEBTS,2,3> Jpar(true);
  LINALG::TMatrix<TYPEBTS,2,2> J(true);
  LINALG::TMatrix<TYPEBTS,2,2> Jinv(true);

  // Initialize vector with parameters xi1, xi2, eta
  TYPEBTS par[3] = {xi1, xi2, eta};

  // Get indices of not fixed parameters
  int par_i[2];
  int j = 0;
  for (int i = 0; i < 3;i++)
  {
    if (i != fixed_par)
    {
      par_i[j] = i;
      j++;
    }
  }

#ifdef FADCHECKLINORTHOGONALITYCONDITION
  // Print information about current parameters
  std::cout  << std::endl << "FAD-Check: Linearization of orthogonality condition with start parameters "
      "xi1: " << xi1.val() << ", xi2: " << xi2.val() << ", eta: " << eta.val() << " and fixed_par: " << fixed_par << std::endl;
#endif

  // Initial scalar residual (L2-norm of f)
  TYPEBTS residual = 0.0;


  // Local netwon iteration
  // -----------------------------------------------------------------

  int iter;
  for (iter = 0; iter < BEAMCONTACTMAXITER; iter++)
  {
#ifdef FADCHECKLINORTHOGONALITYCONDITION
    // Set known parameters xi1, xi2, eta and element positions as primary variables for checking
    // the linearization of the orthogonality conditions (Jacobi-matrix J) borders with FAD
    BEAMCONTACT::SetFADParCoordDofs<numnodessol, numnodes, numnodalvalues>(xi1, xi2, eta);
    BEAMCONTACT::SetFADDispDofs<numnodessol, numnodes, numnodalvalues>(ele1pos_, ele2pos_, 3);
#endif

    // Update shape functions and their derivatives for beam and surface element
    GetBeamShapeFunctions(N1, N1_eta, N1_etaeta, eta);
    GetSurfShapeFunctions(N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1, xi1, xi2);

    // Update coordinates and derivatives for beam and surface element
    ComputeBeamCoordsAndDerivs(r1, r1_eta, r1_etaeta, N1, N1_eta, N1_etaeta);
    ComputeSurfCoordsAndDerivs(x2, x2_xi1, x2_xi2, x2_xi1xi1, x2_xi2xi2, x2_xi1xi2, x2_xi2xi1,
        N2, N2_xi1, N2_xi2, N2_xi1xi1, N2_xi2xi2, N2_xi1xi2, N2_xi2xi1);

    // Compute distance vector rD = r1 - x2
    for (int i = 0; i < 3; i++)
    {
      rD(i) = r1(i) - x2(i);
    }

    // Compute norm of difference vector
    norm_rD = FADUTILS::VectorNorm<3>(rD);

    // If automatic differentiation via FAD is applied, norm_rD has to be of type double since
    // this factor is needed for a pure scaling of the orthogonality conditions and has not to be linearized
    double norm_rD_scale = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(rD));

    // The closer the beam and surface element get, the smaller is norm_rD, but norm_rD is not
    // allowed to be too small, else numerical problems occur. Since in this case |eta| > 1,
    // |xi1| > 1 or |xi2| > 1 they will be sorted out later anyways.
    if (norm_rD_scale < NORMTOL)
    {
      break;
    }

    // Evaluate f at current xi1, xi2, eta
    f.Clear();
    for (int i = 0; i < 3; i++)
    {
      f(0) += rD(i) * x2_xi1(i) / norm_rD_scale;
      f(1) += rD(i) * x2_xi2(i) / norm_rD_scale;
    }

    // Compute scalar residuum
    residual = FADUTILS::VectorNorm<2>(f);

    // Reset matrices
    J.Clear();
    Jpar.Clear();
    Jinv.Clear();

    // Evaluate derivatives of f at current xi1, xi2, eta
    for(int i = 0; i < 3; i++)
    {
      if (fixed_par != 0) // xi1 is not fixed
      {
        // xi1 derivate of f
        Jpar(0,0) += (-x2_xi1(i) * x2_xi1(i) + rD(i) * x2_xi1xi1(i)) / norm_rD_scale;
        Jpar(1,0) += (-x2_xi1(i) * x2_xi2(i) + rD(i) * x2_xi2xi1(i)) / norm_rD_scale;
      }
      if (fixed_par != 1) // xi2 is not fixed
      {
        // xi2 derivate of f
        Jpar(0,1) += (-x2_xi2(i) * x2_xi1(i) + rD(i) * x2_xi1xi2(i)) / norm_rD_scale;
        Jpar(1,1) += (-x2_xi2(i) * x2_xi2(i) + rD(i) * x2_xi2xi2(i)) / norm_rD_scale;
      }
      if (fixed_par != 2) // eta is not fixed
      {
        // eta derivate of f
        Jpar(0,2) += r1_eta(i) * x2_xi1(i) / norm_rD_scale;
        Jpar(1,2) += r1_eta(i) * x2_xi2(i) / norm_rD_scale;
      }
    }

    // Assemble needed derivatives of f in J
    J(0,0) = Jpar(0,par_i[0]);
    J(1,0) = Jpar(1,par_i[0]);
    J(0,1) = Jpar(0,par_i[1]);
    J(1,1) = Jpar(1,par_i[1]);

#ifdef FADCHECKLINORTHOGONALITYCONDITION
    // Print information about current iteration
    std::cout << "Local Newton iteration " << iter + 1 << ":" << std::endl;

    // Check linearization of orthogonality conditions with FAD
    LINALG::TMatrix<TYPEBTS, 2, 2> J_FAD(true);
    FADCheckLinOrthogonalityCondition(fixed_par, rD, x2_xi1, x2_xi2, J_FAD, J);
#endif

    // Inverting 2x2-matrix J by hard coded formula, so that it is possible
    // to handle colinear vectors, because they lead to det(J) = 0
    TYPEBTS det_J = J(0,0) * J(1,1) - J(1,0) * J(0,1);

    // If det_J = 0 we assume, that the beam centerline and the surface edge are parallel.
    // These projection is not needed due the fact that the contact interval can also be
    // identified by two contact interval borders found with the GetContactLines method
    parallel = FADUTILS::CastToDouble(FADUTILS::Norm(det_J)) < COLINEARTOL;

    // Check if the local Newton iteration has converged
    // If the start point fulfills the orthogonalty conditions (residual < BEAMCONTACTTOL), we also check if
    // the beam centerline and the surface edge are parallel. This is done by calculating det_J before checking
    // if the local Newton iteration has converged by fulfilling the condition residual < BEAMCONTACTTOL
    if (FADUTILS::CastToDouble(residual) < BEAMCONTACTTOL && !parallel)
    {
      if (output)
      {
        std::cout << "Local Newton iteration converged after " << iter << " iterations" << std::endl;
        std::cout << "Found point at xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << " with residual: " << residual << std::endl;
      }
      // Local Newton iteration converged
      break;
    }
    else if (output && iter > 0)
    {
      std::cout << "New point at xi1: " << xi1 << ", xi2: " << xi2 << ", eta: " << eta << " with residual: " << residual <<  std::endl;
    }

    // Singular J
    if (parallel)
    {
      // Sort out
      if (output)
      {
        std::cout << "elementscolinear: det_J = " << FADUTILS::CastToDouble(FADUTILS::Norm(det_J)) << std::endl;
      }
      break;
    }
    // Regular J (inversion possible)
    else
    {
      // Do not sort out

      // Invert J
      Jinv(0,0) = J(1,1) / det_J;
      Jinv(0,1) = -J(0,1) / det_J;
      Jinv(1,0) = -J(1,0) / det_J;
      Jinv(1,1) = J(0,0) / det_J;
    }

    par[par_i[0]] += -Jinv(0,0) * f(0) - Jinv(0,1) * f(1);
    par[par_i[1]] += -Jinv(1,0) * f(0) - Jinv(1,1) * f(1);

    xi1 = par[0];
    xi2 = par[1];
    eta = par[2];
  }
  // -----------------------------------------------------------------
  // End: Local Newton iteration


  // Local Newton iteration unconverged after BEAMCONTACTMAXITER
  if (residual > BEAMCONTACTTOL || parallel)
  {
    par[par_i[0]] = 1e+12;
    par[par_i[1]] = 1e+12;

    xi1 = par[0];
    xi2 = par[1];
    eta = par[2];

    if (output)
      std::cout << "Local Newton iteration unconverged (!) after " << iter + 1 << " iterations" << std::endl;
  }
  else if (fabs(xi1) < limit && fabs(xi2) < limit && fabs(eta) < limit)
  {
    LINALG::TMatrix<TYPEBTS, 3, 1> n2(true);
    double sgn = 1.0;

    // Compute distance vector rD, its norm norm_rD and unit distance vector nD
    ComputeDistanceNormal(r1, x2, rD, norm_rD, nD);

    // Compute surface tangent cross product a2, its norm_a2 and surface unit normal vector n2
    ComputeSurfaceNormal(x2_xi1, x2_xi2, a2, norm_a2, n2);

    // Compute sign of projection of distance vector rD on surface unit normal vector n2
    sgn = FADUTILS::CastToDouble(FADUTILS::Signum(FADUTILS::ScalarProduct(nD, n2)));

    if (sgn < 0)
    {
      // If there is no information about the unit distance vector nD of the last time step (e.g. in the first call
      // after the contact pair element has been created), there is no change of the dircetion of nD
      bool normaldir_changed = false;

      // Otherwise check if the direction of nD has changed as follows
      if (normalsets_old_.size() > 0)
      {
        // Find the nearest eta from the last time step to the current eta
        int nearest_index = 0;
        if (normalsets_old_.size() > 1)
        {
          int lower_index = std::lower_bound(normalsets_old_.begin(), normalsets_old_.end(), eta, CompareNormalsets) - normalsets_old_.begin();
          lower_index = std::max(1, std::min(lower_index, (int)normalsets_old_.size() - 1));
          if (fabs(normalsets_old_[lower_index].first - eta) < fabs(normalsets_old_[lower_index-1].first - eta))
          {
            nearest_index = lower_index;
          }
          else
          {
            nearest_index = lower_index - 1;
          }
        }

        // Check if the direction of nD has changed
        if (FADUTILS::ScalarProduct(nD, normalsets_old_[nearest_index].second) < 0)
        {
          normaldir_changed = true;
        }
      }

      // If there is no a change of the direction of the normal vector, the projection is not allowed
      // and the found parameters are not treated in the later caluclation
      if (!normaldir_changed)
      {
        proj_allowed = false;
      }
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Projection                                                      |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate beam shape functions and derivatives                        |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetBeamShapeFunctions(
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
  const DRT::Element::DiscretizationType distype = element1_->Shape();

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
    if (element1_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance())
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    double length = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();

    // Get values and derivatives of shape functions
    DRT::UTILS::shape_function_hermite_1D(N_i, eta, length, distype);
    DRT::UTILS::shape_function_hermite_1D_deriv1(N_i_eta, eta, length, distype);
    DRT::UTILS::shape_function_hermite_1D_deriv2(N_i_etaeta, eta, length, distype);
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  // Assemble the individual shape functions in matrices, such that: r = N * d, r_eta = N_eta * d, r_etaeta = N_etaeta * d
  AssembleBeamShapefunctions(N_i, N_i_eta, N_i_etaeta, N, N_eta, N_etaeta);

  return;
}
/*----------------------------------------------------------------------*
 | End: Evaluate beam shape functions and derivatives                   |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Evaluate solid surface shape functions and derivatives               |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetSurfShapeFunctions(
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi1xi1,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi2xi2,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi1xi2,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodessol>& N_xi2xi1,
    const TYPEBTS& xi1,
    const TYPEBTS& xi2)
{
  // Clear shape functions and derivatives
  N.Clear();
  N_xi1.Clear();
  N_xi2.Clear();
  N_xi1xi1.Clear();
  N_xi2xi2.Clear();
  N_xi1xi2.Clear();
  N_xi2xi1.Clear();

  LINALG::TMatrix<TYPEBTS, 1, numnodessol> N_i(true);
  LINALG::TMatrix<TYPEBTS, 2, numnodessol> N_i_xi(true);
  LINALG::TMatrix<TYPEBTS, 3, numnodessol> N_i_xixi(true);

  DRT::UTILS::shape_function_2D(N_i, xi1, xi2, element2_->Shape());
  DRT::UTILS::shape_function_2D_deriv1(N_i_xi, xi1, xi2, element2_->Shape());
  DRT::UTILS::shape_function_2D_deriv2(N_i_xixi, xi1, xi2, element2_->Shape());

  // Assemble the individual shape functions in matrices, such that: r = N * d, r_xi = N_xi * d, r_xi1xi1 = N_xi1xi1 * d, ...
  AssembleSurfShapefunctions(N_i, N_i_xi, N_i_xixi, N, N_xi1, N_xi2, N_xi1xi1, N_xi2xi2, N_xi1xi2, N_xi2xi1);

  return;
}
/*----------------------------------------------------------------------*
 | End: Evaluate solid surface shape functions and derivatives          |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Assemble beam shape functions                                        |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::AssembleBeamShapefunctions(
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
  for (int i = 0;i < 3*numnodes*numnodalvalues; i++)
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

  for (int i = 0; i < numnodes*numnodalvalues; i++)
  {
    assembly_N[0][3*i] = i + 1;
    assembly_N[1][3*i + 1] = i + 1;
    assembly_N[2][3*i + 2] = i + 1;
  }

  // Assemble the matrices of the shape functions
  for (int i = 0; i < 3*numnodes*numnodalvalues; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if(assembly_N[j][i]==0)
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
  return;
}
/*----------------------------------------------------------------------*
 | End: Assemble beam shape functions                                   |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Assemble solid surface shape functions                               |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::AssembleSurfShapefunctions(
    const LINALG::TMatrix<TYPEBTS,1,numnodessol>& N_i,
    const LINALG::TMatrix<TYPEBTS,2,numnodessol>& N_i_xi,
    const LINALG::TMatrix<TYPEBTS,3,numnodessol>& N_i_xixi,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2xi1)
{
  // assembly_N is just an array to help assemble the matrices of the shape functions
  // it determines, which shape function is used in which column of N
  int assembly_N[3][3*numnodessol];

  // Initialize to zero
  for (int i = 0; i < 3*numnodessol; i++)
    for (int j = 0; j < 3; j++)
      assembly_N[j][i] = 0.0;


  // Set number of shape functions for each 3x3 block:
  // e.g. quad4 surface element (numnodesol = 4)
  // int assembly_N[3][12] = {{1,0,0,2,0,0,3,0,0,4,0,0},
  //                          {0,1,0,0,2,0,0,3,0,0,4,0},
  //                          {0,0,1,0,0,2,0,0,3,0,0,4}};

  for (int i = 0; i < numnodessol; i++)
  {
    assembly_N[0][3*i] = i + 1;
    assembly_N[1][3*i + 1] = i + 1;
    assembly_N[2][3*i + 2] = i + 1;
  }

  // Assemble the matrices of the shape functions
  for (int i = 0; i < 3*numnodessol; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if(assembly_N[j][i] == 0)
      {
        N(j,i) = 0;
        N_xi1(j,i) = 0;
        N_xi2(j,i) = 0;
        N_xi1xi1(j,i) = 0;
        N_xi2xi2(j,i) = 0;
        N_xi1xi2(j,i) = 0;
        N_xi2xi1(j,i) = 0;
      }
      else
      {
        int k = assembly_N[j][i] - 1;
        N(j,i) = N_i(k);
        N_xi1(j,i) = N_i_xi(0,k);
        N_xi2(j,i) = N_i_xi(1,k);
        N_xi1xi1(j,i) = N_i_xixi(0,k);
        N_xi2xi2(j,i) = N_i_xixi(1,k);
        N_xi1xi2(j,i) = N_i_xixi(2,k);
        N_xi2xi1(j,i) = N_i_xixi(2,k); // = N_xi1xi2
      }
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Assemble solid surface shape functions                          |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute beam point coordinates and their derivatives                 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeBeamCoordsAndDerivs(
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
    for (int j = 0; j < 3*numnodes*numnodalvalues; j++)
    {
      r(i) += N(i,j) * ele1pos_(j);
      r_eta(i) += N_eta(i,j) * ele1pos_(j);
      r_etaeta(i) += N_etaeta(i,j) * ele1pos_(j);
    }
  }

  return;

}
/*----------------------------------------------------------------------*
 | End: Compute beam contact point coordinates and their derivative     |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute solid surface point coordinates and their derivatives        |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeSurfCoordsAndDerivs(
    LINALG::TMatrix<TYPEBTS,3,1>& r,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi1,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi1xi1,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi2xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi1xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& r_xi2xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1xi1,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi1xi2,
    LINALG::TMatrix<TYPEBTS,3,3*numnodessol>& N_xi2xi1)
{
  r.Clear();
  r_xi1.Clear();
  r_xi2.Clear();
  r_xi1xi1.Clear();
  r_xi2xi2.Clear();
  r_xi1xi2.Clear();
  r_xi2xi1.Clear();

  // Compute output variable
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3*numnodessol; j++)
    {
      r(i) += N(i,j) * ele2pos_(j);
      r_xi1(i) += N_xi1(i,j) * ele2pos_(j);
      r_xi2(i) += N_xi2(i,j) * ele2pos_(j);
      r_xi1xi1(i) += N_xi1xi1(i,j) * ele2pos_(j);
      r_xi2xi2(i) += N_xi2xi2(i,j) * ele2pos_(j);
      r_xi1xi2(i) += N_xi1xi2(i,j) * ele2pos_(j);
      r_xi2xi1(i) += N_xi2xi1(i,j) * ele2pos_(j);
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute solid surface point coordinates and their derivatives   |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute distance vector                                              |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeDistanceNormal(
    const LINALG::TMatrix<TYPEBTS,3,1>& r1,
    const LINALG::TMatrix<TYPEBTS,3,1>& x2,
    LINALG::TMatrix<TYPEBTS,3,1>& rD,
    TYPEBTS& norm_rD,
    LINALG::TMatrix<TYPEBTS,3,1>& nD)
{
  // Reset variables
  rD.Clear();
  norm_rD = 0.0;
  nD.Clear();

  // Compute distance vector rD = r1 - x2 (not normalized)
  for (int i = 0; i < 3; i++)
  {
    rD(i) = r1(i) - x2(i);
  }

  // Compute norm of distance vector
  norm_rD = FADUTILS::VectorNorm<3>(rD);

  if (FADUTILS::CastToDouble(norm_rD) < NORMTOL)
  {
    dserror("ERROR: Distance vector of length zero! --> Change time step!");
  }

  // Compute unit distance vector
  for (int i = 0; i < 3; i++)
  {
    nD(i) = rD(i) / norm_rD;
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute distance vector                                         |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Compute surface normal vector                                        |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeSurfaceNormal(
    const LINALG::TMatrix<TYPEBTS,3,1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS,3,1>& x2_xi2,
    LINALG::TMatrix<TYPEBTS,3,1>& a2,
    TYPEBTS& norm_a2,
    LINALG::TMatrix<TYPEBTS,3,1>& n2)
{
  // Reset variables
  n2.Clear();
  norm_a2 = 0;
  a2.Clear();

  // Compute surface normal vector (not normalized)
  a2(0) = x2_xi1(1)*x2_xi2(2) - x2_xi1(2)*x2_xi2(1);
  a2(1) = x2_xi1(2)*x2_xi2(0) - x2_xi1(0)*x2_xi2(2);
  a2(2) = x2_xi1(0)*x2_xi2(1) - x2_xi1(1)*x2_xi2(0);

  // Compute norm of surface normal vector
  norm_a2 = FADUTILS::VectorNorm<3>(a2);

  if (FADUTILS::CastToDouble(norm_a2) < NORMTOL)
  {
    dserror("ERROR: Surface normal vector of length zero! --> Change time step!");
  }

  // Compute unit surface normal vector
  for (int i = 0; i < 3; i++)
  {
    n2(i) = a2(i) / norm_a2;
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: Compute surface normal vector                                   |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Check if contact is active or inactive                                |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::CheckContactStatus(
    const double& pp,
    const TYPEBTS& gap,
    bool& contactflag)
{

#ifdef LINPENALTY
  // Linear penalty force law
  if (gap < 0)
    contactflag = true;
  else
    contactflag = false;
#endif

#ifdef QUADPENALTY
  // Linear penalty force law
  if (gap < 0)
    contactflag = true;
  else
    contactflag = false;
#endif

#ifdef ARBITPENALTY
  const double g0 = G0;

  // Linear penalty force law
  if(ARBITPENALTY == 1 or ARBITPENALTY == 3 or ARBITPENALTY == 4 or ARBITPENALTY == 5)
  {
    if (gap < g0)
      contactflag = true;
    else
      contactflag = false;
  }
  else if(ARBITPENALTY == 2)
  {
    if (gap < 0)
      contactflag = true;
    else
      contactflag = false;
  }
#endif

  return;
}
/*----------------------------------------------------------------------*
 | End: Check if contact is active or inactive                           |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Get global dofs of a node                                meier 02/14 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
std::vector<int> CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetGlobalDofs(const DRT::Node* node)
{
  // Get dofs in beam contact discretization
  const std::vector<int> cdofs = ContactDiscret().Dof(node);

  // Get dofs in problem discretization via offset
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k = 0; k < (int)(cdofs.size()); ++k)
  {
    pdofs[k] = (dofoffsetmap_.find(cdofs[k]))->second;
  }

  return pdofs;
}
/*----------------------------------------------------------------------*
 | End: Get global dofs of a node                                       |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Change the sign of the normal vector                   meier 02/2014 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::InvertNormal()
{
//  for (int i=0; i<3;i++)
//    normal_(i) = -normal_(i);
}
/*----------------------------------------------------------------------*
 | End: Change the sign of the old normal vector                        |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Update all class variables at the end of time step                   |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateClassVariablesStep()
{
  // Update history variables
  const bool output = false;
  if (output)
  {
    std::cout << "UpdateClassVariablesStep:" << std::endl;
    for (int i = 0; i < (int)normalsets_.size(); i++)
    {
      std::cout << "normalset " << i << ":" << std::endl;
      std::cout << "  eta: " << normalsets_[i].first << std::endl;
      std::cout << "  nD: " << normalsets_[i].second(0) << ", " << normalsets_[i].second(1) << ", " << normalsets_[i].second(2)<< std::endl;
    }
  }
  normalsets_old_ = normalsets_;
}
/*----------------------------------------------------------------------*
 | End: Update all class variables at the end of time step              |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Shift current normal vector to old normal vector       meier 02/2014 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ShiftNormal()
{
//  for (int j=0;j<3;j++)
//    normal_old_(j) = normal_(j);
}
/*----------------------------------------------------------------------*
 | End: Shift current normal vector to old normal vector                |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Check if there is a difference of old and new gap      meier 02/2014 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
bool CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetNewGapStatus()
{
//  TYPEBTS gap_diff = gap_-gap_original_;
//
//  if (FADUTILS::CastToDouble(FADUTILS::Norm(gap_diff)) < GAPTOL)
//    return false;
//  else
    return true;
}
/*----------------------------------------------------------------------*
 | End: Check if there is a difference of old and new gap               |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Update nodal coordinates (public)                        meier 02/14 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateElePos(
    Epetra_SerialDenseMatrix& newele1pos,
    Epetra_SerialDenseMatrix& newele2pos)
{
  // Beam element positions
  for (int i = 0; i < 3*numnodalvalues; i++)
    for (int j = 0; j < numnodes; j++)
      ele1pos_(3*numnodalvalues*j+i) = newele1pos(i,j);

  // Solid element positions
  for (int i = 0; i < 3; i++)               // Loop over nodal dofs
    for (int j = 0; j < numnodessol; j++)   // Loop over nodes
      ele2pos_(3*j+i) = newele2pos(i,j);

  return;
}
/*----------------------------------------------------------------------*
 | End: Update nodal coordinates (public)                               |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Update nodal tangents for tangent smoothing (public)     meier 02/14 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateEleSmoothTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{
//  //Tangent smoothing is only possible for Reissner beam elements --> dserror() otherwise
//  if (numnodalvalues>1)
//    dserror("Tangent smoothing only possible for Reissner beam elements (numnodalvalues=1)!!!");
//
//  LINALG::Matrix<3*numnodes,1> elepos_aux(true);
//  //Tangent smoothing only possible with data type double (not with Sacado FAD)
//  for (int i=0;i<3*numnodes;i++)
//    elepos_aux(i)=FADUTILS::CastToDouble(ele1pos_(i));
//
//  nodaltangentssmooth1_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element1_,neighbors1_);
//
//  elepos_aux.Clear();
//  //Tangent smoothing only possible with data type double (not with Sacado FAD)
//  for (int i=0;i<3*numnodes;i++)
//    elepos_aux(i)=FADUTILS::CastToDouble(ele2pos_(i));
//
//  nodaltangentssmooth2_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element2_,neighbors2_);
}
/*----------------------------------------------------------------------*
 | End: Update nodal coordinates (public)                               |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | Shift Nodal positions (public)                           meier 02/14 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ShiftNodalPositions()
{
//  //Reissner beams
//  if (numnodalvalues == 1)
//  {
//    for (int i=0; i<numnodes; i++)
//    {
//      for (int j=0;j<3;j++)
//      {
//        ele1pos_(3*i + j) = ele1pos_(3*i + j) + SHIFTVALUE * normal_old_(j);
//      }
//    }
//  }
//  //Kirchhoff beams
//  else if (numnodalvalues == 2)
//  {
//    if (numnodes == 2)
//    {
//      for (int j=0;j<3;j++)
//      {
//        ele1pos_(j) = ele1pos_(j) + SHIFTVALUE * normal_old_(j);
//        ele1pos_(6+j) = ele1pos_(6+j) + SHIFTVALUE * normal_old_(j);
//      }
//    }
//    else
//    {
//      dserror("Only numnodes = 2 possible for Kirchhoff beams!!!");
//    }
//  }
//  else
//  {
//    dserror("The parameter numnodalvalues can only have the values 1 or 2!!!");
//  }
  return;
}
/*----------------------------------------------------------------------*
 | End: Shift Nodal positions (public)                                  |
 *----------------------------------------------------------------------*/


Teuchos::RCP<CONTACT::Beam3tosolidcontactinterface> CONTACT::Beam3tosolidcontactinterface::Impl(
    const int numnodessol,
    const int numnodes,
    const int numnodalvalues,
    const DRT::Discretization& pdiscret,
    const DRT::Discretization& cdiscret,
    const std::map<int,int>& dofoffsetmap,
    DRT::Element* element1,
    DRT::Element* element2,
    Teuchos::ParameterList beamcontactparams)
{

  if (numnodalvalues!=1 and numnodalvalues!=2)
    dserror("Only the values 1 and 2 are valid for numnodalvalues!");

  if (numnodalvalues!=2 and numnodes!=2)
    dserror("Only the values numnodes=2 is possible for Kirchhoff beams, i.e. if numnodalvalues=2!");

  if (numnodes!=2 and numnodes!=3 and numnodes!=4 and numnodes!=5)
    dserror("Only the values 2, 3, 4 and 5 are valid for numnodes!");

  if (numnodessol!=3 and numnodessol!=6 and numnodessol!=4 and numnodessol!=8 and numnodessol!=9)
    dserror("Only the values 3, 4, 6, 8 and 9 are valid for numnodessol!");


  switch (numnodessol)
  {
//    case 3:
//    {
//      switch (numnodalvalues)
//      {
//        case 1:
//        {
//          switch (numnodes)
//          {
//            case 2:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 3:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 4:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 5:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//          }
//          break;
//        }
//        case 2:
//        {
//          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//        }
//      }
//      break;
//    }
    case 4:
    {
      switch (numnodalvalues)
      {
//        case 1:
//        {
//          switch (numnodes)
//          {
//            case 2:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 3:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 4:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 5:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//          }
//          break;
//        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
//    case 6:
//    {
//      switch (numnodalvalues)
//      {
//        case 1:
//        {
//          switch (numnodes)
//          {
//            case 2:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 3:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 4:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 5:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//          }
//          break;
//        }
//        case 2:
//        {
//          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//        }
//      }
//      break;
//    }
    case 8:
    {
      switch (numnodalvalues)
      {
//        case 1:
//        {
//          switch (numnodes)
//          {
//            case 2:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 3:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 4:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 5:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//          }
//          break;
//        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
//    case 9:
//    {
//      switch (numnodalvalues)
//      {
//        case 1:
//        {
//          switch (numnodes)
//          {
//            case 2:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 3:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 4:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//            case 5:
//            {
//              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//            }
//          }
//          break;
//        }
//        case 2:
//        {
//          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
//        }
//      }
//      break;
//    }
  }

  return Teuchos::null;
}

#if defined(FADCHECKLINCONTACTINTERVALBORDER) || defined(FADCHECKLINGAUSSPOINT)
/*----------------------------------------------------------------------*
 | FAD-Check for linearization of element parameters                    |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinParameter(
    const int& fixed_par,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d_FAD,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d_FAD,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d_FAD,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Initialize array with directional derivatives of parameters xi1, xi2, eta for FAD
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> par_d_FAD[3];
  par_d_FAD[0] = xi1_d_FAD;
  par_d_FAD[1] = xi2_d_FAD;
  par_d_FAD[2] = eta_d_FAD;

  // Get indices of not fixed parameters
  int par_i[2];
  int j = 0;
  for (int i = 0; i < 3; i++)
  {
    if (i != fixed_par)
    {
      par_i[j] = i;
      j++;
    }
  }

  // Clear not fixed directional derivatives
  for (int i = 0; i < 2; i++)
  {
    par_d_FAD[par_i[i]].Clear();
  }

  // Compute norm of distance vector rD to scale the equations (this yields better conditioning)
  // NOTE: Even if automatic differentiation via FAD is applied, norm_rD has to be of type double since this
  // factor is needed for a pure scaling of the nonlinear orthogonality conditions and has not to be linearized!
  double norm_rD_scale = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(rD));

  // Evaluate f of orthogonality conditions
  LINALG::TMatrix<TYPEBTS, 2, 1> f(true);
  for (int i = 0; i < 3; i++)
  {
    f(0) += rD(i) * x2_xi1(i) / norm_rD_scale;
    f(1) += rD(i) * x2_xi2(i) / norm_rD_scale;
  }

  // Initialize matrices of system of equations showed in method ComputeLinParameter
  LINALG::TMatrix<TYPEBTS, 2, 2> L(true);
  LINALG::TMatrix<TYPEBTS, 2, 2> L_inv(true);
  LINALG::TMatrix<TYPEBTS, 2, dim1+dim2> B(true);
  LINALG::TMatrix<TYPEBTS, 2, dim1+dim2> D(true);

  // Compute L elementwise for the not fixed parameters
  L(0,0) = f(0).dx(dim1 + dim2 + par_i[0]);
  L(0,1) = f(0).dx(dim1 + dim2 + par_i[1]);
  L(1,0) = f(1).dx(dim1 + dim2 + par_i[0]);
  L(1,1) = f(1).dx(dim1 + dim2 + par_i[1]);

  // Invert L by hand
  TYPEBTS det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
  if (FADUTILS::CastToDouble(FADUTILS::Norm(det_L)) < DETERMINANTTOL)
    dserror("ERROR: Determinant of L = 0");
  L_inv(0,0) = L(1,1) / det_L;
  L_inv(0,1) = -L(0,1) / det_L;
  L_inv(1,0) = -L(1,0) / det_L;
  L_inv(1,1) = L(0,0) / det_L;

  // Compute B
  for (int j = 0; j < dim1+dim2; j++)
  {
    B(0,j) = -f(0).dx(j);
    B(1,j) = -f(1).dx(j);
  }
  if (fixed_par == 2) // eta fixed (eta_d_FAD is known)
  {
    for (int j = 0; j < dim1+dim2; j++)
    {
      B(0,j) += -f(0).dx(dim1 + dim2 + 2) * eta_d_FAD(j);
      B(1,j) += -f(1).dx(dim1 + dim2 + 2) * eta_d_FAD(j);
    }
  }

  // Compute D = L^-1 * B
  D.Multiply(L_inv, B);

  // Finally compute directional derivatives of not fixed parameters
  for (int i = 0; i < dim1+dim2; i++)
  {
    par_d_FAD[par_i[0]](i) = D(0,i);
    par_d_FAD[par_i[1]](i) = D(1,i);
  }

  // Store directional derivatives parameters in original variables
  xi1_d_FAD = par_d_FAD[0];
  xi2_d_FAD = par_d_FAD[1];
  eta_d_FAD = par_d_FAD[2];

  // Store analytical calculated directional derivatives of parameters xi1, xi2, eta in an array for comparison
  LINALG::TMatrix<TYPEBTS, dim1+dim2, 1> par_d[3];
  par_d[0] = xi1_d;
  par_d[1] = xi2_d;
  par_d[2] = eta_d;

  // Compare linearizations
  for (int i = 0; i < 2; i++)
  {
    // Print name of parameter
    switch (par_i[i])
    {
      case 0:
        std::cout << "xi1_d: " << std::endl;
        break;
      case 1:
        std::cout << "xi2_d: " << std::endl;
        break;
      case 2:
        std::cout << "eta_d: " << std::endl;
        break;
    }

    // Print directional derivatives
    std::cout << "Lin: " << std::endl;
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << par_d[par_i[i]](j).val() << "  ";
    std::cout << std::endl;

    // Print directional derivatives calculated via FAD
    std::cout << "FAD: " << std::endl;
    for (int j = 0; j < dim1+dim2; j++)
    {
      std::cout << par_d_FAD[par_i[i]](j).val() << "  ";
    }
    std::cout << std::endl;

    // Calculate maximal difference between par_d and par_d_FAD
    const double tolerance = 1.0e-12;
    double diff = 0.0;
    for (int j = 0; j < dim1+dim2; j++)
    {
      diff = std::max(diff, fabs(par_d[par_i[i]](j).val() - par_d_FAD[par_i[i]](j).val()));
    }
    std::cout << "Maximal difference: " << diff;
    if (diff > tolerance)
    {
      std::cout << " > " << tolerance << " (tolerance)";
    }
    std::cout << std::endl;
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: FAD-Check for linearization of element parameters               |
 *----------------------------------------------------------------------*/
#endif

#ifdef FADCHECKLINGAUSSPOINT
/*----------------------------------------------------------------------*
 | FAD-Check for linearization of gap and distance vector               |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinGapAndDistanceVector(
    const TYPEBTS& gap,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& gap_d_FAD,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& rD_d_FAD,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& gap_d,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& rD_d)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Compute directional derivative of gap with FAD
  for (int i=0;i<dim1+dim2;i++)
  {
    gap_d_FAD(i) = gap.dx(i)
      + gap.dx(dim1+dim2) * xi1_d(i)
      + gap.dx(dim1+dim2+1) * xi2_d(i)
      + gap.dx(dim1+dim2+2) * eta_d(i);
  }

  // Compute directional derivative of distance vector rD with FAD
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j <dim1+dim2; j++)
    {
      rD_d_FAD(i,j) = rD(i).dx(j)
          + rD(i).dx(dim1+dim2) * xi1_d(j)
          + rD(i).dx(dim1+dim2+1) * xi2_d(j)
          + rD(i).dx(dim1+dim2+2) * eta_d(j);
    }
  }

  // Print name of directional derivative of gap
  std::cout << "gap_d: " << std::endl;

  // Print directional derivative
  std::cout << "Lin: " << std::endl;
  for (int i = 0; i < dim1+dim2; i++)
    std::cout << gap_d(i).val() << "  ";
  std::cout << std::endl;
  // Print directional derivative calculated via FAD
  std::cout << "FAD: " << std::endl;
  for (int i = 0; i < dim1+dim2; i++)
    std::cout << gap_d_FAD(i).val() << "  ";
  std::cout << std::endl;

  const double tolerance = 1.0e-12;

  // Calculate maximal difference between gap_d and gap_d_FAD
  double diff = 0.0;
  for (int i = 0; i < dim1+dim2; i++)
    diff = std::max(diff, fabs(gap_d(i).val() - gap_d_FAD(i).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;


  // Print name of directional derivative of distance vector rD
  std::cout << "rD_d: " << std::endl;

  // Print directional derivative
  std::cout << "Lin: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << rD_d(i,j).val() << "  ";
    std::cout << std::endl;
  }
  // Print directional derivative calculated via FAD
  std::cout << "FAD: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << rD_d_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Calculate maximal difference between rD_d and rD_d_FAD
  diff = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1+dim2; j++)
      diff = std::max(diff, fabs(rD_d(i,j).val() - rD_d_FAD(i,j).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;
}
/*----------------------------------------------------------------------*
 | End: FAD-Check for linearization of gap and distance vector          |
 *----------------------------------------------------------------------*/
#endif

#ifdef FADCHECKLINGAUSSPOINT
/*----------------------------------------------------------------------*
 | FAD-Check for linearization of normal vector                         |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinNormal(
    const LINALG::TMatrix<TYPEBTS, 3, 1>& nD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& n2,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi1_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& xi2_d,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues+3*numnodessol, 1>& eta_d,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& nD_d_FAD,
    LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& n2_d_FAD,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& nD_d,
    const LINALG::TMatrix<TYPEBTS, 3, 3*numnodes*numnodalvalues+3*numnodessol>& n2_d)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Calculate directional derivative of unit distance vector nD with FAD
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
    {
      nD_d_FAD(i,j) = nD(i).dx(j)
          + nD(i).dx(dim1+dim2) * xi1_d(j)
          + nD(i).dx(dim1+dim2+1) * xi2_d(j)
          + nD(i).dx(dim1+dim2+2) * eta_d(j);
    }
  }

  // Calculate directional derivative of unit surface normal vector n2 with FAD
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
    {
      n2_d_FAD(i,j) = n2(i).dx(j)
          + n2(i).dx(dim1+dim2) * xi1_d(j)
          + n2(i).dx(dim1+dim2+1) * xi2_d(j)
          + n2(i).dx(dim1+dim2+2) * eta_d(j);
    }
  }

  // Print name of directional derivative of unit distance vector nD
  std::cout << "nD_d: " << std::endl;

  // Print directional derivative
  std::cout << "Lin: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << nD_d(i,j).val() << "  ";
    std::cout << std::endl;
  }
  // Print directional derivative calculated via FAD
  std::cout << "FAD: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << nD_d_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  const double tolerance = 1.0e-12;

  // Calculate maximal difference between nD_d and nD_d_FAD
  double diff = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1+dim2; j++)
      diff = std::max(diff, fabs(nD_d(i,j).val() - nD_d_FAD(i,j).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;


  // Print name of directional derivative of unit surface normal vector n2
  std::cout << "n2_d: " << std::endl;

  // Print directional derivative
  std::cout << "Lin: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << n2_d(i,j).val() << "  ";
    std::cout << std::endl;
  }
  // Print directional derivative calculated via FAD
  std::cout << "FAD: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
      std::cout << n2_d_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Calculate maximal difference between n2_d and n2_d_FAD
  diff = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < dim1+dim2; j++)
      diff = std::max(diff, fabs(n2_d(i,j).val() - n2_d_FAD(i,j).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;

}
/*----------------------------------------------------------------------*
 | End: FAD-Check for linearization of normal vector                    |
 *----------------------------------------------------------------------*/
#endif

#ifdef FADCHECKLINORTHOGONALITYCONDITION
/*----------------------------------------------------------------------*
 | FAD-Check for linearizations of orthogonality conditions             |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinOrthogonalityCondition(
    const int& fixed_par,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& rD,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi1,
    const LINALG::TMatrix<TYPEBTS, 3, 1>& x2_xi2,
    LINALG::TMatrix<TYPEBTS, 2, 2>& J_FAD,
    const LINALG::TMatrix<TYPEBTS, 2, 2>& J)
{

  // Get indices of not fixed parameters
  int par_i[2];
  int j = 0;
  for (int i = 0; i < 3; i++)
  {
    if (i != fixed_par)
    {
      par_i[j] = i;
      j++;
    }
  }

  // Compute norm of distance vector rD to scale the equations (this yields better conditioning)
  // NOTE: Even if automatic differentiation via FAD is applied, norm_rD has to be of type double since this
  // factor is needed for a pure scaling of the nonlinear orthogonality conditions and has not to be linearized!
  double norm_rD_scale = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(rD));

  // Evaluate f of orthogonality conditions
  LINALG::TMatrix<TYPEBTS, 2, 1> f(true);
  for (int i = 0; i < 3; i++)
  {
    f(0) += rD(i) * x2_xi1(i) / norm_rD_scale;
    f(1) += rD(i) * x2_xi2(i) / norm_rD_scale;
  }

  // Calculate linearization of orthogonality conditions (Jacobi-matrix J)
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      J_FAD(i,j) = f(i).dx(3*numnodes*numnodalvalues+3*numnodessol + par_i[j]);

  // Print linearization of f
  std::cout << "Lin: " << std::endl;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
      std::cout << J(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Print linearization of f calculated via FAD
  std::cout << "FAD: " << std::endl;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
      std::cout << J_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  return;
}
/*----------------------------------------------------------------------*
 | End: FAD-Check for linearizations of orthogonality conditions        |
 *----------------------------------------------------------------------*/
#endif

#ifdef FDCHECKSTIFFNESS
/*----------------------------------------------------------------------*
 | Differentation with finite difference for contact stiffness          |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FDCheckStiffness(
    const double& pp,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 1>& fc1,
    const LINALG::TMatrix<TYPEBTS, 3*numnodessol, 1>& fc2,
    LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1,
    LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Options for finite difference check
  const double h = 1.0e-9;
  const bool use_stiffc_FD = false;
  const bool output = true;

#ifdef GMSHDEBUG
  // Store class variables before adding h
  std::vector<gmshDebugPoint> gmshDebugPoints = gmshDebugPoints_;
  std::vector<std::pair<TYPEBTS, LINALG::TMatrix<TYPEBTS, 3, 1> > > normalsets = normalsets_;
#endif

  // Intialize matrices for contact stiffness
  LINALG::TMatrix<TYPEBTS, dim1, dim1+dim2> stiffc1_FD(true);
  LINALG::TMatrix<TYPEBTS, dim2, dim1+dim2> stiffc2_FD(true);

  // Loop over all coloums, size of displacement vector d = [d1, d2]^T
  for (int col = 0; col < dim1+dim2; col++)
  {
    // Initialize temporary varaibles
    TYPEBTS ele1pos_col = 0.0;
    TYPEBTS ele2pos_col = 0.0;

    // Store current element position before adding h and then add h to d(col)
    if (col < dim1)
    {
      // Beam element positions
      ele1pos_col = ele1pos_(col);
      ele1pos_(col) += h;
    }
    else
    {
      // Solid surface element positions
      ele2pos_col = ele2pos_(col-dim1);
      ele2pos_(col-dim1) += h;
    }


    // Find contact interval borders
    // -----------------------------------------------------------------

    // Vector for contact interval border parameter sets (xi1, xi2, eta and index of fixed paramater)
    std::vector<std::pair<LINALG::TMatrix<TYPEBTS, 3, 1>, int> > parsets;

    // Find contact interval borders
    GetContactIntervalBorders(parsets);

    // -----------------------------------------------------------------
    // End: Find contact interval borders


    // Compute contact forces and stiffness for all contact intervals
    // -----------------------------------------------------------------

    // Temporary vectors for contact forces at position with added h
    LINALG::TMatrix<TYPEBTS, dim1, 1> fc1_h(true);
    LINALG::TMatrix<TYPEBTS, dim2, 1> fc2_h(true);

    // Temporary matrices for contact stiffness at position with added h
    LINALG::TMatrix<TYPEBTS, dim1, dim1+dim2> stiffc1_h(true);
    LINALG::TMatrix<TYPEBTS, dim2, dim1+dim2> stiffc2_h(true);

    // Temporary matrices for contact stiffness calculated via FAD
    // TODO: Declare stiffc1_FAD and stiffc2_FAD only if needed (#ifdef AUTOMATICDIFFBTS)
    LINALG::TMatrix<TYPEBTS, dim1, dim1+dim2> stiffc1_FAD(true);
    LINALG::TMatrix<TYPEBTS, dim2, dim1+dim2> stiffc2_FAD(true);

#ifdef GMSHDEBUG
    // Clear class variables in each loop, because the information at positions with added h is not needed
    gmshDebugPoints_.clear();
    gmshDebugPoints_.resize(0);
    normalsets_.clear();
    normalsets_.resize(0);
#endif

    // Loop over all contact intervals, calculate and sum up contact forces (fc1, fc2) and stiffness (stiffc1, stiffc2)
    for (int i = 0; i < (int)parsets.size() - 1; i += 2)
    {
      bool doAssembleContactInterval = false;

      // Two parameter sets indicate the beginning (a) and the end (b) of the contact interval
      // TODO: Use stiffc1_FAD and stiffc2_FAD only if needed (#ifdef AUTOMATICDIFFBTS)
      EvaluateContactInterval(pp, parsets[i], parsets[i + 1], fc1_h, fc2_h, stiffc1_h, stiffc2_h,
          doAssembleContactInterval, stiffc1_FAD, stiffc2_FAD);
    }

    // -----------------------------------------------------------------
    // End: Compute contact forces and stiffness for all contact intervals


    // Restore element positions after the contact forces and stiffness at the position with added h are calculated
    if (col < dim1)
      ele1pos_(col) = ele1pos_col;
    else
      ele2pos_(col-dim1) = ele2pos_col;

    // Compute contact stiffness with finite difference
    for (int i = 0; i < dim1; i++)
    {
      // Beam element
      stiffc1_FD(i,col) = (fc1_h(i) - fc1(i)) / h;
    }
    for (int i = 0; i < dim2; i++)
    {
      // Solid element
      stiffc2_FD(i,col) = (fc2_h(i) - fc2(i)) / h;
    }
  }

  if (output)
  {
    // Print contact force fc1, stiffc1 and stiffc1_FD calculated via finite difference
    std::cout << "fc1: " << fc1 << std::endl;
    std::cout << "stiffc1: " << stiffc1 << std::endl;
    std::cout << "stiffc1_FD (h: " << h << "): " << stiffc1_FD << std::endl;

    // Set tolerance for comparison
    const double tolerance = 1.0e-12;

    // Calculate maximal difference between stiffc1 and stiffc1_FD
    double diff = 0.0;
    for (int i = 0; i < dim1; i++)
      for (int j = 0; j < dim1+dim2; j++)
        diff = std::max(diff, fabs(stiffc1(i,j) - stiffc1_FD(i,j)));
    std::cout << "Maximal difference: " << diff;
    if (diff > tolerance)
      std::cout << " > " << tolerance << " (tolerance)";
    std::cout << std::endl;


    // Print contact force fc2, stiffc2 and stiffc2_FD calculated via finite difference
    std::cout << "fc2: " << fc2 << std::endl;
    std::cout << "stiffc2: " << stiffc2 << std::endl;
    std::cout << "stiffc2_FD (h: " << h << "): " << stiffc2_FD << std::endl;

    // Calculate maximal difference between stiffc2 and stiffc2_FD
    diff = 0.0;
    for (int i = 0; i < dim2; i++)
      for (int j = 0; j < dim1+dim2; j++)
        diff = std::max(diff, fabs(stiffc2(i,j) - stiffc2_FD(i,j)));
    std::cout << "Maximal difference: " << diff;
    if (diff > tolerance)
      std::cout << " > " << tolerance << " (tolerance)";
    std::cout << std::endl;
  }

  // Assemble stiffc1_FD and stiffc2_FD instead of stiffc1 and stiffc2
  if (use_stiffc_FD)
  {
    stiffc1 = stiffc1_FD;
    stiffc2 = stiffc2_FD;
  }

#ifdef GMSHDEBUG
  // Restore class variables at position without added h
  gmshDebugPoints_ = gmshDebugPoints;
  normalsets_ = normalsets;
#endif

  return;
}
/*----------------------------------------------------------------------*
 | End: Differentation with finite difference for contact stiffness     |
 *----------------------------------------------------------------------*/
#endif

#ifdef FADCHECKSTIFFNESS
/*----------------------------------------------------------------------*
 | Check difference between stiffness and FAD stiffness                 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::FADCheckStiffness(
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1,
    const LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2,
    const LINALG::TMatrix<TYPEBTS, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc1_FAD,
    const LINALG::TMatrix<TYPEBTS, 3*numnodessol, 3*numnodes*numnodalvalues+3*numnodessol>& stiffc2_FAD)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodessol;

  // Print information
  std::cout << std::endl << "FAD-Check: Stiffmatrix:" << std::endl;

  // Print stiffc1 and stiffc1_FAD calculated via FAD
  std::cout << "stiffc1:" << std::endl;

  // Print stiffmatrix
  std::cout << "Lin:" << std::endl;
  for (int i = 0; i < dim1; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
        std::cout << stiffc1(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Print stiffmatrix calculated via FAD
  std::cout << "FAD:" << std::endl;
  for (int i = 0; i < dim1; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
        std::cout << stiffc1_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  const double tolerance = 1.0e-12;

  // Calculate maximal difference between stiffc1 and stiffc1_FAD
  double diff = 0.0;
  for (int i = 0; i < dim1; i++)
    for (int j = 0; j < dim1+dim2; j++)
      diff = std::max(diff, fabs(stiffc1(i,j).val() - stiffc1_FAD(i,j).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;


  // Print stiffc2 and stiffc2_FAD calculated via FAD
  std::cout << "stiffc2:" << std::endl;

  // Print stiffmatrix
  std::cout << "Lin:" << std::endl;
  for (int i = 0; i < dim2; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
        std::cout << stiffc2(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Print stiffmatrix calculated via FAD
  std::cout << "FAD:" << std::endl;
  for (int i = 0; i < dim2; i++)
  {
    for (int j = 0; j < dim1+dim2; j++)
        std::cout << stiffc2_FAD(i,j).val() << "  ";
    std::cout << std::endl;
  }

  // Calculate maximal difference between stiffc2 and stiffc2_FAD
  diff = 0.0;
  for (int i = 0; i < dim2; i++)
    for (int j = 0; j < dim1+dim2; j++)
      diff = std::max(diff, fabs(stiffc2(i,j).val() - stiffc2_FAD(i,j).val()));
  std::cout << "Maximal difference: " << diff;
  if (diff > tolerance)
    std::cout << " > " << tolerance << " (tolerance)";
  std::cout << std::endl;
}
/*----------------------------------------------------------------------*
 | End: Check difference between stiffness and FAD stiffness            |
 *----------------------------------------------------------------------*/
#endif

// Possible template cases: this is necessary for the compiler
//template class CONTACT::Beam3tosolidcontact<3,2,1>;
//template class CONTACT::Beam3tosolidcontact<3,3,1>;
//template class CONTACT::Beam3tosolidcontact<3,4,1>;
//template class CONTACT::Beam3tosolidcontact<3,5,1>;
//template class CONTACT::Beam3tosolidcontact<3,2,2>;

//template class CONTACT::Beam3tosolidcontact<4,2,1>;
//template class CONTACT::Beam3tosolidcontact<4,3,1>;
//template class CONTACT::Beam3tosolidcontact<4,4,1>;
//template class CONTACT::Beam3tosolidcontact<4,5,1>;
template class CONTACT::Beam3tosolidcontact<4,2,2>;     // Hermite beam element, quad4 surface element

//template class CONTACT::Beam3tosolidcontact<6,2,1>;
//template class CONTACT::Beam3tosolidcontact<6,3,1>;
//template class CONTACT::Beam3tosolidcontact<6,4,1>;
//template class CONTACT::Beam3tosolidcontact<6,5,1>;
//template class CONTACT::Beam3tosolidcontact<6,2,2>;

//template class CONTACT::Beam3tosolidcontact<8,2,1>;
//template class CONTACT::Beam3tosolidcontact<8,3,1>;
//template class CONTACT::Beam3tosolidcontact<8,4,1>;
//template class CONTACT::Beam3tosolidcontact<8,5,1>;
template class CONTACT::Beam3tosolidcontact<8,2,2>;     // Hermite beam element, quad8 surface element

//template class CONTACT::Beam3tosolidcontact<9,2,1>;
//template class CONTACT::Beam3tosolidcontact<9,3,1>;
//template class CONTACT::Beam3tosolidcontact<9,4,1>;
//template class CONTACT::Beam3tosolidcontact<9,5,1>;
//template class CONTACT::Beam3tosolidcontact<9,2,2>;
