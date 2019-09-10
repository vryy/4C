/*---------------------------------------------------------------------------*/
/*! \file
\brief Integrator class for the inequality constraint level-set approach
       (originally created by student Michael Hofer)

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------------*/

#include "../drt_contact_xcontact/xcontact_integrator.H"

#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_element.H"
#include "../drt_contact_xcontact/xcontact_debug_utils.H"
#include "../drt_contact_xcontact/xcontact_utils.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_contact/contact_paramsinterface.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>::Integrator(
    Teuchos::ParameterList& params, const Epetra_Comm& comm)
    : CONTACT::CoIntegrator(params, eletype, comm),
      is_const_normal_(
          DRT::INPUT::IntegralValue<bool>(params.sublist("XCONTACT", true), "CONST_CPP_NORMAL")),
      is_l2_var_jacobi_(true),
      is_h1_(
          DRT::INPUT::IntegralValue<bool>(params.sublist("XCONTACT", true), "H1_DUALITY_PAIRING"))
{
  // empty
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
Teuchos::RCP<const XCONTACT::GenericEvaluator<probdim, eletype, dim, numNodeEle>>
XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>::BuildEvaluator(
    const CONTACT::ParamsInterface& cparams) const
{
  Teuchos::RCP<const XCONTACT::GenericEvaluator<probdim, eletype>> xevaluator = Teuchos::null;

  switch (cparams.GetActionType())
  {
    case MORTAR::eval_force:
    {
      xevaluator = Teuchos::rcp(new const XCONTACT::ForceEvaluator<probdim, eletype>());
      break;
    }
    case MORTAR::eval_force_stiff:
    {
      xevaluator = Teuchos::rcp(new const XCONTACT::ForceTangentEvaluator<probdim, eletype>());
      break;
    }
    case MORTAR::eval_weighted_gap:
      xevaluator = Teuchos::rcp(new const XCONTACT::GapEvaluator<probdim, eletype>());
      break;
    default:
    {
      dserror("Unsupported action type in integrator call! (action = \"%s\")",
          MORTAR::ActionType2String(cparams.GetActionType()).c_str());
      exit(EXIT_FAILURE);
    }
  }

  return xevaluator;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
void XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>::IntegrateDerivEle2D(
    MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
    const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // ==========================================================================
  // Check integrator input for non-reasonable quantities
  // ==========================================================================

  // Check contact parameter interface pointer
  if (cparams_ptr.is_null())
    dserror("XContact integrator: Contact parameter interface pointer is undefined.");

  // Check problem dimension
  if (probdim != 2)
    dserror("XContact integrator: 2D integration method called for non-2D problem.");

  // Check contact element pairs for right definition of slave and master element
  for (unsigned me = 0; me < meles.size(); ++me)
  {
    if (!sele.IsSlave() or meles[me]->IsSlave())
      dserror(
          "XContact integrator: Integration called on wrong type of mortar "
          "element pair.");
  }

  //  // Check numerical integration type
  //  INPAR::MORTAR::IntType inttype =
  //      DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");
  //  if (inttype != INPAR::MORTAR::inttype_elements)
  //  {
  //    dserror("XContact integrator: Only element based integration implemented yet. "
  //        "Use INTTYPE Elements.");
  //  }

  // Check constant normal on master side
  if (not IsConstantNormal())
    dserror(
        "XContact integrator: Variation and linearization of non-constant "
        "closest point normal on master side not implemented yet.\n"
        "In case of constant normal on master side use CONST_CPP_NORMAL Yes.");

  // ==========================================================================
  // Integrate contact residuum and tangent matrix
  // ==========================================================================


  // Get the desired evaluator based on the given MORTAR::ActionType
  Teuchos::RCP<const XCONTACT::GenericEvaluator<probdim, eletype, dim, numNodeEle>> xevaluator =
      BuildEvaluator(*cparams_ptr);

  // --------------------------------------------------------------------------
  // Define slave element quantities
  // --------------------------------------------------------------------------

  // Get number of element nodes
  const int snnodes = numNodeEle;

  // Initialize shape functions and first parametric derivatives
  LINALG::Matrix<1, snnodes> sval(true);
  LINALG::Matrix<1, snnodes> lmval(true);
  LINALG::Matrix<dim, snnodes> sderiv(true);
  LINALG::Matrix<dim, snnodes> lmderiv(true);

  // Determine boundaries of element contact zone [sxia, sxib]
  const double sxi1a = -1.0;
  const double sxi1b = 1.0;

  // Compute constant element contact zone Jacobi determinant
  const double sjc = 0.5 * (sxi1b - sxi1a);

  // --------------------------------------------------------------------------
  // Integrate on slave element contact zone
  // --------------------------------------------------------------------------

  // Define Gauss point iterator
  typedef DRT::UTILS::GaussIntegration::const_iterator GI;

  // Get Gauss points
  const DRT::UTILS::GaussIntegration gaussPoints(eletype);

  // Loop over all Gauss points for integration
  for (GI gp = gaussPoints.begin(); gp != gaussPoints.end(); ++gp)
  {
    // Get current Gauss point coordinate and weight in unit interval [-1, 1]
    const double nu = gp.Point()[0];
    const double wgt = gp.Weight();

    // ------------------------------------------------------------------------
    // Evaluate surface geometry of slave element at Gauss point
    // ------------------------------------------------------------------------

    // Compute parametric coordinate
    LINALG::Matrix<dim, 1> sxi(true);
    sxi(0) = 0.5 * (sxi1a + sxi1b) + 0.5 * (sxi1b - sxi1a) * nu;

    // Get spatial coorinates of the slave element nodes
    LINALG::Matrix<probdim, numNodeEle> sxyze;
    SpatialCoordinates(sele, sxyze);
    // Evaluate coordinates of spatial point and shape functions
    LINALG::Matrix<probdim, 1> sx(true);
    SpatialPoint<probdim, eletype>(sxyze, sxi, sval, sx);

    // Evaluate covariant basis vectors and first parametric derivatives of shape functions
    LINALG::Matrix<dim, probdim> sg(true);
    CovariantBasisVectors<probdim, eletype>(sxyze, sxi, sderiv, sg);

    // Compute unit normal and length of non unit normal to compute element Jacobi determinant
    double sln = 0;
    LINALG::Matrix<probdim, 1> sn(true);
    Normal<probdim, eletype>(sele, sxi, sg, sn, sln);

    // Compute element Jacobi determinant and inverse
    const double sj = sln;
    const double sjinv = 1 / sj;

    // Evaluate contravariant basis vectors
    LINALG::Matrix<dim, probdim> sgc(true);
    InvMetrics<probdim, eletype>(sg, sj, sgc);

    // Evaluate surface gradient of spatial point
    LINALG::Matrix<probdim, probdim> sx_sx(true);
    GradSpatialPoint<probdim, eletype>(sg, sgc, sx_sx);

    // Get nodal Lagrange multipliers
    LINALG::Matrix<1, numNodeEle> lme;
    NodalLMValues(sele, lme);
    // Evaluate Lagrange multiplier in normal direction and shape functions
    // (defined on slave element)
    double lm = 0.0;
    LagrangeMultiplier<probdim, eletype>(lme, sxi, lmval, lm);

    /* Evaluate first parametric derivative of Lagrange multiplier in normal
     * direction and shape functions (defined on slave element) */
    LINALG::Matrix<1, dim> lm_sxi(true);
    LagrangeMultiplierParDeriv<probdim, eletype>(lme, sxi, lmderiv, lm_sxi);

    // ------------------------------------------------------------------------
    // Find corresponding master elements to current slave element
    // ------------------------------------------------------------------------

    // Loop over all master elements
    for (unsigned me = 0; me < meles.size(); ++me)
    {
      // ----------------------------------------------------------------------
      // Define master element quantities
      // ----------------------------------------------------------------------

      // Get number of element nodes
      const int mnnodes = numNodeEle;

      // Initialize shape functions and first parametric derivatives
      LINALG::Matrix<1, mnnodes> mval(true);
      LINALG::Matrix<dim, mnnodes> mderiv(true);

      // ----------------------------------------------------------------------
      // Solve closest point projection and evaluate surface geometry of
      // master element at projected Gauss point
      // ----------------------------------------------------------------------

      // Initialize parametric coordinate
      LINALG::Matrix<dim, 1> mxi(true);

      // Initialize coordinates of spatial point
      LINALG::Matrix<probdim, 1> mx(true);

      // Initialize covariant basis vectors
      LINALG::Matrix<dim, probdim> mg(true);

      // Initialize unit normal and length of non-unit normal
      LINALG::Matrix<probdim, 1> mn(true);

      /* Project Gauss point on slave element onto master element via closest
       * point projection */
      double binv = 0.0;
      xevaluator->ClosestPointProjection(sx, sg, *meles[me], mxi, mval, mderiv, mx, mg, binv);

      // Check if projected Gauss point is on master element
      if (mxi(0) >= -1.0 && mxi(0) <= 1.0)
      {
        // Evaluate unit normal and length of non-unit normal
        double mln = 0.0;
        LINALG::Matrix<probdim, 1> mn(true);
        Normal<probdim, eletype>(*meles[me], mxi, mg, mn, mln);

        // --------------------------------------------------------------------
        // Evaluate Gap in normal direction at Gauss point based on closest
        // point projection
        // --------------------------------------------------------------------
        double gN = 0.0;
        LINALG::Matrix<1, dim> gN_sxi(true);
        xevaluator->NormalGapCPP(sx, mx, mn, sg, gN, gN_sxi);

        // --------------------------------------------------------------------
        // Add contribution of Gauss point
        // --------------------------------------------------------------------
        xevaluator->Evaluate(sele, *meles[me], sval, mval, lmval, sderiv, mderiv, lmderiv, sxi, mxi,
            sx, mx, sx_sx, mg, mn, gN, gN_sxi, lm, lm_sxi, sj, sjinv, sjc, binv, wgt, *this);

        // Output Gauss point data
        const bool output = true;
        if (output)
        {
          DEBUG::OutputGaussPoint<probdim, eletype>(sx, mx, gN, gN_sxi);
        }
      }
    }
  }

  // Output element geometry
  DEBUG::OutputElement(sele);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
void XCONTACT::GenericEvaluator<probdim, eletype, dim, numNodeEle>::ClosestPointProjection(
    const LINALG::Matrix<probdim, 1>& sx, const LINALG::Matrix<dim, probdim>& sg,
    MORTAR::MortarElement& mele, LINALG::Matrix<dim, 1>& mxi, LINALG::Matrix<1, numNodeEle>& mval,
    LINALG::Matrix<dim, numNodeEle>& mderiv, LINALG::Matrix<probdim, 1>& mx,
    LINALG::Matrix<dim, probdim>& mg, double& binv) const
{
  // TODO: Extend 2D closest point projection to 3D
  if (probdim != 2)
  {
    dserror("Closest point projection only implemented for 2D yet.");
  }

  // ==========================================================================
  // Start point for local Newton iteration
  // ==========================================================================

  // Define start point for local Newton iteration in the element middle
  mxi.Scale(0.0);


  // ==========================================================================
  // Solve projection condition with local Newton iteration
  // ==========================================================================

  // Initialize projection condition
  double f = 1.0e12;

  // Get spatial coordinates of the master element nodes
  LINALG::Matrix<probdim, numNodeEle> mxyze;
  SpatialCoordinates(mele, mxyze);

  // Iterate until solution found or maximal number of iterations reached
  for (int iter = 0; iter < MORTARMAXITER; ++iter)
  {
    // ------------------------------------------------------------------------
    // Evaluate surface geometry of master element at current solution
    // ------------------------------------------------------------------------

    // Evaluate coordinates of spatial point and shape functions
    SpatialPoint<probdim, eletype>(mxyze, mxi, mval, mx);

    // Evaluate covariant basis vectors and first parametric derivatives of shape functions
    CovariantBasisVectors<probdim, eletype>(mxyze, mxi, mderiv, mg);

    // Compute parametric derivatives of covariant basis vectors
    LINALG::Matrix<dim, 1, LINALG::Matrix<probdim, 1>> mg_xi(true);
    ParDerivMetrics<probdim, eletype>(mxyze, mxi, mg_xi);


    // ------------------------------------------------------------------------
    // Update solution
    // ------------------------------------------------------------------------

    // Evaluate projection condition (residual)
    f = 0.0;
    for (unsigned i = 0; i < probdim; ++i)
    {
      f += mg(0, i) * (sx(i) - mx(i));
    }

    // Evaluate derivative of projection condition with respect to parameter coordinate of master
    // element (tangent matrix, Jacobi matrix)
    double f_mxi1 = 0.0;
    for (unsigned i = 0; i < probdim; ++i)
    {
      f_mxi1 += -mg(0, i) * mg(0, i) + (sx(i) - mx(i)) * mg_xi(0)(i, 0);
    }

    // Check for singularity of Jacobi matrix
    if (abs(f_mxi1) < 1.0e-12)
    {
      dserror("Closest point projection: Singular Jacobi matrix in local Newton iteration.");
    }

    // Invert tangent matrix
    binv = 1 / f_mxi1;

    // Check for fulfillment of projection condition
    if (abs(f) < MORTARCONVTOL)
    {
      const bool output = false;
      if (output)
      {
        std::cout << "f:    " << f << std::endl;
        std::cout << "iter: " << iter << std::endl;
        std::cout << "sx:   " << sx;
        std::cout << "mx:   " << mx;
        std::cout << "mg:   " << mg;
        std::cout << "------------------------------" << std::endl;
      }
      break;
    }

    // Update solution for projected parameter coordinate on master element
    mxi(0) += -binv * f;
  }

  // Check for unconverged local Newton iteration
  if (abs(f) > MORTARCONVTOL)
  {
    for (unsigned i = 0; i < dim; ++i)
    {
      mxi(i) = 1.0e12;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
void XCONTACT::GenericEvaluator<probdim, eletype, dim, numNodeEle>::NormalGapCPP(
    const LINALG::Matrix<probdim, 1>& sx, const LINALG::Matrix<probdim, 1>& mx,
    const LINALG::Matrix<probdim, 1>& mn, const LINALG::Matrix<dim, probdim>& sg, double& gN,
    LINALG::Matrix<1, dim>& gN_sxi) const
{
  // Compute gap in normal direction
  for (unsigned i = 0; i < probdim; ++i)
  {
    gN += mn(i) * (sx(i) - mx(i));
  }

  // Compute parametric derivative of gap
  for (unsigned a = 0; a < dim; ++a)
  {
    for (unsigned i = 0; i < probdim; ++i)
    {
      gN_sxi(a) += mn(i) * sg(a, i);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, bool wgap_only, unsigned dim,
    unsigned numNodeEle>
void XCONTACT::ForceEvaluator<probdim, eletype, wgap_only, dim, numNodeEle>::Evaluate(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<1, numNodeEle>& sval, const LINALG::Matrix<1, numNodeEle>& mval,
    const LINALG::Matrix<1, numNodeEle>& lmval, const LINALG::Matrix<dim, numNodeEle>& sderiv,
    const LINALG::Matrix<dim, numNodeEle>& mderiv, const LINALG::Matrix<dim, numNodeEle>& lmderiv,
    const LINALG::Matrix<dim, 1>& sxi, const LINALG::Matrix<dim, 1>& mxi,
    const LINALG::Matrix<probdim, 1>& sx, const LINALG::Matrix<probdim, 1>& mx,
    const LINALG::Matrix<probdim, probdim>& sx_sx, const LINALG::Matrix<dim, probdim>& mg,
    const LINALG::Matrix<probdim, 1>& mn, const double& gN, const LINALG::Matrix<1, dim>& gN_sxi,
    const double& lm, const LINALG::Matrix<1, dim>& lm_sxi, const double& sj, const double& sjinv,
    const double& sjc, const double& binv, const double& wgt,
    const XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>& integrator) const
{
  // Get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
  {
    dserror("Null pointer to element nodes.");
  }
  // ==========================================================================
  // Add contribution of Gauss point to constraint residual
  // ==========================================================================

  // Loop over all Lagrange multiplier nodes (slave element nodes)
  for (unsigned s = 0; s < integrator.SlNumNodeEle(); ++s)
  {
    CONTACT::CoNode* lmnode = dynamic_cast<CONTACT::CoNode*>(snodes[s]);

    double val = 0.0;

    // (L2-1) TODO
    val = lmval(s) * gN * sj * sjc * wgt;

    lmnode->AddgValue(val);

    // --- if we need only the weighted gap, we can skip the remaining part ---
    if (not wgap_only)
    {
      if (integrator.IsH1())
      {
        // (H1-1) TODO
        val += lmderiv(0, s) * gN_sxi(0) * sjinv * sjc * wgt;
      }

      lmnode->AddWcLm(val);
    }  // not wgap_only
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
    unsigned numNodeEle>
void XCONTACT::ForceTangentEvaluator<probdim, eletype, dim, numNodeEle>::Evaluate(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<1, numNodeEle>& sval, const LINALG::Matrix<1, numNodeEle>& mval,
    const LINALG::Matrix<1, numNodeEle>& lmval, const LINALG::Matrix<dim, numNodeEle>& sderiv,
    const LINALG::Matrix<dim, numNodeEle>& mderiv, const LINALG::Matrix<dim, numNodeEle>& lmderiv,
    const LINALG::Matrix<dim, 1>& sxi, const LINALG::Matrix<dim, 1>& mxi,
    const LINALG::Matrix<probdim, 1>& sx, const LINALG::Matrix<probdim, 1>& mx,
    const LINALG::Matrix<probdim, probdim>& sx_sx, const LINALG::Matrix<dim, probdim>& mg,
    const LINALG::Matrix<probdim, 1>& mn, const double& gN, const LINALG::Matrix<1, dim>& gN_sxi,
    const double& lm, const LINALG::Matrix<1, dim>& lm_sxi, const double& sj, const double& sjinv,
    const double& sjc, const double& binv, const double& wgt,
    const XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>& integrator) const
{
  // ==========================================================================
  // Add contribution of Gauss point to constraint residual
  // ==========================================================================
  XCONTACT::ForceEvaluator<probdim, eletype>::Evaluate(sele, mele, sval, mval, lmval, sderiv,
      mderiv, lmderiv, sxi, mxi, sx, mx, sx_sx, mg, mn, gN, gN_sxi, lm, lm_sxi, sj, sjinv, sjc,
      binv, wgt, integrator);

  // ==========================================================================
  // Define slave and master element quantities
  // ==========================================================================

  // Get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
  {
    dserror("Null pointer to element nodes.");
  }

  // Get master element nodes
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
  {
    dserror("Null pointer to element nodes.");
  }

  // Define map iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CIM;

  // Define number of slave element DOFs
  const int sndof = integrator.SlNumNodeEle() * probdim;

  // Define number of master element DOFs
  const int mndof = integrator.MaNumNodeEle() * probdim;


  // ==========================================================================
  // Evaluate variations/linearizations
  // ==========================================================================

  // --------------------------------------------------------------------------
  // Variation/linearization of gap in normal direction
  // --------------------------------------------------------------------------

  GEN::pairedvector<int, double> gN_u(sndof + mndof);

  for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int sDofId = snode->Dofs()[i];
      gN_u[sDofId] += sval(k) * mn(i);
    }
  }

  for (unsigned k = 0; k < integrator.MaNumNodeEle(); ++k)
  {
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int mDofId = mnode->Dofs()[i];
      gN_u[mDofId] += -mval(k) * mn(i);
    }
  }


  // --------------------------------------------------------------------------
  // Variation/linearization of slave element Jacobi determinant
  // --------------------------------------------------------------------------

  GEN::pairedvector<int, double> sj_su(sndof);
  double sxi_tmp[2] = {sxi(0), 0.0};
  sele.DerivJacobian(sxi_tmp, sj_su);  // TODO


  // --------------------------------------------------------------------------
  // Variation/linearization of projected parameter coordinate on
  // master element
  // --------------------------------------------------------------------------

  GEN::pairedvector<int, double> mxi_u(sndof + mndof);

  for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int sDofId = snode->Dofs()[i];
      mxi_u[sDofId] += -binv * mg(0, i) * sval(k);
    }
  }

  for (unsigned k = 0; k < integrator.MaNumNodeEle(); ++k)
  {
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int mDofId = mnode->Dofs()[i];
      mxi_u[mDofId] += -binv * (-mg(0, i) * mval(k) + (sx(i) - mx(i)) * mderiv(0, k));
    }
  }


  // --------------------------------------------------------------------------
  // Variation/linearization of parametric derivative of gap in normal
  // direction
  // --------------------------------------------------------------------------

  GEN::pairedvector<int, double> gN_sxi1_u(sndof + mndof);

  for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int sDofId = snode->Dofs()[i];
      gN_sxi1_u[sDofId] += sderiv(0, k) * mn(i);
    }
  }


  // --------------------------------------------------------------------------
  // Variation/linearization of inverse slave element Jacobi determinant
  // --------------------------------------------------------------------------

  GEN::pairedvector<int, double> sjinv_su(sndof);

  for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int sDofId = snode->Dofs()[i];
      sjinv_su[sDofId] += -sjinv * sjinv * sj_su[sDofId];
    }
  }


  // ==========================================================================
  // Evaluate linearizations of variations
  // ==========================================================================

  // Define identity matrix for later usage
  LINALG::Matrix<probdim, probdim> I(true);
  for (unsigned i = 0; i < probdim; ++i)
  {
    I(i, i) = 1;
  }

  // --------------------------------------------------------------------------
  // Linearization of variation of gap in normal direction
  // --------------------------------------------------------------------------

  // TODO: Nested pairedvector instead of nested map
  std::map<int, std::map<int, double>> gN_uu;
  std::map<int, std::map<int, double>> gN_sxi1_uu;

  for (unsigned k = 0; k < integrator.MaNumNodeEle(); ++k)
  {
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int mDofId = mnode->Dofs()[i];

      const double val = -mderiv(0, k) * mn(i);
      for (CI p = mxi_u.begin(); p != mxi_u.end(); ++p)
      {
        gN_uu[mDofId][p->first] += val * p->second;
      }
    }
  }

  // TODO: Test
  //  typedef std::map<int, std::map<int, double> >::const_iterator CIMM;
  //  for (CIMM p = gN_uu.begin(); p != gN_uu.end(); ++p)
  //  {
  //    for (CIM q = gN_uu[p->first].begin(); q != gN_uu[p->first].end(); ++q)
  //    {
  //      std::cout << p->first << "  " << q->first << "  " << gN_uu[p->first][q->first] <<
  //      std::endl;
  //    }
  //  }
  //  dserror("test");


  // --------------------------------------------------------------------------
  // Linearization of variation of slave element Jacobi determinant
  // --------------------------------------------------------------------------

  // TODO: Nested pairedvector instead of nested map
  std::map<int, std::map<int, double>> sj_susu;
  std::map<int, std::map<int, double>> sjinv_susu;

  for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
  {
    CONTACT::CoNode* snodek = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

    for (unsigned i = 0; i < probdim; ++i)
    {
      const int sDofIdk = snodek->Dofs()[i];

      for (unsigned l = 0; l < integrator.SlNumNodeEle(); ++l)
      {
        CONTACT::CoNode* snodel = dynamic_cast<CONTACT::CoNode*>(snodes[l]);

        for (unsigned j = 0; j < probdim; ++j)
        {
          const int sDofIdl = snodel->Dofs()[j];

          sj_susu[sDofIdk][sDofIdl] +=
              sderiv(0, k) * (I(i, j) - sx_sx(i, j)) * sjinv * sderiv(0, l);

          sjinv_susu[sDofIdk][sDofIdl] +=
              -sderiv(0, k) * (I(i, j) - 3 * sx_sx(i, j)) * sjinv * sjinv * sjinv * sderiv(0, l);
        }
      }
    }
  }

  // ==========================================================================
  // Add contribution of Gauss point to mortar matrices
  // ==========================================================================

  // Loop over all Lagrange multiplier nodes (slave element nodes)
  for (unsigned s = 0; s < integrator.SlNumNodeEle(); ++s)
  {
    CONTACT::CoNode* lmnode = dynamic_cast<CONTACT::CoNode*>(snodes[s]);

    // ------------------------------------------------------------------------
    // Add contribution of slave element
    // ------------------------------------------------------------------------

    // Loop over all slave element nodes
    for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
    {
      CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

      for (unsigned i = 0; i < probdim; ++i)
      {
        const int sDofId = snode->Dofs()[i];

        double val = 0.0;

        // (L2-1) Variation of gap in normal direction
        val = gN_u[sDofId] * lmval(s) * sj * sjc * wgt;

        if (integrator.IsL2VarJacobi())
        {
          // (L2-2) Variation of slave element Jacobi determinant
          val += sj_su[sDofId] * lmval(s) * gN * sjc * wgt;
        }

        if (integrator.IsH1())
        {
          // (H1-1) Variation of parametric derivative of gap in normal direction
          val += gN_sxi1_u[sDofId] * lmderiv(0, s) * sjinv * sjc * wgt;

          // (H1-2) Variation of inverse slave element Jacobi determinant
          val += sjinv_su[sDofId] * lmderiv(0, s) * gN_sxi(0) * sjc * wgt;
        }

        // Add contribution of slave element
        lmnode->AddWcSuLm(sDofId, val);
      }
    }

    // ------------------------------------------------------------------------
    // Add contribution of master element
    // ------------------------------------------------------------------------

    // Loop over all master element nodes
    for (unsigned k = 0; k < integrator.MaNumNodeEle(); ++k)
    {
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

      for (unsigned i = 0; i < probdim; ++i)
      {
        const int mDofId = mnode->Dofs()[i];

        // (L2-1) Variation of gap
        const double val = gN_u[mDofId] * lmval(s) * sj * sjc * wgt;

        // (L2-2) Variation of slave element Jacobi determinant
        // No contribution

        // (H1-1) Variation of parametric derivative of gap in normal direction
        // No contribution in case of constant normal

        // (H1-2) Variation of inverse slave element Jacobi determinant
        // No contribution

        // Add contribution of master element
        lmnode->AddWcMuLm(mDofId, val);
      }
    }
  }


  // ==========================================================================
  // Add contribution of Gauss point to structural contact tangent matrix
  // ==========================================================================

  // Loop over all Lagrange multiplier nodes (slave element nodes)
  for (unsigned s = 0; s < integrator.SlNumNodeEle(); ++s)
  {
    CONTACT::CoNode* lmnode = dynamic_cast<CONTACT::CoNode*>(snodes[s]);

    // ------------------------------------------------------------------------
    // Add contribution of slave element
    // ------------------------------------------------------------------------

    // Loop over all slave element nodes
    for (unsigned k = 0; k < integrator.SlNumNodeEle(); ++k)
    {
      CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[k]);

      for (unsigned i = 0; i < probdim; ++i)
      {
        const int sDofId = snode->Dofs()[i];
        std::map<int, double>& Wc_su_u = lmnode->CoData().GetWcSuU()[sDofId];

        double val = 0.0;

        // (L2-1) Linearization of variation of gap
        // No contribution in case of constant normal

        // (L2-2) Variation of gap and linearization of slave element Jacobi determinant
        val = gN_u[sDofId] * lmval(s) * sjc * wgt;
        for (CI p = sj_su.begin(); p != sj_su.end(); ++p)
        {
          Wc_su_u[p->first] += val * p->second;
        }

        if (integrator.IsL2VarJacobi())
        {
          // (L2-3) Variation of slave element Jacobi determinant and linearization of gap
          val = sj_su[sDofId] * lmval(s) * sjc * wgt;
          for (CI p = gN_u.begin(); p != gN_u.end(); ++p)
          {
            Wc_su_u[p->first] += val * p->second;
          }

          // (L2-4) Variation of linearization of slave element Jacobi determinant
          val = lmval(s) * gN * sjc * wgt;
          for (CIM p = sj_susu[sDofId].begin(); p != sj_susu[sDofId].end(); ++p)
          {
            Wc_su_u[p->first] += val * p->second;
          }
        }

        if (integrator.IsH1())
        {
          // (H1-1) Linearization of variation of parametric derivative of gap
          // No contribution in case of constant normal

          // (H1-2) Variation of parametric derivative of gap and linearization of inverse slave
          // element Jacobi determinant
          val = gN_sxi1_u[sDofId] * lmderiv(0, s) * sjc * wgt;
          for (CI p = sjinv_su.begin(); p != sjinv_su.end(); ++p)
          {
            Wc_su_u[p->first] += val * p->second;
          }

          // (H1-3) Variation of inverse slave element Jacobi determinant and linearization of
          // parametric derivative of gap
          val = sjinv_su[sDofId] * lmderiv(0, s) * sjc * wgt;
          for (CI p = gN_sxi1_u.begin(); p != gN_sxi1_u.end(); ++p)
          {
            Wc_su_u[p->first] += val * p->second;
          }

          // (H1-4) Variation of linearization of inverse slave element Jacobi determinant
          val = lmderiv(0, s) * gN_sxi(0) * sjc * wgt;
          for (CIM p = sjinv_susu[sDofId].begin(); p != sjinv_susu[sDofId].end(); ++p)
          {
            Wc_su_u[p->first] += val * p->second;
          }
        }
      }
    }

    // ------------------------------------------------------------------------
    // Add contribution of master element
    // ------------------------------------------------------------------------

    // Loop over all master element nodes
    for (unsigned k = 0; k < integrator.MaNumNodeEle(); ++k)
    {
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(mnodes[k]);

      for (unsigned i = 0; i < probdim; ++i)
      {
        const int mDofId = mnode->Dofs()[i];
        std::map<int, double>& Wc_mu_u = lmnode->CoData().GetWcMuU()[mDofId];

        double val = 0.0;

        // (L2-1) Linearization of variation of gap
        val = lmval(s) * sj * sjc * wgt;
        for (CIM p = gN_uu[mDofId].begin(); p != gN_uu[mDofId].end(); ++p)
        {
          Wc_mu_u[p->first] += val * p->second;
        }

        // (L2-2) Variation of gap and linearization of slave element Jacobi determinant
        val = gN_u[mDofId] * lmval(s) * sjc * wgt;
        for (CI p = sj_su.begin(); p != sj_su.end(); ++p)
        {
          Wc_mu_u[p->first] += val * p->second;
        }

        // (L2-3) Variation of slave element Jacobi determinant and linearization of gap
        // No contribution

        // (L2-4) Variation of linearization of slave element Jacobi determinant
        // No contribution

        if (integrator.IsH1())
        {
          // (H1-1) Linearization of variation of paramatric derivative of gap
          val = lmderiv(0, s) * sjinv * sjc * wgt;
          for (CIM p = gN_sxi1_uu[mDofId].begin(); p != gN_sxi1_uu[mDofId].end(); ++p)
          {
            Wc_mu_u[p->first] += val * p->second;
          }

          // (H1-2) Variation of parametric derivative of gap and linearization of inverse slave
          // element Jacobi determinant
          // No contribution in case of constant normal

          // (H1-3) Variation of inverse slave element Jacobi determinant and linearization of
          // parametric derivative of gap
          // No contribution

          // (H1-4) Variation of linearization of inverse slave element Jacobi determinant
          // No contribution
        }
      }
    }
  }
}

// Template class with possible problem dimension and slave element type
template class XCONTACT::Integrator<2, DRT::Element::line2>;

template class XCONTACT::GenericEvaluator<2, DRT::Element::line2>;

template class XCONTACT::ForceEvaluator<2, DRT::Element::line2>;
template class XCONTACT::ForceEvaluator<2, DRT::Element::line2, true>;

template class XCONTACT::ForceTangentEvaluator<2, DRT::Element::line2>;
