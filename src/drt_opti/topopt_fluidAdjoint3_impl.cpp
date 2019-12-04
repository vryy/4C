/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality of the element level of the fluid adjoint equations

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#include "topopt_fluidAdjoint3_impl.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/newtonianfluid.H"



#include "../linalg/linalg_utils_sparse_algebra_math.H"


//----------------------------------------------------------------------*
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidAdjoint3ImplInterface* DRT::ELEMENTS::FluidAdjoint3ImplInterface::Impl(
    DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::hex8:
    {
      return FluidAdjoint3Impl<DRT::Element::hex8>::Instance();
    }
    case DRT::Element::hex20:
    {
      return FluidAdjoint3Impl<DRT::Element::hex20>::Instance();
    }
    case DRT::Element::hex27:
    {
      return FluidAdjoint3Impl<DRT::Element::hex27>::Instance();
    }
    case DRT::Element::quad4:
    {
      return FluidAdjoint3Impl<DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return FluidAdjoint3Impl<DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return FluidAdjoint3Impl<DRT::Element::quad9>::Instance();
    }
    default:  // no 1D elements
    {
      dserror("Element shape %s not activated. Just do it.", DRT::DistypeToString(distype).c_str());
      break;
    }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype>* DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Instance(
    bool create)
{
  static FluidAdjoint3Impl<distype>* instance;
  if (create)
  {
    if (instance == NULL) instance = new FluidAdjoint3Impl<distype>();
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype>::FluidAdjoint3Impl()
    : xyze_(true),
      funct_(true),
      deriv_(true),
      deriv2_(true),
      xjm_(true),
      xji_(true),
      derxy_(true),
      derxy2_(true),
      visc_shp_(true),
      velint_(true),
      vderxy_(true),
      pres_(0.0),
      gradp_(true),
      visc_(true),
      fluidvelint_(true),
      fluidvelxy_(true),
      fluidpres_(0.0),
      fluidgradp_(true),
      fluidvisc_(true),
      fluidbodyforce_(true),
      bodyforce_(true),
      contforce_(0.0),
      vdiv_(0.0),
      conv1_(true),
      conv2_(true),
      velint_old_(true),
      vderxy_old_(true),
      pres_old_(0.0),
      gradp_old_(true),
      visc_old_(true),
      fluidvelint_old_(true),
      fluidvelxy_old_(true),
      bodyforce_old_(true),
      contforce_old_(0.0),
      vdiv_old_(0.0),
      conv1_old_(true),
      conv2_old_(true),
      fluidvelint_new_(true),
      fluidvelxy_new_(true),
      fluidgradp_new_(true),
      fluidvisc_new_(true),
      fluidbodyforce_new_(true),
      tau_(true),
      tau_old_(true),
      intpoints_(distype),
      xsi_(true),
      det_(0.0),
      fac_(0.0),
      reacoeff_(0.0),
      ele_(NULL),
      is_higher_order_ele_(false)
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fldAdPara_ = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elesysmat,
    Epetra_SerialDenseVector& elerhs)
{
  return Evaluate(ele, discretization, lm, params, mat, elesysmat, elerhs, intpoints_);
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elesysmat,
    Epetra_SerialDenseVector& elerhs, const DRT::UTILS::GaussIntegration& intpoints)
{
  ele_ = ele;

  // construct views
  LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_> elemat(elesysmat, true);
  LINALG::Matrix<(nsd_ + 1) * nen_, 1> elevec(elerhs, true);

  // ---------------------------------------------------------------------
  // get all general state vectors: fluid/adjoint velocity/pressure
  // velocity/pressure values are at time n/n+1
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_, nen_> eveln(true);
  LINALG::Matrix<nen_, 1> epren(true);
  ExtractValuesFromGlobalVector(discretization, lm, &eveln, &epren, "veln");

  LINALG::Matrix<nsd_, nen_> evelnp(true);
  LINALG::Matrix<nen_, 1> eprenp(true);
  ExtractValuesFromGlobalVector(discretization, lm, &evelnp, &eprenp, "velnp");

  LINALG::Matrix<nsd_, nen_> efluidveln(true);
  LINALG::Matrix<nen_, 1> efluidpren(true);
  ExtractValuesFromGlobalVector(discretization, lm, &efluidveln, &efluidpren, "fluidveln");

  LINALG::Matrix<nsd_, nen_> efluidvelnp(true);
  LINALG::Matrix<nen_, 1> efluidprenp(true);
  ExtractValuesFromGlobalVector(discretization, lm, &efluidvelnp, &efluidprenp, "fluidvelnp");

  LINALG::Matrix<nsd_, nen_> efluidvelnpp(true);
  LINALG::Matrix<nen_, 1> efluidprenpp(true);
  ExtractValuesFromGlobalVector(discretization, lm, &efluidvelnpp, &efluidprenpp, "fluidvelnpp");

  // evaluate nodal porosities
  LINALG::Matrix<nen_, 1> edens(true);
  if (params.get<INPAR::TOPOPT::DensityField>("dens_type") == INPAR::TOPOPT::dens_node_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    for (int nn = 0; nn < nen_; ++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();
      edens(nn, 0) = (*topopt_density)[lid];
    }
  }
  else if (params.get<INPAR::TOPOPT::DensityField>("dens_type") == INPAR::TOPOPT::dens_ele_based)
  {
    Teuchos::RCP<const Epetra_Vector> topopt_density =
        params.get<Teuchos::RCP<const Epetra_Vector>>("topopt_density");

    int lid = ele->LID();
    for (int nn = 0; nn < nen_; ++nn)  // set all values equal to hack a constant element porosity
                                       // on element level -> inefficient, but not relevant
      edens(nn, 0) = (*topopt_density)[lid];
  }
  else
    dserror("not implemented type of density function");

  LINALG::Matrix<nsd_, nen_> efluidbodyforcenp(true);
  LINALG::Matrix<nsd_, nen_> efluidbodyforcenpp(true);
  FluidBodyForce(fldAdPara_->Time(), efluidbodyforcenp);
  FluidBodyForce(fldAdPara_->Time() - fldAdPara_->Dt(), efluidbodyforcenpp);

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, xyze_);

  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;

  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (fldAdPara_->IsInconsistent() == true) is_higher_order_ele_ = false;
  // TODO deactivate this maybe?!

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(eveln, evelnp, epren, eprenp, efluidveln, efluidvelnp, efluidvelnpp, efluidpren,
      efluidprenp, efluidprenpp, efluidbodyforcenp, efluidbodyforcenpp, elemat, elevec, edens, mat,
      intpoints);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side        winklmaier 03/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Sysmat(const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nen_, 1>& epren,
    const LINALG::Matrix<nen_, 1>& eprenp, const LINALG::Matrix<nsd_, nen_>& efluidveln,
    const LINALG::Matrix<nsd_, nen_>& efluidvelnp, const LINALG::Matrix<nsd_, nen_>& efluidvelnpp,
    const LINALG::Matrix<nen_, 1>& efluidpren, const LINALG::Matrix<nen_, 1>& efluidprenp,
    const LINALG::Matrix<nen_, 1>& efluidprenpp,
    const LINALG::Matrix<nsd_, nen_>& efluidbodyforcenp,
    const LINALG::Matrix<nsd_, nen_>& efluidbodyforcenpp,
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const LINALG::Matrix<nen_, 1>& edens,
    Teuchos::RCP<const MAT::Material> material, const DRT::UTILS::GaussIntegration& intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  LINALG::Matrix<nen_ * nsd_, nen_ * nsd_> estif_w_v(true);
  LINALG::Matrix<nen_ * nsd_, nen_> estif_w_q(true);
  LINALG::Matrix<nen_, nen_ * nsd_> estif_r_v(true);
  LINALG::Matrix<nen_, nen_> estif_r_q(true);

  // definition of vectors
  LINALG::Matrix<nen_, 1> preforce(true);
  LINALG::Matrix<nsd_, nen_> velforce(true);

  // definition of velocity-based momentum residual vectors
  LINALG::Matrix<nsd_ * nsd_, nen_> lin_resM_Du(true);
  LINALG::Matrix<nsd_, 1> resM_Du(true);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter();

  // calculate subgrid viscosity and/or stabilization parameter at element center
  if (not fldAdPara_->EvalTauAtGP())
  {
    // get velocity at element center
    fluidvelint_.Multiply(efluidvelnp, funct_);
    fluidvelint_old_.Multiply(efluidveln, funct_);

    // calculate stabilization parameters at element center
    CalcStabParameter(fac_, fluidvelint_, tau_);
    CalcStabParameter(fac_, fluidvelint_old_, tau_old_);
  }

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad);

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives)
    //  2) fluid velocity (including derivatives)
    //  3) pressure (including derivatives)
    //  4) body-force vector
    //  5) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    velint_.Multiply(evelnp, funct_);
    velint_old_.Multiply(eveln, funct_);

    // get velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vderxy_.MultiplyNT(evelnp, derxy_);
    vderxy_old_.MultiplyNT(eveln, derxy_);

    // get fluid velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelint_.Multiply(efluidvelnp, funct_);
    fluidvelint_old_.Multiply(efluidveln, funct_);
    fluidvelint_new_.Multiply(efluidvelnpp, funct_);

    // get fluid velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelxy_.MultiplyNT(efluidvelnp, derxy_);
    fluidvelxy_old_.MultiplyNT(efluidveln, derxy_);
    fluidvelxy_new_.MultiplyNT(efluidvelnpp, derxy_);

    // get pressure at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    pres_ = funct_.Dot(eprenp);
    pres_old_ = funct_.Dot(epren);

    // get pressure gradient at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    gradp_.Multiply(derxy_, eprenp);
    gradp_old_.Multiply(derxy_, epren);

    // get pressure at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidpres_ = funct_.Dot(efluidprenp);

    // get pressure gradient at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidgradp_.Multiply(derxy_, efluidprenp);
    fluidgradp_new_.Multiply(derxy_, efluidprenpp);

    // get reaction coefficient due to porosity for topology optimization
    // !do this only at gauss point!
    {
      double densint = funct_.Dot(edens);
      const double* params = fldAdPara_->TopoptParams();

      if (params[2] > -1.0e-15)  // >=0 as it should be -> standard case
      {
        reacoeff_ =
            params[1] + (params[0] - params[1]) * densint * (1 + params[2]) / (densint + params[2]);
      }
      else  // special cases
      {
        INPAR::TOPOPT::OptiCase testcase = (INPAR::TOPOPT::OptiCase)round(-params[2]);

        switch (testcase)
        {
          case INPAR::TOPOPT::optitest_channel:
          {
            LINALG::Matrix<nsd_, 1> gp(true);
            gp.Multiply(xyze_, funct_);

            if (gp(1) < -0.1 || gp(1) > 0.1)  // wall area
              reacoeff_ = params[1];
            else
              reacoeff_ = params[0];
            break;
          }
          case INPAR::TOPOPT::optitest_channel_with_step:
          {
            LINALG::Matrix<nsd_, 1> gp(true);
            gp.Multiply(xyze_, funct_);

            if ((gp(0) > 1.5 and gp(0) < 1.9) and (gp(1) < 0.4))  // step -> wall
              reacoeff_ = params[1];
            else
              reacoeff_ = params[0];
            break;
          }
          case INPAR::TOPOPT::optitest_cornerflow:
          {
            LINALG::Matrix<nsd_, 1> gp(true);
            gp.Multiply(xyze_, funct_);

            if ((gp(0) > 0.875) or (gp(1) > 0.875) or (gp(0) < 0.625 and gp(1) < 0.625))  // wall
              reacoeff_ = params[1];
            else
              reacoeff_ = params[0];
            break;
          }
          case INPAR::TOPOPT::optitest_lin_poro:
          {
            double diff = params[1] - params[0];
            double pmax = params[1];

            reacoeff_ = -diff * densint + pmax;
            break;
          }
          case INPAR::TOPOPT::optitest_quad_poro:
          {
            double diff = params[1] - params[0];
            double pmax = params[1];

            double k = 0.1;

            reacoeff_ =
                (diff - k * pmax) * densint * densint + (-2 * diff + k * pmax) * densint + pmax;
            break;
          }
          case INPAR::TOPOPT::optitest_cub_poro:
          {
            double diff = params[1] - params[0];
            double pmax = params[1];

            double k1 = -50000.0;
            double k2 = 0.1;

            reacoeff_ = (2 * diff + k1 - k2 * pmax) * densint * densint * densint +
                        (-3 * diff - 2 * k1 + k2 * pmax) * densint * densint + k1 * densint + pmax;
            break;
          }
          default:
          {
            dserror("you should not be here with a testcase not handled above");
            break;
          }
        }
      }
    }

    // calculate fluid bodyforce on gausspoint
    fluidbodyforce_.Multiply(efluidbodyforcenp, funct_);
    fluidbodyforce_new_.Multiply(efluidbodyforcenpp, funct_);

    // calculate stabilization parameter at integration point
    if (fldAdPara_->EvalTauAtGP())
    {
      CalcStabParameter(fac_, fluidvelint_, tau_);
      CalcStabParameter(fac_, fluidvelint_old_, tau_old_);
    }

    BodyForce(efluidveln, efluidvelnp);
    ContForce();

    // get first convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv1_.Multiply(vderxy_, fluidvelint_);
    conv1_old_.Multiply(vderxy_old_, fluidvelint_old_);

    // get second convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv2_.MultiplyTN(fluidvelxy_, velint_);
    conv2_old_.MultiplyTN(fluidvelxy_old_, velint_old_);

    // get divergence at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vdiv_ = vdiv_old_ = 0.0;
    for (int idim = 0; idim < nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim, idim);
      vdiv_old_ += vderxy_old_(idim, idim);
    }

    CalcDivEps(eveln, evelnp, efluidvelnp, efluidvelnpp);

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = fldAdPara_->Timefac() * fac_;
    const double timefacfacrhs = fldAdPara_->TimefacRhs() * fac_;

    // TODO change these four factors maybe for genalpha
    const double timefacfacpre = fldAdPara_->Timefac() * fac_;
    const double timefacfacprerhs = fldAdPara_->TimefacRhs() * fac_;

    const double timefacfacdiv = fldAdPara_->Timefac() * fac_;
    const double timefacfacdivrhs = fldAdPara_->TimefacRhs() * fac_;

    /* ------------------------------------------------------------------------ *
     * standard terms                                                            *
     * ------------------------------------------------------------------------- */

    // 1) mass matrix + reactive term
    MassReactionGalPart(estif_w_v, velforce, timefacfac, timefacfacrhs);

    // 2) convection terms
    ConvectionGalPart(estif_w_v, velforce, timefacfac, timefacfacrhs);

    // 3) viscous terms
    ViscousGalPart(estif_w_v, velforce, timefacfac, timefacfacrhs);

    // 4) pressure term
    PressureGalPart(estif_w_q, velforce, timefacfacpre, timefacfacprerhs);

    // 5) continuity term
    ContinuityGalPart(estif_r_v, preforce, timefacfacdiv, timefacfacdivrhs);

    // 6) standard Galerkin bodyforce term on right-hand side
    BodyForceGalPart(velforce, timefacfac, timefacfacrhs);

    // 7) standard right-hand side term of continuity equation forces
    ContForceGalPart(preforce, timefacfacdiv, timefacfacdivrhs);
    /* ------------------------------------------------------------------------ *
     * standard terms done                                                       *
     * ------------------------------------------------------------------------- */



    /* ------------------------------------------------------------------------ *
     * stabilization part                                                        *
     * ------------------------------------------------------------------------- */

    if ((fldAdPara_->PSPG()) or (fldAdPara_->SUPG()))
    {
      if (fldAdPara_->AdjointType() == INPAR::TOPOPT::discrete_adjoint)
      {
        // prework for supg/psgp - stabilization: evaluate strong residual

        /* order of the derivatives in GalMomResnU is:
         * from 1 to nsd:       col-dim = x, row-dim = 1-nsd
         * from nsd+1 to 2*nsd: col-dim = y, row-dim = 1-nsd
         * and so on. so the outer loop is the column dimension
         * and the inner loop the row dimension */
        LINALG::Matrix<nsd_ * nsd_, nen_> GalMomTestStat(true);

        DiscreteGalMom(GalMomTestStat, timefacfac, timefacfacrhs, timefacfacpre, timefacfacprerhs);

        // 8) PSPG term
        if (fldAdPara_->PSPG())
        {
          DiscretePSPG(estif_w_q, estif_r_q, velforce, preforce, GalMomTestStat, timefacfac,
              timefacfacrhs, timefacfacpre, timefacfacprerhs);
        }

        // 9) SUPG term
        if (fldAdPara_->SUPG())
        {
          DiscreteSUPG(estif_w_v, estif_r_v, velforce, preforce, GalMomTestStat, timefacfac,
              timefacfacrhs, timefacfacpre, timefacfacprerhs);
        }
      }
      else if (fldAdPara_->AdjointType() == INPAR::TOPOPT::cont_adjoint)
      {
        // prework for supg/psgp - stabilization: evaluate strong residual

        /* order of the derivatives in GalMomResnU is:
         * from 1 to nsd:       col-dim = x, row-dim = 1-nsd
         * from nsd+1 to 2*nsd: col-dim = y, row-dim = 1-nsd
         * and so on. so the outer loop is the column dimension
         * and the inner loop the row dimension */
        LINALG::Matrix<nsd_ * nsd_, nen_> GalMomResnU(true);

        // strong residual of momentum equation of last iteration, scaled with fac*dt/rho
        LINALG::Matrix<nsd_, 1> StrongResMomScaled(true);

        MomRes(GalMomResnU, StrongResMomScaled, timefacfac, timefacfacrhs, timefacfacpre,
            timefacfacprerhs);

        // 8) PSPG term
        if (fldAdPara_->PSPG())
        {
          PSPG(estif_r_v, estif_r_q, preforce, GalMomResnU, StrongResMomScaled, timefacfac,
              timefacfacrhs, timefacfacpre, timefacfacprerhs);
        }

        // 9) SUPG term
        if (fldAdPara_->SUPG())
        {
          SUPG(estif_w_v, estif_w_q, velforce, GalMomResnU, StrongResMomScaled, timefacfac,
              timefacfacrhs, timefacfacpre, timefacfacprerhs);
        }
      }
      else
        dserror("not implemented type of adjoint approach");
    }

    // 10) continuity stabilization
    if (fldAdPara_->CStab())
    {
      if (fldAdPara_->AdjointType() == INPAR::TOPOPT::discrete_adjoint)
      {
        DiscreteContStab(estif_w_v, velforce, timefacfacdiv, timefacfacdivrhs);
      }
      else if (fldAdPara_->AdjointType() == INPAR::TOPOPT::cont_adjoint)
      {
        ContStab(estif_w_v, velforce, timefacfacdiv, timefacfacdivrhs);
      }
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    eforce(numdofpernode_ * vi + nsd_) += preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      eforce(numdofpernode_ * vi + idim) += velforce(idim, vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fuippp = numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_ * vi + nsd_;

      estif(numdof_vi_p_nsd, fuippp) += estif_r_q(vi, ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_ * vi;
        const int nsd_vi = nsd_ * vi;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif(numdof_vi + idim, numdof_ui_jdim) += estif_w_v(nsd_vi + idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int nsd_vi = nsd_ * vi;
      const int numdof_vi = numdofpernode_ * vi;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif(numdof_vi + idim, numdof_ui_nsd) += estif_w_q(nsd_vi + idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
        estif(numdofpernode_ * vi + nsd_, numdof_ui_jdim) += estif_r_v(vi, nsd_ui_jdim);
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at element center     winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtEleCenter()
{
  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_stab(
      DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim = 0; idim < nsd_; idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }
  const double wquad = intpoints_stab.IP().qwgt[0];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_, funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_, deriv2_);
  }


  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_, xyze_);
  det_ = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", ele_->Id(), det_);

  // compute integration factor
  fac_ = wquad * det_;

  // compute global first derivates
  derxy_.Multiply(xji_, deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype, nen_>(xjm_, derxy_, deriv2_, xyze_, derxy2_);
  }
  else
    derxy2_.Clear();

  return;
}



/*-------------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point    winklmaier 02/12 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::GaussIntegration::iterator& iquad  // actual integration point
)
{
  // coordinates of the current integration point
  const double* gpcoord = iquad.Point();
  for (int idim = 0; idim < nsd_; idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_, funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);
  derxy2_.Clear();
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_, deriv2_);
  }

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  xjm_.MultiplyNT(deriv_, xyze_);
  det_ = xji_.Invert(xjm_);

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", ele_->Id(), det_);

  // compute integration factor
  fac_ = iquad.Weight() * det_;

  // compute global first derivates
  derxy_.Multiply(xji_, deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype, nen_>(xjm_, derxy_, deriv2_, xyze_, derxy2_);
  }
  else
    derxy2_.Clear();

  return;
}



/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter             winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcStabParameter(
    const double vol, const LINALG::Matrix<nsd_, 1>& fluidvel, LINALG::Matrix<3, 1>& tau)
{
  //---------------------------------------------------------------------
  // preliminary definition of values which will already be computed for
  // tau_M and later be used for tau_C again by some of the subsequent
  // stabilization parameter definitions
  //---------------------------------------------------------------------
  double traceG = 0.0;
  double Gnormu = 0.0;
  double Gvisc = 0.0;

  double strle = 0.0;
  double hk = 0.0;
  double fluidvel_norm = 0.0;
  double re12 = 0.0;
  double c3 = 0.0;

  // material parameters
  const double dens = fldAdPara_->Density();
  const double visc = fldAdPara_->Viscosity();

  //---------------------------------------------------------------------
  // first step: computation of tau_M with the following options
  // (both with or without inclusion of dt-part):
  // A) definition according to Taylor et al. (1998)
  //    -> see also Gravemeier and Wall (2010) for version for
  //       variable-density flow at low Mach number
  // B) combined definition according to Franca and Valentin (2000) as
  //    well as Barrenechea and Valentin (2002)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // C) definition according to Shakib (1989) / Shakib and Hughes (1991)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // D) definition according to Codina (1998)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // E) definition according to Franca et al. (2005) as well as Badia
  //    and Codina (2010)
  //    -> only for Darcy or Darcy-Stokes/Brinkman flow, hence only
  //       tau_Mp for this definition
  //---------------------------------------------------------------------
  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // computation depending on which parameter definition is used
  switch (fldAdPara_->TauType())
  {
    case INPAR::FLUID::tau_taylor_hughes_zarins:
    case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
    {
      /*

      literature:
      1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
         of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
         (1998) 155-196.
      2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
         multigrid method for large-eddy simulation of turbulent variable-
         density flow at low Mach number, J. Comput. Phys. 229 (2010)
         6047-6070.
         -> version for variable-density low-Mach-number flow as implemented
            here, which corresponds to version for incompressible flow as
            given in the previous publications when density is constant

                                                                             1
                       +-                                               -+ - -
                       |        2                                        |   2
                       | c_1*rho                                  2      |
            tau  = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
               M       |     2                                           |
                       |   dt                                            |
                       +-                                               -+

            with the constants and covariant metric tensor defined as follows:

            C   = 1.0 (not explicitly defined here),
            c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
            c_2 = 1.0 (not explicitly defined here),
            c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)

                    +-           -+   +-           -+   +-           -+
                    |             |   |             |   |             |
                    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
              G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
               ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                    |    i     j  |   |    i     j  |   |    i     j  |
                    +-           -+   +-           -+   +-           -+

                    +----
                     \
            G : G =   +   G   * G
                     /     ij    ij
                    +----
                     i,j
                               +----
                               \
            rho*u*G*rho*u  =   +   rho*u * G  *rho*u
                               /        i   ij      j
                              +----
                                i,j
      */

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = reacoeff_;
      if (fldAdPara_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins or
          fldAdPara_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
          fldAdPara_->TauType() == INPAR::FLUID::tau_taylor_hughes_zarins_scaled)
        sigma_tot += 1.0 / fldAdPara_->Dt();

      // definition of constants as described above
      const double c1 = 4.0;
      c3 = 12.0 / mk;

      // computation of various values derived from covariant metric tensor
      // (trace of covariant metric tensor required for computation of tau_C below)
      double G;
      double normG = 0.0;
      const double dens_sqr = dens * dens;
      for (int nn = 0; nn < nsd_; ++nn)
      {
        const double dens_sqr_velint_nn = dens_sqr * fluidvel(nn);
        for (int mm = 0; mm < nsd_; ++mm)
        {
          traceG += xji_(nn, mm) * xji_(nn, mm);
        }
        for (int rr = 0; rr < nsd_; ++rr)
        {
          G = xji_(nn, 0) * xji_(rr, 0);
          for (int mm = 1; mm < nsd_; ++mm)
          {
            G += xji_(nn, mm) * xji_(rr, mm);
          }
          normG += G * G;
          Gnormu += dens_sqr_velint_nn * G * fluidvel(rr);
        }
      }

      // compute viscous part
      Gvisc = c3 * visc * visc * normG;

      // computation of stabilization parameters tau_Mu and tau_Mp
      // -> identical for the present definitions
      tau(0) = 1.0 / (sqrt(c1 * dens_sqr * DSQR(sigma_tot) + Gnormu + Gvisc));
      tau(1) = tau(0);
      break;
    }

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
    {
      /*

      literature:
      1) L.P. Franca, F. Valentin, On an improved unusual stabilized
         finite element method for the advective-reactive-diffusive
         equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.
      2) G.R. Barrenechea, F. Valentin, An unusual stabilized finite
         element method for a generalized Stokes problem, Numer. Math.
         92 (2002) 652-677.


                    xi1,xi2 ^
                            |      /
                            |     /
                            |    /
                          1 +---+
                            |
                            |
                            |
                            +--------------> re1,re2
                                1

      */
      // get velocity norm
      fluidvel_norm = fluidvel.Norm2();

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      const double sigma_tot = 1.0 / fldAdPara_->Timefac() + reacoeff_;

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel, fluidvel_norm, strle, hk);

      // various parameter computations for case with dt:
      // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
      const double re01 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(strle));
      const double re11 = 4.0 * visc / (mk * dens * sigma_tot * DSQR(hk));

      // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
      const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
      re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

      // respective "switching" parameters
      const double xi01 = std::max(re01, 1.0);
      const double xi11 = std::max(re11, 1.0);
      const double xi02 = std::max(re02, 1.0);
      const double xi12 = std::max(re12, 1.0);

      tau(0) = DSQR(strle) / (DSQR(strle) * dens * sigma_tot * xi01 + (4.0 * visc / mk) * xi02);
      tau(1) = DSQR(hk) / (DSQR(hk) * dens * sigma_tot * xi11 + (4.0 * visc / mk) * xi12);
      break;
    }

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
    {
      /*

       stabilization parameter as above without inclusion of dt-part

      */
      // get velocity norm
      fluidvel_norm = fluidvel.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel, fluidvel_norm, strle, hk);

      // various parameter computations for case without dt:
      // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
      double re01 = 4.0 * visc / (mk * dens * reacoeff_ * DSQR(strle));
      double re11 = 4.0 * visc / (mk * dens * reacoeff_ * DSQR(hk));

      // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
      const double re02 = mk * dens * fluidvel_norm * strle / (2.0 * visc);
      re12 = mk * dens * fluidvel_norm * hk / (2.0 * visc);

      // respective "switching" parameters
      const double xi01 = std::max(re01, 1.0);
      const double xi11 = std::max(re11, 1.0);
      const double xi02 = std::max(re02, 1.0);
      const double xi12 = std::max(re12, 1.0);

      tau(0) = DSQR(strle) / (DSQR(strle) * dens * reacoeff_ * xi01 + (4.0 * visc / mk) * xi02);
      tau(1) = DSQR(hk) / (DSQR(hk) * dens * reacoeff_ * xi11 + (4.0 * visc / mk) * xi12);
      break;
    }

    case INPAR::FLUID::tau_shakib_hughes_codina:
    case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
    {
      /*

      literature:
      1) F. Shakib, Finite element analysis of the compressible Euler and
         Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
         Stanford University, Stanford, CA, USA, 1989.
      2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
         computational fluid dynamics: IX. Fourier analysis of space-time
         Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
         Engrg. 87 (1991) 35-58.
      3) R. Codina, Stabilized finite element approximation of transient
         incompressible flows using orthogonal subscales, Comput. Methods
         Appl. Mech. Engrg. 191 (2002) 4295-4321.

         constants defined as in Shakib (1989) / Shakib and Hughes (1991),
         merely slightly different with respect to c_3:

         c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
         c_2 = 4.0,
         c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

         Codina (2002) proposed present version without dt and explicit
         definition of constants
         (condition for constants as defined here: c_2 <= sqrt(c_3)).

      */
      // get velocity norm
      fluidvel_norm = fluidvel.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel, fluidvel_norm, strle, hk);

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = reacoeff_;
      if (fldAdPara_->TauType() == INPAR::FLUID::tau_shakib_hughes_codina)
        sigma_tot += 1.0 / fldAdPara_->Dt();

      // definition of constants as described above
      const double c1 = 4.0;
      const double c2 = 4.0;
      c3 = 4.0 / (mk * mk);
      // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

      tau(0) = 1.0 / (sqrt(c1 * DSQR(dens) * DSQR(sigma_tot) +
                           c2 * DSQR(dens) * DSQR(fluidvel_norm) / DSQR(strle) +
                           c3 * DSQR(visc) / (DSQR(strle) * DSQR(strle))));
      tau(1) = 1.0 / (sqrt(c1 * DSQR(dens) * DSQR(sigma_tot) +
                           c2 * DSQR(dens) * DSQR(fluidvel_norm) / DSQR(hk) +
                           c3 * DSQR(visc) / (DSQR(hk) * DSQR(hk))));
      break;
    }
    case INPAR::FLUID::tau_codina:
    case INPAR::FLUID::tau_codina_wo_dt:
    {
      /*

        literature:
           R. Codina, Comparison of some finite element methods for solving
           the diffusion-convection-reaction equation, Comput. Methods
           Appl. Mech. Engrg. 156 (1998) 185-210.

           constants:
           c_1 = 1.0 (for version with dt), 0.0 (for version without dt),
           c_2 = 2.0,
           c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

           Codina (1998) proposed present version without dt.

      */
      // get velocity norm
      fluidvel_norm = fluidvel.Norm2();

      // calculate characteristic element length
      CalcCharEleLength(vol, fluidvel, fluidvel_norm, strle, hk);

      // total reaction coefficient sigma_tot: sum of "artificial" reaction
      // due to time factor and reaction coefficient (reaction coefficient
      // ensured to remain zero in GetMaterialParams for non-reactive material)
      double sigma_tot = reacoeff_;
      if (fldAdPara_->TauType() == INPAR::FLUID::tau_codina) sigma_tot += 1.0 / fldAdPara_->Dt();

      // definition of constants as described above
      const double c1 = 1.0;
      const double c2 = 2.0;
      c3 = 4.0 / mk;

      tau(0) = 1.0 / (sqrt(c1 * dens * sigma_tot + c2 * dens * fluidvel_norm / strle +
                           c3 * visc / DSQR(strle)));
      tau(1) =
          1.0 /
          (sqrt(c1 * dens * sigma_tot + c2 * dens * fluidvel_norm / hk + c3 * visc / DSQR(hk)));
      break;
    }
    default:
    {
      dserror("unknown definition for tau_M\n %i  ", fldAdPara_->TauType());
      break;
    }
  }  // end switch (fldAdPara_->whichtau_)


  //---------------------------------------------------------------------
  // second step: computation of tau_C with the following options:
  // A) definition according to Taylor et al. (1998)
  // B) definition according to Whiting (1999)/Whiting and Jansen (2001)
  // C) scaled version of definition according to Taylor et al. (1998)
  // D) definition according to Wall (1999)
  // E) definition according to Codina (2002)
  // F) definition according to Badia and Codina (2010)
  //    (only for Darcy or Darcy-Stokes/Brinkman flow)
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (fldAdPara_->TauType())
  {
    case INPAR::FLUID::tau_taylor_hughes_zarins:
    case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
    {
      /*

      literature:
         C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
         of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
         (1998) 155-196.

                                                1/2
                             (c_2*rho*u*G*rho*u)
                      tau  = -------------------
                         C       trace (G)


         -> see respective definitions for computation of tau_M above

      */

      tau(2) = sqrt(Gnormu) / traceG;
    }
    break;

    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
    case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
    {
      /*

      literature:
      1) C.H. Whiting, Stabilized finite element methods for fluid dynamics
         using a hierarchical basis, PhD thesis, Rensselaer Polytechnic
         Institute, Troy, NY, USA, 1999.
      2) C.H. Whiting, K.E. Jansen, A stabilized finite element method for
         the incompressible Navier-Stokes equations using a hierarchical
         basis, Int. J. Numer. Meth. Fluids 35 (2001) 93-116.

                                    1.0
                      tau  = ------------------
                         C    tau  * trace (G)
                                 M

         -> see respective definitions for computation of tau_M above

      */

      tau(2) = 1.0 / (tau_(0) * traceG);
    }
    break;

    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
    case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
    {
      /*

        Caution: This is an experimental version of a stabilization
                 parameter definition which scales the definition
                 for tau_C by Taylor et al. (1998) in a similar
                 way as proposed below by Franca and Frey (1992)
                 and Wall (1999) by appropriately defining an
                 element Reynolds number based on the covariant
                 metric tensor.

                    /                        1/2    \
                    |  /                    \       |                       1/2
                    | |  c_2*rho*u*G*rho*u  |       |    (c_2*rho*u*G*rho*u)
        tau  =  MIN | | ------------------- | | 1.0 | *  -------------------
           C        | |          2          |       |         trace (G)
                    | \    c_3*mu *G:G      /       |
                    \                               /
                      |                     |
                      -----------------------
                      element Reynolds number
                        based on covariant
                          metric tensor

         -> see respective definitions for computation of tau_M above

      */

      // element Reynolds number based on covariant metric tensor
      const double reG = std::sqrt(Gnormu / Gvisc);

      // "switching" parameter
      const double xi_tau_c = std::min(reG, 1.0);

      tau(2) = xi_tau_c * sqrt(Gnormu) / traceG;
    }
    break;

    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
    case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
    {
      /*

      literature:
      1) L.P. Franca, S.L. Frey, Stabilized finite element methods:
         II. The incompressible Navier-Stokes equations, Comput. Methods
         Appl. Mech. Engrg. 99 (1992) 209-293.
      2) W.A. Wall, Fluid-Struktur-Interaktion mit stabilisierten Finiten
         Elementen, Dissertation, Universitaet Stuttgart, 1999.

                   xi_tau_c ^
                            |
                          1 |   +-----------
                            |  /
                            | /
                            |/
                            +--------------> re12
                                1

         -> see respective definitions for computation of tau_M above

      */

      // "switching" parameter
      const double xi_tau_c = std::min(re12, 1.0);

      tau(2) = 0.5 * dens * fluidvel_norm * hk * xi_tau_c;
    }
    break;

    case INPAR::FLUID::tau_shakib_hughes_codina:
    case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
    case INPAR::FLUID::tau_codina:
    case INPAR::FLUID::tau_codina_wo_dt:
    {
      /*

      literature:
         R. Codina, Stabilized finite element approximations of transient
         incompressible flows using orthogonal subscales, Comput. Methods
         Appl. Mech. Engrg. 191 (2002) 4295-4321.

         -> see respective definitions for computation of tau_M above

      */

      tau(2) = DSQR(hk) / (sqrt(c3) * tau(1));
    }
    break;
    default:
    {
      dserror("unknown definition for tau_C\n %i  ", fldAdPara_->TauType());
      break;
    }
  }  // end switch (fldAdPara_->whichtau_)

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of characteristic element length       winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcCharEleLength(const double vol,
    const LINALG::Matrix<nsd_, 1> fluidvel, const double fluidvel_norm, double& strle,
    double& hk) const
{
  // cast dimension to a double varibale -> pow()
  const double dim = double(nsd_);

  //! direction of flow (normed velocity vector)
  LINALG::Matrix<nsd_, 1> fluidvelino;

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mu
  //---------------------------------------------------------------------
  // a) streamlength due to Tezduyar et al. (1992) -> default
  // normed velocity vector
  if (fluidvel_norm >= 1e-6)
    fluidvelino.Update(1.0 / fluidvel_norm, fluidvel);
  else
  {
    fluidvelino.Clear();
    fluidvelino(0, 0) = 1.0;
  }

  LINALG::Matrix<nen_, 1> tmp;
  tmp.MultiplyTN(derxy_, fluidvelino);
  const double val = tmp.Norm1();
  strle = 2.0 / val;

  // b) volume-equivalent diameter (warning: 3-D formula!)
  // strle = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // c) cubic/square root of element volume/area
  // strle = std::pow(vol,1/dim);

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mp
  //---------------------------------------------------------------------
  // a) volume-equivalent diameter -> default for 3-D computations
  if (nsd_ == 3) hk = std::pow((6. * vol / M_PI), (1.0 / 3.0)) / sqrt(3.0);

  // b) square root of element area -> default for 2-D computations,
  // may also alternatively be used for 3-D computations
  else if (nsd_ == 2)
    hk = std::pow(vol, 1 / dim);
  // check for potential 1-D computations
  else
    dserror("element length calculation not implemented for 1-D computation!");

  return;
}



/*---------------------------------------------------------------------------------*
 | compute bodyforce of momentum equation                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForce(
    const LINALG::Matrix<nsd_, nen_>& efluidveln, const LINALG::Matrix<nsd_, nen_>& efluidvelnp)
{
  bodyforce_.Clear();
  bodyforce_old_.Clear();

  if (fldAdPara_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    const double dissipation =
        fldAdPara_->ObjDissipationFac();  // zero if no dissipation part in obj-fcn
    /* ------------------------------------------------------------------------ *
     * 1) evaluate bodyforce at new time step                                   *
     * ------------------------------------------------------------------------ */

    // dissipation term due to reaction
    if (fldAdPara_->ObjDissipationTerm() == INPAR::TOPOPT::obj_diss_yes)
      bodyforce_.Update(-2 * dissipation * reacoeff_, fluidvelint_);

    // dissipation term due to viscosity
    if (is_higher_order_ele_)  // TODO check this
    {
      LINALG::Matrix<nsd_, numderiv2_> fluidvelxy2(true);
      fluidvelxy2.MultiplyNT(efluidvelnp, derxy2_);

      LINALG::Matrix<nsd_, 1> laplaceU(true);
      for (int idim = 0; idim < nsd_; ++idim)
      {
        for (int jdim = 0; jdim < nsd_; ++jdim) laplaceU(idim) += fluidvelxy2(idim, jdim);
      }

      bodyforce_.Update(dissipation * fldAdPara_->Viscosity(), laplaceU, 1.0);
    }


    /* ------------------------------------------------------------------------ *
     * 2) evaluate bodyforce at old time step in instationary case              *
     * ------------------------------------------------------------------------ */
    if (not fldAdPara_->IsStationary())
    {
      if (fldAdPara_->ObjDissipationTerm() == INPAR::TOPOPT::obj_diss_yes)
        bodyforce_old_.Update(-2 * dissipation * reacoeff_, fluidvelint_old_);

      // dissipation term due to viscosity
      if (is_higher_order_ele_)  // TODO check this
      {
        LINALG::Matrix<nsd_, numderiv2_> fluidvelxy2_old(true);
        fluidvelxy2_old.MultiplyNT(efluidveln, derxy2_);

        LINALG::Matrix<nsd_, 1> laplaceU_old(true);
        for (int idim = 0; idim < nsd_; ++idim)
        {
          for (int jdim = 0; jdim < nsd_; ++jdim) laplaceU_old(idim) += fluidvelxy2_old(idim, jdim);
        }

        bodyforce_old_.Update(dissipation * fldAdPara_->Viscosity(), laplaceU_old, 1.0);
      }
    }
  }
  else  // special cases
  {
    // get global coordinates of gauss point
    double x = 0.0;
    double y = 0.0;
    LINALG::Matrix<nsd_, 1> coords(true);
    coords.Multiply(xyze_, funct_);
    x = coords(0);
    y = coords(1);  // z-component currently not required in tests


    switch (fldAdPara_->TestCase())
    {
      case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
      {
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
      {
        bodyforce_(0) = bodyforce_old_(0) = 4 * y;
        bodyforce_(1) = bodyforce_old_(1) = 24994 * x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
      {
        bodyforce_(0) = bodyforce_old_(0) = 24999 * x * x - 4 * x * y;
        bodyforce_(1) = bodyforce_old_(1) = -74989 * x * x + 50002 * y * y + 12 * x * y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
      {
        bodyforce_(0) = bodyforce_old_(0) = 24998 * x * x + 24994 * x * y - 4 * y * y;
        bodyforce_(1) = bodyforce_old_(1) = -74978 * x * x + 50004 * y * y + 28 * x * y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_varying_theta:
      {
        double t = fldAdPara_->Time();
        bodyforce_(0) = 5 * x;
        bodyforce_(1) = -5 * y;

        t += fldAdPara_->Dt();  // old time = t + dt
        bodyforce_old_(0) = 5 * x;
        bodyforce_old_(1) = -5 * y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
      {
        double t = fldAdPara_->Time();
        bodyforce_(0) = -10 * x * y + 5 * x * x - 4 * y * y * t + 9 * x * y * t;
        bodyforce_(1) = -8 * y * y * t + x * x + 24 * x * y + 18 * y * y * t * t + 4 * x * y * t;

        t += fldAdPara_->Dt();  // old time = t + dt
        bodyforce_old_(0) = -10 * x * y + 5 * x * x - 4 * y * y * t + 9 * x * y * t;
        bodyforce_old_(1) =
            -8 * y * y * t + x * x + 24 * x * y + 18 * y * y * t * t + 4 * x * y * t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
      {
        double t = fldAdPara_->Time();
        bodyforce_(0) = -3 * x * y - 3 * y * t * t - 6 * y * y * t - 6 * y + 8 * x * y * t;
        bodyforce_(1) = 9 * x + x * t + 6 * x * y * t;

        t += fldAdPara_->Dt();  // old time = t + dt
        bodyforce_old_(0) = -3 * x * y - 3 * y * t * t - 6 * y * y * t - 6 * y + 8 * x * y * t;
        bodyforce_old_(1) = 9 * x + x * t + 6 * x * y * t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_primal:
        break;
      default:
      {
        dserror("no dirichlet condition implemented for special test case");
        break;
      }
    }
  }
}



/*------------------------------------------------------------------------------*
 |  compute fluid body force at element nodes (public)         winklmaier 04/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::FluidBodyForce(
    const double time, LINALG::Matrix<nsd_, nen_>& efluidbodyforce)
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  if (nsd_ == 3)
    DRT::UTILS::FindElementConditions(ele_, "VolumeNeumann", myneumcond);
  else if (nsd_ == 2)
    DRT::UTILS::FindElementConditions(ele_, "SurfaceNeumann", myneumcond);
  else
    dserror("Body force for 1D problem not yet implemented!");

  if (myneumcond.size() > 1) dserror("More than one Neumann condition on one node!");

  if (myneumcond.size() == 1)
  {
    const std::string* condtype = myneumcond[0]->Get<std::string>("type");

    // get values and switches from the condition
    const std::vector<int>* onoff = myneumcond[0]->Get<std::vector<int>>("onoff");
    const std::vector<double>* val = myneumcond[0]->Get<std::vector<double>>("val");
    const std::vector<int>* functions = myneumcond[0]->Get<std::vector<int>>("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the ebofoaf array
    for (int isd = 0; isd < nsd_; isd++)
    {
      // get factor given by spatial function
      if (functions)
        functnum = (*functions)[isd];
      else
        functnum = -1;

      double num = (*onoff)[isd] * (*val)[isd];

      for (int jnode = 0; jnode < nen_; ++jnode)
      {
        if (functnum > 0)
        {
          // evaluate function at the position of the current node
          // ------------------------------------------------------
          // comment: this introduces an additional error compared to an
          // evaluation at the integration point. However, we need a node
          // based element bodyforce vector for prescribed pressure gradients
          // in some fancy turbulance stuff.
          functionfac = DRT::Problem::Instance()
                            ->Funct(functnum - 1)
                            .Evaluate(isd, (ele_->Nodes()[jnode])->X(), time);
        }
        else
          functionfac = 1.0;

        // get usual body force
        if (*condtype == "neum_dead" or *condtype == "neum_live")
          efluidbodyforce(isd, jnode) = num * functionfac;
        // get prescribed pressure gradient
        else if (*condtype == "neum_pgrad")
          dserror("not implemented for adjoints in several parts of the code");
        else
          dserror("Unknown Neumann condition");
      }
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute continuity force                                       winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContForce()
{
  contforce_ = 0.0;
  contforce_old_ = 0.0;

  if (fldAdPara_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    ;  // currently no domain pressure entries in objective -> no entry here
  }
  else
  {
    // get global coordinates of gauss point
    double x = 0.0;
    double y = 0.0;
    LINALG::Matrix<nsd_, 1> coords(true);
    coords.Multiply(xyze_, funct_);
    x = coords(0);
    y = coords(1);  // z-component currently not required in tests

    switch (fldAdPara_->TestCase())
    {
      case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
      {
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
      {
        contforce_ = contforce_old_ = 12.0;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
      {
        contforce_ = contforce_old_ = 2 * x + 4 * y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
      {
        contforce_ = contforce_old_ = 2 * x + 5 * y;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_varying_theta:
      {
        contforce_ = contforce_old_ = 0;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
      {
        double t = fldAdPara_->Time();
        contforce_ = 2 * x + y * t + 4 * y * t * t;

        t += fldAdPara_->Dt();  // old time = t + dt
        contforce_old_ = 2 * x + y * t + 4 * y * t * t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
      {
        double t = fldAdPara_->Time();
        contforce_ = y * t + 1 + t;

        t += fldAdPara_->Dt();  // old time = t + dt
        contforce_old_ = y * t + 1 + t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_primal:
        break;
      default:
      {
        dserror("no dirichlet condition implemented for special test case");
        break;
      }
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute mass and reaction terms of galerkin part               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MassReactionGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfac, const double& timefacfacrhs) const
{
  /* inertia (contribution to mass matrix) if not is_stationary */
  /*
            /              \
           |                |
           |    rho*Dv , w  |
           |                |
            \              /
  */
  /*  reaction */
  /*
            /                \
           |                  |
           |    sigma*Dv , w  |
           |                  |
            \                /
  */
  double massreacfac = 0.0;  // factor summing up coefficients of reactive term and mass-matrix

  if (fldAdPara_->IsStationary())
    massreacfac = reacoeff_ * timefacfac;
  else
    massreacfac = fldAdPara_->Density() * fac_ +
                  reacoeff_ * timefacfac;  // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;

    const double uifunct = massreacfac * funct_(ui);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = nsd_ * vi;

      const double value = funct_(vi) * uifunct;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_w_v(fvi + idim, fui + idim) += value;
      }  // end for (idim)
    }    // vi
  }      // ui

  // rhs at new time step
  LINALG::Matrix<nsd_, 1> scaled_vel(true);
  scaled_vel.Update(massreacfac, velint_);
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      velforce(jdim, vi) -= funct_(vi) * scaled_vel(jdim);
    }
  }

  // rhs at old time step
  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
  {
    double massreacfacrhs =
        -fldAdPara_->Density() * fac_ +
        reacoeff_ * timefacfacrhs;  // fac -> mass matrix // reac*timefacfac/dens -> reactive
    scaled_vel.Update(massreacfacrhs, velint_old_);

    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        velforce(jdim, vi) -= funct_(vi) * scaled_vel(jdim);
      }
    }
  }
  return;
}



/*---------------------------------------------------------------------------------*
 | compute (two!) convection terms of galerkin part               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ConvectionGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfac, const double& timefacfacrhs) const
{
  /* convection, part 1 */
  /*
            /                          \
           |       /  n         \       |
     rho * | Dv , |  u   o nabla | w    |
           |       \            /   (i) |
            \                          /
  */

  double value = 0.0;  // helper

  for (int vi = 0; vi < nen_; ++vi)
  {
    const int fvi = nsd_ * vi;

    value = 0.0;
    for (int dim = 0; dim < nsd_; ++dim) value += derxy_(dim, vi) * fluidvelint_(dim);
    value *= fldAdPara_->Density() * timefacfac;

    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui = nsd_ * ui;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_w_v(fvi + idim, fui + idim) += value * funct_(ui);
      }  // end for (idim)
    }    // ui
  }      // vi

  // rhs at new time step
  for (int vi = 0; vi < nen_; ++vi)
  {
    value = 0.0;
    for (int dim = 0; dim < nsd_; ++dim) value += derxy_(dim, vi) * fluidvelint_(dim);
    value *= fldAdPara_->Density() * timefacfac;

    for (int jdim = 0; jdim < nsd_; ++jdim) velforce(jdim, vi) -= value * velint_(jdim);
  }  // vi

  // rhs at old time step
  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      value = 0.0;
      for (int dim = 0; dim < nsd_; ++dim) value += derxy_(dim, vi) * fluidvelint_(dim);
      value *= fldAdPara_->Density() * timefacfacrhs;

      for (int jdim = 0; jdim < nsd_; ++jdim) velforce(jdim, vi) -= value * velint_old_(jdim);
    }  // vi
  }

  /*  convection, part 2 */
  /*
            /                           \
           |  /               \   n      |
     rho * | |  Dv o nabla     | u  , w  |
           |  \           (i) /          |
            \                           /
  */
  for (int ui = 0; ui < nen_; ++ui)
  {
    value = fldAdPara_->Density() * timefacfac * funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fui = nsd_ * ui + idim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          estif_w_v(fvi + jdim, fui) += funct_(vi) * fluidvelxy_(idim, jdim) * value;
        }
      }
    }
  }

  // rhs at new time step
  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    value = 0.0;  // product of fluid and adjoint velocity
    for (int dim = 0; dim < nsd_; ++dim) value += fluidvelxy_(dim, jdim) * velint_(dim);
    value *= fldAdPara_->Density() * timefacfac;

    for (int vi = 0; vi < nen_; ++vi) velforce(jdim, vi) -= funct_(vi) * value;
  }

  // rhs at new and old time step
  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      value = 0.0;  // product of fluid and adjoint velocity
      for (int dim = 0; dim < nsd_; ++dim) value += fluidvelxy_(dim, jdim) * velint_old_(dim);
      value *= fldAdPara_->Density() * timefacfacrhs;

      for (int vi = 0; vi < nen_; ++vi) velforce(jdim, vi) -= funct_(vi) * value;
    }
  }
}



/*---------------------------------------------------------------------------------*
 | compute viscous terms of galerkin part                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ViscousGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfac, const double& timefacfacrhs) const
{
  /* viscosity term: overall */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Dv | , eps | w |  |
                  |       \  /         \ /   |
                   \                        /
  */

  /* split in two parts: */

  /* viscosity term, part 1 */
  /*
                   /                                    \
                  |         /      \           /     \   |
              mu  |  nabla | Dv     | , nabla | w     |  |
                  |         \  (i) /           \ (i) /   |
                   \                                    /
  */

  /* viscosity term, part 2 */
  /*
                   /                                  \
                  |           /  \           /     \   |
              mu  |  nabla   | Dv | , nabla | w     |  |
                  |       (i) \  /           \ (i) /   |
                   \                                  /
  */

  double value = 0.0;  // helper
  const double viscdenstimefac = timefacfac * fldAdPara_->Viscosity();

  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = nsd_ * vi;

      value = 0.0;

      for (int dim = 0; dim < nsd_; ++dim) value += derxy_(dim, vi) * derxy_(dim, ui);

      value *= viscdenstimefac;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_w_v(fvi + idim, fui + idim) += value;
      }  // end for (idim)
    }    // vi
  }      // ui

  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;  // shp fcn of u known, derivative not

    for (int jdim = 0; jdim < nsd_; ++jdim)  // derivative of u known by jdim
    {
      value = viscdenstimefac * derxy_(jdim, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;  // shp fcn of v known, derivative not

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_w_v(fvi + jdim, fui + idim) += derxy_(idim, vi) * value;
        }
      }  // vi
    }
  }  // ui

  // viscosity at new and old time step (if instationary)
  LINALG::Matrix<nsd_, nsd_> viscstress(true);

  double viscdenstimefacrhs = timefacfacrhs * fldAdPara_->Viscosity();  // for instationary problems
  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      viscstress(idim, jdim) += viscdenstimefac * (vderxy_(jdim, idim) + vderxy_(idim, jdim));

      if (not fldAdPara_->IsStationary())
        viscstress(idim, jdim) +=
            viscdenstimefacrhs * (vderxy_old_(jdim, idim) + vderxy_old_(idim, jdim));
    }
  }


  // computation of right-hand-side viscosity term
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
        velforce(idim, vi) -= viscstress(idim, jdim) * derxy_(jdim, vi);
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute pressure term of galerkin part                         winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PressureGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_>& estif_w_q, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfacpre, const double& timefacfacprerhs) const
{
  /* pressure term */
  /*
       /                \
      |                  |
    + | Dq  , nabla o w  |
      |                  |
       \                /
  */

  double value = 0.0;  // helper

  for (int ui = 0; ui < nen_; ++ui)
  {
    value = timefacfacpre * funct_(ui);

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * nsd_;

      for (int jdim = 0; jdim < nsd_; ++jdim) estif_w_q(fvi + jdim, ui) += derxy_(jdim, vi) * value;
    }  // vi
  }    // ui

  // rhs at new time step
  value = timefacfacpre * pres_;
  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
    value += timefacfacprerhs * pres_old_;
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)  // derivative of u known by jdim
      velforce(jdim, vi) -= derxy_(jdim, vi) * value;
  }  // vi

  return;
}



/*---------------------------------------------------------------------------------*
 | compute continuity term of galerkin part                       winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContinuityGalPart(
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_r_v, LINALG::Matrix<nen_, 1>& preforce,
    const double& timefacfacdiv, const double& timefacfacdivrhs) const
{
  double value = 0.0;  // helper
  /* continuity term */
  /*
       /                \
      |                  |
    - | nabla o Dv  , r  |
      |                  |
       \                /
  */
  for (int vi = 0; vi < nen_; ++vi)
  {
    value = timefacfacdiv * funct_(vi);

    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * nsd_;

      for (int idim = 0; idim < nsd_; ++idim) estif_r_v(vi, fui + idim) -= value * derxy_(idim, ui);
    }  // vi
  }    // ui

  // rhs at new and old time step
  value = timefacfacdiv * vdiv_;
  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
    value += timefacfacdivrhs * vdiv_old_;

  for (int vi = 0; vi < nen_; ++vi) preforce(vi) += funct_(vi) * value;


  return;
}


/*---------------------------------------------------------------------------------*
 | compute body force term of galerkin part                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForceGalPart(
    LINALG::Matrix<nsd_, nen_>& velforce, const double& timefacfac,
    const double& timefacfacrhs) const
{
  double value = 0.0;

  if (fldAdPara_->TestCase() == INPAR::TOPOPT::adjointtest_no)
  {
    double objfac = fac_ * fldAdPara_->Dt();
    if (abs(fldAdPara_->Time()) < 1.0e-10)
      objfac *= 1.0 - fldAdPara_->ThetaObj();  // time 0
    else if (fldAdPara_->IsInitInstatStep())
      objfac *= fldAdPara_->ThetaObj();  // time T
    else
      objfac *= 1.0;

    const double dissipation =
        fldAdPara_->ObjDissipationFac();  // zero if no dissipation part in obj-fcn

    /*
     *  d   /             \                    /         \
     * --- |   reac*u*u   | (w)   =   2*reac |   u , w   |
     *  du  \             /                    \         /
     */
    if (fldAdPara_->ObjDissipationTerm() == INPAR::TOPOPT::obj_diss_yes)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        value = 2 * objfac * dissipation * reacoeff_ * fluidvelint_(idim);

        for (int vi = 0; vi < nen_; ++vi)
        {
          velforce(idim, vi) -= value * funct_(vi);
        }
      }  // end for(idim)
    }


    /*
     *  d   /                    \                  /                    \
     * --- |  2*mu*eps(u)*eps(u)  | (w)   =   4*mu |   eps(u) , nabla w   |
     *  du  \                    /                  \                    /
     */
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        value = 2 * objfac * fldAdPara_->Viscosity() *
                (fluidvelxy_(idim, jdim) + fluidvelxy_(jdim, idim));

        for (int vi = 0; vi < nen_; ++vi)
        {
          velforce(jdim, vi) -= derxy_(idim, vi) * value;
        }
      }
    }  // end for(idim)
  }
  else  // special cases -> no partial integration
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      value = timefacfac * bodyforce_(idim);

      if (not fldAdPara_->IsStationary()) value += timefacfacrhs * bodyforce_old_(idim);

      for (int vi = 0; vi < nen_; ++vi)
      {
        velforce(idim, vi) += value * funct_(vi);
      }
    }  // end for(idim)
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute continuity force term of galerkin part                 winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContForceGalPart(LINALG::Matrix<nen_, 1>& preforce,
    const double& timefacfacdiv, const double& timefacfacdivrhs) const
{
  double value = timefacfacdiv * contforce_;

  if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
    value += timefacfacdivrhs * contforce_old_;

  for (int vi = 0; vi < nen_; ++vi)
  {
    preforce(vi, 0) -= value * funct_(vi);
  }


  return;
}



/*---------------------------------------------------------------------------------*
 | compute momentum residuum                                      winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscreteGalMom(
    LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomTestStat, const double& timefacfac,
    const double& timefacfacrhs, const double& timefacfacpre, const double& timefacfacprerhs) const
{
  /*
//      Left hand side terms of Galerkin part for PSPG/SUPG with Dv
//
//    instationary + reactive         /
//                                   |
//      (rho + alpha) w +  dt*Theta  |
//                                   |
//                                    \
//
//                convective term 1             convective term 2 TODO hier anders!!!
//     VZ!!!     /  n             \           /               n    \
//      + rho * |  u o nabla w     | + rho * |   w o nabla   u      |
//               \            (i) /           \               (i)  /
//
//                     viscous term            \
//     VZ!!!    /                         \     |
//      + 2 mu |  nabla o epsilon   ( w )  |    |
//              \                (i)      /     |
//                                             /
  */
  int idim_nsd_p_idim[nsd_];

  for (int idim = 0; idim < nsd_; ++idim)
  {
    idim_nsd_p_idim[idim] = idim * nsd_ + idim;
  }

  const double timefacfac_densaf = timefacfac * fldAdPara_->Density();

  for (int ui = 0; ui < nen_; ++ui)
  {
    double value = 0.0;

    for (int dim = 0; dim < nsd_; ++dim) value += fluidvelint_(dim) * derxy_(dim, ui);

    value *= fldAdPara_->Density() * timefacfac;

    for (int idim = 0; idim < nsd_; ++idim)
    {
      GalMomTestStat(idim_nsd_p_idim[idim], ui) += value;
    }
  }


  // dr_j   d    /    du_j \          du_j         dN_B
  // ----= ---- | u_i*----  | = N_B * ---- + u_i * ---- * d_jk
  // du_k  du_k  \    dx_i /          dx_k         dx_i

  for (int ui = 0; ui < nen_; ++ui)
  {
    const double temp = timefacfac_densaf * funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int idim_nsd = idim * nsd_;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        GalMomTestStat(idim_nsd + jdim, ui) += temp * fluidvelxy_(idim, jdim);
      }
    }
  }


  const double fac_reac = timefacfac * reacoeff_;

  for (int ui = 0; ui < nen_; ++ui)
  {
    const double v = fac_reac * funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      GalMomTestStat(idim_nsd_p_idim[idim], ui) += v;
    }
  }

  // viscous
  if (is_higher_order_ele_)
  {
    const double v = -2.0 * fldAdPara_->Viscosity() * timefacfac;
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int nsd_idim = nsd_ * jdim;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim_p_jdim = nsd_idim + idim;

        for (int ui = 0; ui < nen_; ++ui)
        {
          GalMomTestStat(nsd_idim_p_jdim, ui) += v * visc_shp_(nsd_idim_p_jdim, ui);
        }
      }
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute momentum residuum                                      winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MomRes(
    LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomResnU, LINALG::Matrix<nsd_, 1>& StrongResMomScaled,
    const double& timefacfac, const double& timefacfacrhs, const double& timefacfacpre,
    const double& timefacfacprerhs) const
{
  /*
      Left hand side terms of Galerkin part for PSPG/SUPG with Dv

    instationary + reactive         /
                                   |
      (rho + alpha) Du +  dt*Theta |
                                   |
                                    \

                convective term 1             convective term 2
               /  n              \           /               n \
      - rho * |  u o nabla Dv     | + rho * |  Dv o nabla   u   |
               \             (i) /           \           (i)   /

                     viscous term            \
              /                         \     |
      - 2 mu |  nabla o epsilon   ( Dv ) |    |
              \                (i)      /     |
                                             /
  */

  GalMomResnU.Clear();

  // mass matrix + reaction
  double massreacfac = 0.0;  // factor summing up coefficients of reactive term and mass-matrix
  if (fldAdPara_->IsStationary())
    massreacfac = reacoeff_ * timefacfac;
  else
    massreacfac = fldAdPara_->Density() * fac_ +
                  reacoeff_ * timefacfac;  // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui = 0; ui < nen_; ++ui)
  {
    const double uifunct = massreacfac * funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim) GalMomResnU(idim * nsd_ + idim, ui) += uifunct;
  }  // ui

  //   convection
  for (int ui = 0; ui < nen_; ++ui)
  {
    double value = 0.0;

    for (int dim = 0; dim < nsd_; ++dim) value += fluidvelint_(dim) * derxy_(dim, ui);

    value *= fldAdPara_->Density() * timefacfac;

    for (int idim = 0; idim < nsd_; ++idim) GalMomResnU(idim * nsd_ + idim, ui) -= value;
  }  // ui

  for (int ui = 0; ui < nen_; ++ui)
  {
    const double uifunct = fldAdPara_->Density() * timefacfac * funct_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        GalMomResnU(jdim + idim * nsd_, ui) += fluidvelxy_(jdim, idim) * uifunct;
      }
    }
  }

  // viscous
  if (is_higher_order_ele_)
  {
    // add viscous part
    GalMomResnU.Update(-2.0 * timefacfac * fldAdPara_->Viscosity(), visc_shp_, 1.0);
  }


  // residuum of momentum equation in strong form
  if (not fldAdPara_->IsStationary())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      StrongResMomScaled(idim) =
          fldAdPara_->Density() * fac_ *
              (velint_(idim) - velint_old_(idim))  // mass term last iteration
          + timefacfac *  // velocity part of last iteration (at t^n) coming
                (fldAdPara_->Density() * (-conv1_(idim) + conv2_(idim)) -
                    2 * fldAdPara_->Viscosity() * visc_(idim) + reacoeff_ * velint_(idim) -
                    bodyforce_(idim)) -
          timefacfacpre * gradp_(idim)  // pressure part of last iteration (at t^n)
          + timefacfacrhs *             // last time step (= t^n+1) coming
                (fldAdPara_->Density() * (-conv1_old_(idim) + conv2_old_(idim)) -
                    2 * fldAdPara_->Viscosity() * visc_old_(idim) + reacoeff_ * velint_old_(idim) -
                    bodyforce_old_(idim)) -
          timefacfacprerhs * gradp_old_(idim);
    }
  }
  else
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      StrongResMomScaled(idim) =
          timefacfac * (fldAdPara_->Density() * (-conv1_(idim) + conv2_(idim)) -
                           2 * fldAdPara_->Viscosity() * visc_(idim) + reacoeff_ * velint_(idim) -
                           gradp_(idim) - bodyforce_(idim));
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute divergence of epsilon of v                             winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcDivEps(const LINALG::Matrix<nsd_, nen_>& eveln,
    const LINALG::Matrix<nsd_, nen_>& evelnp, const LINALG::Matrix<nsd_, nen_>& efluidvelnp,
    const LINALG::Matrix<nsd_, nen_>& efluidvelnpp)
{
  /*--- viscous term: div(epsilon(u)) --------------------------------*/
  /*   /                                                \
       |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
     1 |                                                |
     - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
     2 |                                                |
       |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
       \                                                /

       with N_x .. x-line of N
       N_y .. y-line of N                                             */

  visc_.Clear();
  visc_old_.Clear();
  fluidvisc_.Clear();
  fluidvisc_new_.Clear();

  if (nsd_ == 3)
  {
    for (int inode = 0; inode < nen_; ++inode)
    {
      double sum = (derxy2_(0, inode) + derxy2_(1, inode) + derxy2_(2, inode));
      visc_shp_(0, inode) = 0.5 * (sum + derxy2_(0, inode));
      visc_shp_(1, inode) = 0.5 * derxy2_(3, inode);
      visc_shp_(2, inode) = 0.5 * derxy2_(4, inode);
      visc_shp_(3, inode) = 0.5 * derxy2_(3, inode);
      visc_shp_(4, inode) = 0.5 * (sum + derxy2_(1, inode));
      visc_shp_(5, inode) = 0.5 * derxy2_(5, inode);
      visc_shp_(6, inode) = 0.5 * derxy2_(4, inode);
      visc_shp_(7, inode) = 0.5 * derxy2_(5, inode);
      visc_shp_(8, inode) = 0.5 * (sum + derxy2_(2, inode));
    }
  }
  else if (nsd_ == 2)
  {
    for (int inode = 0; inode < nen_; ++inode)
    {
      double sum = (derxy2_(0, inode) + derxy2_(1, inode));
      visc_shp_(0, inode) = 0.5 * (sum + derxy2_(0, inode));
      visc_shp_(1, inode) = 0.5 * derxy2_(2, inode);
      visc_shp_(2, inode) = 0.5 * derxy2_(2, inode);
      visc_shp_(3, inode) = 0.5 * (sum + derxy2_(1, inode));
    }
  }
  else
    dserror("Epsilon(N) is not implemented for the 1D case");

  for (int inode = 0; inode < nen_; ++inode)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int nsd_idim = idim * nsd_;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        visc_(idim) += visc_shp_(nsd_idim + jdim, inode) * evelnp(jdim, inode);
        fluidvisc_(idim) += visc_shp_(nsd_idim + jdim, inode) * efluidvelnp(jdim, inode);
      }
    }
  }

  if (not fldAdPara_->IsStationary())
  {
    for (int inode = 0; inode < nen_; ++inode)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          visc_old_(idim) += visc_shp_(nsd_idim + jdim, inode) * eveln(jdim, inode);
          fluidvisc_new_(idim) += visc_shp_(nsd_idim + jdim, inode) * efluidvelnpp(jdim, inode);
        }
      }
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscretePSPG(
    LINALG::Matrix<nen_ * nsd_, nen_>& estif_w_q, LINALG::Matrix<nen_, nen_>& estif_r_q,
    LINALG::Matrix<nsd_, nen_>& velforce, LINALG::Matrix<nen_, 1>& preforce,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomTestStat, const double& timefacfac,
    const double& timefacfacrhs, const double& timefacfacpre, const double& timefacfacprerhs) const
{
  const double tau = tau_(1);
  const double tau_old = tau_old_(1);

  /* pressure stabilisation: inertia if not stationary*/
  /*
              /                  \
             |                    |
             |  rho*Du , nabla q  |
             |                    |
              \                  /
   */
  /* pressure stabilisation: convection, convective part */
  /*
              /                                   \
             |  /       n+1       \                |
             | |   rho*u   o nabla | Du , nabla q  |
             |  \      (i)        /                |
              \                                   /
   */
  /* pressure stabilisation: convection, reactive part if Newton */
  /*
              /                                   \
             |  /                \   n+1           |
             | |   rho*Du o nabla | u     , grad q |
             |  \                /   (i)           |
              \                                   /
   */
  /* pressure stabilisation: reaction if included */
  /*
              /                     \
             |                      |
             |  sigma*Du , nabla q  |
             |                      |
              \                    /
   */
  /* pressure stabilisation: viscosity (-L_visc_u) */
  /*
              /                              \
             |               /  \             |
         mu  |  nabla o eps | Du | , nabla q  |
             |               \  /             |
              \                              /
   */

  // stationary part
  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui_p_jdim = nsd_ * ui + jdim;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = nsd_ * idim;

        velforce(jdim, ui) -= tau * GalMomTestStat(nsd_idim + jdim, ui) * gradp_(idim);
        if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
          velforce(jdim, ui) -= tau_old * timefacfacrhs / timefacfac *
                                GalMomTestStat(nsd_idim + jdim, ui) * gradp_old_(idim);

        for (int vi = 0; vi < nen_; ++vi)
        {
          const double temp_vi_idim = derxy_(idim, vi) * tau;

          estif_w_q(fui_p_jdim, vi) += GalMomTestStat(nsd_idim + jdim, ui) * temp_vi_idim;

        }  // jdim
      }    // vi
    }      // ui
  }        // idim


  // instationary part
  if (not fldAdPara_->IsStationary())
  {
    double fac = fldAdPara_->Density() * fac_ * tau;          // fac -> mass matrix
    double fac_old = fldAdPara_->Density() * fac_ * tau_old;  // fac -> mass matrix

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double uifunct = fac * funct_(ui);
      const double uifunctold = fac_old * funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        velforce(idim, ui) -= uifunct * gradp_(idim);
        if (fldAdPara_->IsInitInstatStep() == false)
          velforce(idim, ui) += uifunctold * gradp_old_(idim);

        for (int vi = 0; vi < nen_; vi++)
        {
          estif_w_q(ui * nsd_ + idim, vi) += uifunct * derxy_(idim, vi);
        }
      }
    }  // ui
  }


  /* pressure stabilisation: pressure( L_pres_p) */
  /*
               /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
               \                    /
   */
  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double v = timefacfacpre * derxy_(idim, ui) * tau;

      preforce(ui) -= v * gradp_(idim);
      if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
        preforce(ui) -= tau_old * timefacfacprerhs * derxy_(idim, ui) * gradp_old_(idim);

      for (int vi = 0; vi < nen_; ++vi)
      {
        estif_r_q(ui, vi) += v * derxy_(idim, vi);
      }  // vi
    }    // end for(idim)
  }      // ui


  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscreteSUPG(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v,
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_r_v, LINALG::Matrix<nsd_, nen_>& velforce,
    LINALG::Matrix<nen_, 1>& preforce, const LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomTestStat,
    const double& timefacfac, const double& timefacfacrhs, const double& timefacfacpre,
    const double& timefacfacprerhs) const
{
  const double tau = tau_(0);
  const double tau_old = tau_old_(0);


  LINALG::Matrix<nen_, 1> supg_test(true);
  LINALG::Matrix<nsd_, 1> supg_vel(true);
  LINALG::Matrix<nsd_, 1> supg_vel_old(true);

  supg_test.MultiplyTN(derxy_, fluidvelint_);
  supg_vel.Multiply(vderxy_, fluidvelint_);
  supg_vel_old.Multiply(vderxy_old_, fluidvelint_old_);

  supg_test.Scale(fldAdPara_->Density() * tau);
  supg_vel.Scale(fldAdPara_->Density() * tau);
  supg_vel_old.Scale(fldAdPara_->Density() * tau_old);

  LINALG::Matrix<nsd_, 1> conv_fluidvel(true);
  LINALG::Matrix<nsd_, 1> conv_fluidvel_new(true);
  conv_fluidvel.Multiply(fluidvelxy_, fluidvelint_);
  conv_fluidvel_new.Multiply(fluidvelxy_new_, fluidvelint_new_);


  LINALG::Matrix<nsd_, 1> momres(true);
  momres.Update(1.0, fluidgradp_, reacoeff_, fluidvelint_);
  momres.Update(
      fldAdPara_->Density(), conv_fluidvel, -2.0 * fldAdPara_->Viscosity(), fluidvisc_, 1.0);
  momres.Update(-fldAdPara_->Density(), fluidbodyforce_, 1.0);

  if (fldAdPara_->IsStationary() == false)
  {
    const double fac1 = fldAdPara_->Density() / (fldAdPara_->Dt() * fldAdPara_->Theta());
    const double fac2 = fldAdPara_->OmTheta() / fldAdPara_->Theta();

    momres.Update(fac1, fluidvelint_, -fac1, fluidvelint_new_, 1.0);
    momres.Update(fac2, fluidgradp_new_, fac2 * reacoeff_, fluidvelint_new_, 1.0);
    momres.Update(fac2 * fldAdPara_->Density(), conv_fluidvel_new,
        -2.0 * fac2 * fldAdPara_->Viscosity(), fluidvisc_new_, 1.0);
    momres.Update(-fldAdPara_->Density() * fac2, fluidbodyforce_new_, 1.0);
  }

  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double w = timefacfac * fldAdPara_->Density() * tau * momres(idim) * funct_(ui);
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        velforce(jdim, ui) -= w * vderxy_(idim, jdim);

        velforce(jdim, ui) -= GalMomTestStat(nsd_ * idim + jdim, ui) * supg_vel(idim);
        if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
          velforce(jdim, ui) -= timefacfacrhs / timefacfac *
                                GalMomTestStat(nsd_ * idim + jdim, ui) * supg_vel_old(idim);

        for (int vi = 0; vi < nen_; ++vi)
        {
          estif_w_v(nsd_ * ui + jdim, nsd_ * vi + idim) +=
              GalMomTestStat(nsd_ * idim + jdim, ui) * supg_test(vi) + derxy_(jdim, vi) * w;
        }  // vi
      }    // jdim
    }      // idim
  }        // ui

  // instationary part
  if (fldAdPara_->IsStationary() == false)
  {
    double fac = fldAdPara_->Density() * fac_;  // fac -> mass matrix

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double uifunct = fac * funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        velforce(idim, ui) -= uifunct * supg_vel(idim);
        if (fldAdPara_->IsInitInstatStep() == false)
          velforce(idim, ui) += uifunct * supg_vel_old(idim);

        for (int vi = 0; vi < nen_; vi++)
        {
          estif_w_v(nsd_ * ui + idim, nsd_ * vi + idim) += uifunct * supg_test(vi);
        }
      }
    }  // ui
  }



  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double v = timefacfacpre * derxy_(idim, ui);

      preforce(ui) -= v * supg_vel(idim);
      if ((fldAdPara_->IsStationary() == false) and (fldAdPara_->IsInitInstatStep() == false))
        preforce(ui) -= v * timefacfacprerhs / timefacfacpre * supg_vel_old(idim);

      for (int vi = 0; vi < nen_; ++vi)
      {
        estif_r_v(ui, nsd_ * vi + idim) += v * supg_test(vi);
      }
    }
  }  // end for(idim)

  return;
}



/*---------------------------------------------------------------------------------*
 | compute Grad-Div stabilization terms                           winklmaier 02/14 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::DiscreteContStab(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfacdiv, const double& timefacfacdivrhs) const
{
  double cont_stab_fac = timefacfacdiv * tau_(2);

  /* continuity stabilisation */
  /*
              /                        \
             |                          |
        tauC | nabla o Du  , nabla o v  |
             |                          |
              \                        /
  */
  // continuity stabilization is symmetric, so that primal entries are equal to dual entries
  for (int vi = 0; vi < nen_; ++vi)
  {
    const int fvi = nsd_ * vi;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int fvi_jdim = fvi + jdim;
      const double v0 = cont_stab_fac * derxy_(jdim, vi);

      /* viscosity term on right-hand side */
      velforce(jdim, vi) -= cont_stab_fac * vdiv_ * derxy_(jdim, vi);
      velforce(jdim, vi) -= timefacfacdivrhs * tau_old_(2) * vdiv_old_ * derxy_(jdim, vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        const int fui = nsd_ * ui;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_w_v(fvi_jdim, fui + idim) += v0 * derxy_(idim, ui);
        }
      }
    }  // end for(idim)
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | compute PSPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PSPG(LINALG::Matrix<nen_, nen_ * nsd_>& estif_r_v,
    LINALG::Matrix<nen_, nen_>& estif_r_q, LINALG::Matrix<nen_, 1>& preforce,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomResnU,
    const LINALG::Matrix<nsd_, 1>& StrongResMomScaled, const double& timefacfac,
    const double& timefacfacrhs, const double& timefacfacpre, const double& timefacfacprerhs) const
{
  const double tau = tau_(1);

  /*
      pressure stabilization

        instationary + reactive
     /                            \
    |                              |
  - |  (rho + alpha) Dv , nabla r  |
    |                              |
     \                            /

                   convective term 1                      convective term 2
             /                          \           /                            \
            |   n                        |         |                 n            |
      + rho |  u o nabla Dv   , nabla r  | - rho * |  Dv o nabla    u  , nabla r  |
            |              (i)           |         |            (i)               |
             \                          /           \                            /

              /      viscous term                  \
             |                                      |
      + 2 mu |  nabla o epsilon   ( Dv ) , nabla r  |
             |                 (i)                  |
              \                                    /
  */

  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui_p_jdim = nsd_ * ui + jdim;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = nsd_ * idim;

        for (int vi = 0; vi < nen_; ++vi)
        {
          estif_r_v(vi, fui_p_jdim) -= tau * derxy_(idim, vi) * GalMomResnU(nsd_idim + jdim, ui);
        }  // jdim
      }    // vi
    }      // ui
  }        // idim

  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const double v = tau * timefacfacpre * derxy_(idim, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        /* pressure stabilisation: pressure( L_pres_p) */
        /*
               /                    \
              |                      |
            + |  nabla Dq , nabla r  |
              |                      |
               \                    /
         */
        estif_r_q(vi, ui) += v * derxy_(idim, vi);
      }  // vi
    }    // end for(idim)
  }      // ui

  // rhs for new and old time step
  for (int idim = 0; idim < nsd_; ++idim)
  {
    const double resmom_scaled = -tau * StrongResMomScaled(idim);

    for (int vi = 0; vi < nen_; ++vi)
    {
      // pressure stabilization
      preforce(vi) -= derxy_(idim, vi) * resmom_scaled;
    }
  }  // end for(idim)
  return;
}



/*---------------------------------------------------------------------------------*
 | compute SUPG stabilization terms                               winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::SUPG(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v,
    LINALG::Matrix<nen_ * nsd_, nen_>& estif_w_q, LINALG::Matrix<nsd_, nen_>& velforce,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& GalMomResnU,
    const LINALG::Matrix<nsd_, 1>& StrongResMomScaled, const double& timefacfac,
    const double& timefacfacrhs, const double& timefacfacpre, const double& timefacfacprerhs) const
{
  /*
     test function
                     /  n             \
  supg_test =   rho |  u o nabla w     |
                     \            (i) /
   */

  double supgfac = -fldAdPara_->Density() * tau_(0);

  LINALG::Matrix<nen_, 1> supg_test(true);
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int dim = 0; dim < nsd_; ++dim)
    {
      supg_test(vi) += supgfac * derxy_(dim, vi) * fluidvelint_(dim);
    }
  }

  /*
      supg stabilization

        instationary + reactive
     /                              \
    |                                |
    |  (rho + alpha) Dv , supg_test  |
    |                                |
     \                              /

                   convective term 1                      convective term 2
             /                            \           /                              \
            |   n                          |         |                 n              |
      - rho |  u o nabla Dv   , supg_test  | + rho * |  Dv o nabla    u  , supg_test  |
            |              (i)             |         |            (i)                 |
             \                            /           \                              /

              /      viscous term                    \
             |                                        |
      - 2 mu |  nabla o epsilon   ( Dv ) , supg_test  |
             |                 (i)                    |
              \                                      /
  */

  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int nsd_idim = nsd_ * idim;

      const int fvi_p_idim = nsd_ * vi + idim;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const int nsd_idim_p_jdim = nsd_idim + jdim;
        for (int ui = 0; ui < nen_; ++ui)
        {
          const int fui_p_jdim = nsd_ * ui + jdim;

          estif_w_v(fvi_p_idim, fui_p_jdim) += supg_test(vi) * GalMomResnU(nsd_idim_p_jdim, ui);
        }  // jdim
      }    // vi
    }      // ui
  }        // idim

  /* supg stabilisation: pressure part  ( L_pres_p) */
  /*
              /                      \
             |                        |
           - |  nabla Dq , supg_test  |
             |                        |
              \                      /
   */
  for (int vi = 0; vi < nen_; ++vi)
  {
    const double v = timefacfacpre * supg_test(vi);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fvi = nsd_ * vi + idim;

      for (int ui = 0; ui < nen_; ++ui)
      {
        estif_w_q(fvi, ui) -= v * derxy_(idim, ui);
      }
    }
  }  // end for(idim)


  // rhs for new and old time step
  for (int idim = 0; idim < nsd_; ++idim)
  {
    for (int vi = 0; vi < nen_; ++vi)
      velforce(idim, vi) -= supg_test(vi) * StrongResMomScaled(idim);
  }  // end for(idim)
  return;
}



/*---------------------------------------------------------------------------------*
 | compute residual of continuity equation                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContRes(
    double& StrongResContScaled, const double& timefacfacdiv, const double& timefacfacdivrhs) const
{
  StrongResContScaled = timefacfacdiv * (vdiv_ - contforce_);

  if (not fldAdPara_->IsStationary())
    StrongResContScaled += timefacfacdivrhs * (vdiv_old_ - contforce_old_);
}



/*---------------------------------------------------------------------------------*
 | compute div-grad (=continuity) stabilization term              winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContStab(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_w_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfacdiv, const double& timefacfacdivrhs) const
{
  double graddivfac = timefacfacdiv * tau_(2);
  double value = 0.0;

  /* continuity stabilisation on left hand side */
  /*
              /                        \
             |                          |
        tauC | nabla o Dv  , nabla o w  |
             |                          |
              \                        /
  */

  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fui_p_idim = fui + idim;

      value = graddivfac * derxy_(idim, ui);

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          estif_w_v(fvi + jdim, fui_p_idim) += value * derxy_(jdim, vi);
        }
      }
    }  // end for(idim)
  }

  double StrongResContScaled = 0.0;

  ContRes(StrongResContScaled, timefacfacdiv, timefacfacdivrhs);

  // computation of rhs viscosity term at new time step
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      /* viscosity term on right-hand side */
      velforce(idim, vi) -= tau_(2) * derxy_(idim, vi) * StrongResContScaled;
    }
  }

  return;
}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization& discretization,  ///< discretization
    const std::vector<int>& lm,                 ///<
    LINALG::Matrix<nsd_, nen_>* matrixtofill,   ///< vector field
    LINALG::Matrix<nen_, 1>* vectortofill,      ///< scalar field
    const std::string state                     ///< state of the global vector
    ) const
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

  for (int inode = 0; inode < nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for (int idim = 0; idim < nsd_; ++idim)  // number of dimensions
      {
        (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * numdofpernode_)];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL) (*vectortofill)(inode, 0) = mymatrix[nsd_ + (inode * numdofpernode_)];
  }
}
