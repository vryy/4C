/*----------------------------------------------------------------------*/
/*! \file

 \brief routines for calculation of off diagonal terms of scatra element

 \level 3

 *----------------------------------------------------------------------*/


#include "../drt_geometry/position_array.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele.H"
#include "scatra_ele_action.H"

#include "scatra_ele_calc.H"


/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvaluateOD(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  if (SetupCalc(ele, discretization) == -1) return 0;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele, params, discretization, la);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix
  //--------------------------------------------------------------------------------

  // extract action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params, "action");

  // evaluate action
  EvaluateActionOD(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvaluateActionOD(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::calc_scatra_mono_odblock_mesh:
    {
      SysmatODMesh(ele, elemat1_epetra, nsd_);

      break;
    }

    case SCATRA::calc_scatra_mono_odblock_fluid:
    {
      SysmatODFluid(ele, elemat1_epetra, nsd_ + 1);

      break;
    }

    default:
    {
      dserror("Not acting on action %i. Forgot implementation?", action);
      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::SysmatODMesh(
    DRT::Element* ele,               ///< the element those matrix is calculated
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    const int ndofpernodemesh        ///< number of DOF of mesh displacement field
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // const double vol=EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP())
  {
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);
  }

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_, 0.0);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // velocity divergence required for conservative form
    double vdiv(0.0);
    if (scatrapara_->IsConservative()) GetDivergence(vdiv, evelnp_);

    //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} =
    // J * N_x
    // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
    // i.e. det(dx/ds)
    static LINALG::Matrix<1, nsd_ * nen_> dJ_dmesh(false);
    const double J = xjm_.Determinant();
    for (unsigned i = 0; i < nen_; i++)
      for (unsigned j = 0; j < nsd_; j++) dJ_dmesh(j + i * nsd_) = J * derxy_(j, i);

    // loop all scalars
    for (int k = 0; k < numscal_; ++k)  // deal with a system of transported scalars
    {
      // reactive part of the form: (reaction coefficient)*phi
      double rea_phi(0.0);
      rea_phi = densnp[k] * scatravarmanager_->Phinp(k) * reamanager_->GetReaCoeff(k);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      GetRhsInt(rhsint, densnp[k], k);


      // subgrid-scale convective term
      LINALG::Matrix<nen_, 1> sgconv(true);
      // subgrid-scale velocity vector in gausspoint
      LINALG::Matrix<nsd_, 1> sgvelint(true);

      // residual of convection-diffusion-reaction eq
      double scatrares(0.0);

      // compute residual of scalar transport equation and
      // subgrid-scale part of scalar
      CalcStrongResidual(k, scatrares, densam[k], densnp[k], rea_phi, rhsint, tau[k]);

      double rhsfac = scatraparatimint_->TimeFacRhs() * fac;

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsStationary())
        CalcLinMassODMesh(emat, k, ndofpernodemesh, rhsfac, fac, densam[k], densnp[k],
            scatravarmanager_->Phinp(k), scatravarmanager_->Hist(k), J, dJ_dmesh);

      // the order of the following three functions is important
      // and must not be changed
      ComputeRhsInt(rhsint, densam[k], densnp[k], scatravarmanager_->Hist(k));

      // diffusive part used in stabilization terms
      LINALG::Matrix<nen_, 1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
      }

      RecomputeScatraResForRhs(scatrares, k, diff, densn[k], densnp[k], rea_phi, rhsint);

      RecomputeConvPhiForRhs(k, sgvelint, densnp[k], densn[k], vdiv);

      //----------------------------------------------------------------
      // standard Galerkin transient, old part of rhs and bodyforce term
      //----------------------------------------------------------------
      CalcHistAndSourceODMesh(emat, k, ndofpernodemesh, fac, rhsint, J, dJ_dmesh, densnp[k]);

      //----------------------------------------------------------------
      // standard Galerkin terms - convective term
      //----------------------------------------------------------------
      if (not scatrapara_->IsConservative())
        CalcConvODMesh(emat, k, ndofpernodemesh, fac, rhsfac, densnp[k], J,
            scatravarmanager_->GradPhi(k), scatravarmanager_->ConVel(k));
      else
        CalcConvConsODMesh(emat, k, fac, scatraparatimint_->TimeFac() * fac, densnp[k], J);

      //----------------------------------------------------------------
      // standard Galerkin terms  --  diffusive term
      //----------------------------------------------------------------
      CalcDiffODMesh(emat, k, ndofpernodemesh, diffmanager_->GetIsotropicDiff(k), fac, rhsfac, J,
          scatravarmanager_->GradPhi(k), scatravarmanager_->ConVel(k), dJ_dmesh);

      //----------------------------------------------------------------
      // standard Galerkin terms  -- "shapederivatives" reactive term
      //----------------------------------------------------------------
      CalcReactODMesh(emat, k, ndofpernodemesh, rhsfac, rea_phi, J, dJ_dmesh);

    }  // end loop all scalars

  }  // end loop Gauss points


  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::SysmatODFluid(
    DRT::Element* ele,               ///< the element those matrix is calculated
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    const int numdofpernode_fluid    ///< number of DOF of fluid field
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // const double vol=EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP())
  {
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);
  }

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_, 0.0);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // velocity divergence required for conservative form
    double vdiv(0.0);
    if (scatrapara_->IsConservative()) GetDivergence(vdiv, evelnp_);

    // loop all scalars
    for (int k = 0; k < numscal_; ++k)  // deal with a system of transported scalars
    {
      // reactive part of the form: (reaction coefficient)*phi
      double rea_phi(0.0);
      rea_phi = densnp[k] * scatravarmanager_->Phinp(k) * reamanager_->GetReaCoeff(k);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      GetRhsInt(rhsint, densnp[k], k);

      const double timefacfac = scatraparatimint_->TimeFac() * fac;

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsStationary())
        CalcLinMassODFluid(emat, k, numdofpernode_fluid, timefacfac, fac, densam[k], densnp[k],
            scatravarmanager_->Phinp(k), scatravarmanager_->Hist(k));

      ComputeRhsInt(rhsint, densam[k], densnp[k], scatravarmanager_->Hist(k));

      //----------------------------------------------------------------
      // standard Galerkin transient, old part of rhs and bodyforce term
      //----------------------------------------------------------------
      CalcHistAndSourceODFluid(emat, k, numdofpernode_fluid, fac, rhsint, densnp[k]);

      //----------------------------------------------------------------
      // standard Galerkin terms - convective term
      //----------------------------------------------------------------
      // calculation of convective element matrix in convective form
      CalcMatConvODFluid(
          emat, k, numdofpernode_fluid, timefacfac, densnp[k], scatravarmanager_->GradPhi(k));

      // add conservative contributions
      if (scatrapara_->IsConservative())
        CalcMatConvAddConsODFluid(
            emat, k, numdofpernode_fluid, timefacfac, densnp[k], scatravarmanager_->Phinp(k));

      //----------------------------------------------------------------
      // standard Galerkin terms  --  diffusive term
      //----------------------------------------------------------------
      CalcDiffODFluid(emat, k, numdofpernode_fluid, timefacfac, scatravarmanager_->GradPhi(k));

      //----------------------------------------------------------------
      // standard Galerkin terms  -- "shapederivatives" reactive term
      //----------------------------------------------------------------
      CalcReactODFluid(emat, k, numdofpernode_fluid, timefacfac, rea_phi);

    }  // end loop all scalars

  }  // end loop Gauss points

  return;
}


/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix                                    |
 |  in convective form (OD fluid)                                   vuong 08/14 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatConvODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodefluid,
    const double timefacfac, const double densnp, const LINALG::Matrix<nsd_, 1>& gradphi)
{
  // convective term in convective form
  const double densfac = timefacfac * densnp;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = densfac * funct_(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * ndofpernodefluid;

      for (unsigned udim = 0; udim < nsd_; ++udim)
        emat(fvi, fui + udim) += v * funct_(ui) * gradphi(udim);
    }
  }
  return;
}

/*-----------------------------------------------------------------------------*
 |   calculation of convective element matrix: add conservative                 |
 |   contributions (OD fluid)                                       vuong 08/14 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatConvAddConsODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodefluid,
    const double timefacfac, const double densnp, const double phinp)
{
  const double consfac = timefacfac * densnp * phinp;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = consfac * funct_(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * ndofpernodefluid;

      for (unsigned udim = 0; udim < nsd_; ++udim) emat(fvi, fui) += v * derxy_(udim, udim);
    }
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of linearized mass (OD mesh)             vuong 08/14 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcLinMassODMesh(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double rhsfac,
    const double fac, const double densam, const double densnp, const double phinp,
    const double hist, const double J, const LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  double vtrans = 0.0;

  if (scatraparatimint_->IsGenAlpha())
    vtrans = rhsfac * densam * hist / J;
  else
  {
    // compute scalar at integration point
    vtrans = fac * densnp * phinp / J;
  }

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    const double val = vtrans * funct_(vi);
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * ndofpernodemesh;

      for (unsigned udim = 0; udim < nsd_; ++udim)
        emat(fvi, fui + udim) += val * dJ_dmesh(fui + udim);
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term  (OD mesh)   vuong 08/14 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcHistAndSourceODMesh(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double fac,
    const double rhsint, const double J, const LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh,
    const double densnp)
{
  double vrhs = -1.0 * fac / J * rhsint;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;
    const double val = vrhs * funct_(vi);

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * ndofpernodemesh;

      for (unsigned udim = 0; udim < nsd_; ++udim)
        emat(fvi, fui + udim) += val * dJ_dmesh(fui + udim);
    }
  }
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term (OD mesh)            vuong 08/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcConvODMesh(Epetra_SerialDenseMatrix& emat,
    const int k, const int ndofpernodemesh, const double fac, const double rhsfac,
    const double densnp, const double J, const LINALG::Matrix<nsd_, 1>& gradphi,
    const LINALG::Matrix<nsd_, 1>& convelint)
{
  if (not scatraparatimint_->IsStationary())
  {
    // convective term in convective form
    const double densfac = fac * densnp;
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = densfac * funct_(vi);
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * ndofpernodemesh;

        for (unsigned udim = 0; udim < nsd_; ++udim)
          emat(fvi, fui + udim) += -1.0 * v * funct_(ui) * gradphi(udim);
      }
    }
  }

  //----------------------------------------------------------------
  // standard Galerkin terms  -- "shapederivatives" convective term
  //----------------------------------------------------------------
  ApplyShapeDerivsConv(emat, k, rhsfac, densnp, J, gradphi, convelint);

  return;
}


/*------------------------------------------------------------------------------ *
 | standard Galerkin convective term in conservative form (OD mesh)   fang 08/17 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcConvConsODMesh(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    const int k,                     //!< species index
    const double fac,                //!< domain-integration factor
    const double timefacfac,         //!< time-integration factor times domain-integration factor
    const double densnp,             //!< density
    const double J  //!< Jacobian determinant of mapping between spatial and parameter coordinates
    ) const
{
  if (not scatraparatimint_->IsStationary())
  {
    // d (phi div v_s)
    // --------------- = phi*derxy*(d_n+1)/theta/dt (plus shape derivatives, see below)
    //    d (d_n+1)
    //
    // prefactor is fac, since timefacfac/theta/dt = fac
    const double densfac = scatravarmanager_->Phinp(k) * fac * densnp;
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;
      const double v = densfac * funct_(vi);

      for (unsigned ui = 0; ui < nen_; ++ui)
        for (unsigned idim = 0; idim < nsd_; ++idim)
          emat(fvi, ui * nsd_ + idim) += v * derxy_(idim, ui);
    }
  }

  // shape derivatives associated with divergence operator
  LINALG::Matrix<nsd_, nsd_> gridvelderiv(true);
  gridvelderiv.MultiplyNT(evelnp_, deriv_);

  if (nsd_ == 3)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double v0 =
          gridvelderiv(1, 0) * (deriv_(2, ui) * xjm_(1, 2) - deriv_(1, ui) * xjm_(2, 2)) +
          gridvelderiv(1, 1) * (deriv_(0, ui) * xjm_(2, 2) - deriv_(2, ui) * xjm_(0, 2)) +
          gridvelderiv(1, 2) * (deriv_(1, ui) * xjm_(0, 2) - deriv_(0, ui) * xjm_(1, 2)) +
          gridvelderiv(2, 0) * (deriv_(1, ui) * xjm_(2, 1) - deriv_(2, ui) * xjm_(1, 1)) +
          gridvelderiv(2, 1) * (deriv_(2, ui) * xjm_(0, 1) - deriv_(0, ui) * xjm_(2, 1)) +
          gridvelderiv(2, 2) * (deriv_(0, ui) * xjm_(1, 1) - deriv_(1, ui) * xjm_(0, 1));
      const double v1 =
          gridvelderiv(0, 0) * (deriv_(1, ui) * xjm_(2, 2) - deriv_(2, ui) * xjm_(1, 2)) +
          gridvelderiv(0, 1) * (deriv_(2, ui) * xjm_(0, 2) - deriv_(0, ui) * xjm_(2, 2)) +
          gridvelderiv(0, 2) * (deriv_(0, ui) * xjm_(1, 2) - deriv_(1, ui) * xjm_(0, 2)) +
          gridvelderiv(2, 0) * (deriv_(2, ui) * xjm_(1, 0) - deriv_(1, ui) * xjm_(2, 0)) +
          gridvelderiv(2, 1) * (deriv_(0, ui) * xjm_(2, 0) - deriv_(2, ui) * xjm_(0, 0)) +
          gridvelderiv(2, 2) * (deriv_(1, ui) * xjm_(0, 0) - deriv_(0, ui) * xjm_(1, 0));
      const double v2 =
          gridvelderiv(0, 0) * (deriv_(2, ui) * xjm_(1, 1) - deriv_(1, ui) * xjm_(2, 1)) +
          gridvelderiv(0, 1) * (deriv_(0, ui) * xjm_(2, 1) - deriv_(2, ui) * xjm_(0, 1)) +
          gridvelderiv(0, 2) * (deriv_(1, ui) * xjm_(0, 1) - deriv_(0, ui) * xjm_(1, 1)) +
          gridvelderiv(1, 0) * (deriv_(1, ui) * xjm_(2, 0) - deriv_(2, ui) * xjm_(1, 0)) +
          gridvelderiv(1, 1) * (deriv_(2, ui) * xjm_(0, 0) - deriv_(0, ui) * xjm_(2, 0)) +
          gridvelderiv(1, 2) * (deriv_(0, ui) * xjm_(1, 0) - deriv_(1, ui) * xjm_(0, 0));

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const int fvi = vi * numdofpernode_ + k;
        const double v = scatravarmanager_->Phinp(k) * timefacfac / J * funct_(vi);

        emat(fvi, ui * 3) += v * v0;
        emat(fvi, ui * 3 + 1) += v * v1;
        emat(fvi, ui * 3 + 2) += v * v2;
      }
    }
  }

  else if (nsd_ == 2)
  {
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;
      const double v = scatravarmanager_->Phinp(k) * timefacfac / J * funct_(vi);

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        emat(fvi, ui * 2) +=
            v * (-gridvelderiv(1, 0) * deriv_(1, ui) + gridvelderiv(1, 1) * deriv_(0, ui));
        emat(fvi, ui * 2 + 1) +=
            v * (gridvelderiv(0, 0) * deriv_(1, ui) - gridvelderiv(0, 1) * deriv_(0, ui));
      }
    }
  }

  else
    dserror("Shape derivatives not implemented for 1D problems!");

  return;
}


/*-------------------------------------------------------------------- *
 |  "shapederivatives" convective term (OD mesh)           vuong 08/14 |
 |  put it into its own function                      kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ApplyShapeDerivsConv(
    Epetra_SerialDenseMatrix& emat, const int k, const double rhsfac, const double densnp,
    const double J, const LINALG::Matrix<nsd_, 1>& gradphi,
    const LINALG::Matrix<nsd_, 1>& convelint)
{
  //----------------------------------------------------------------
  // standard Galerkin terms  -- "shapederivatives" convective term
  //----------------------------------------------------------------

  if (nsd_ == 3)
  {
    const double xjm_0_0 = xjm_(0, 0);
    const double xjm_0_1 = xjm_(0, 1);
    const double xjm_0_2 = xjm_(0, 2);
    const double xjm_1_0 = xjm_(1, 0);
    const double xjm_1_1 = xjm_(1, 1);
    const double xjm_1_2 = xjm_(1, 2);
    const double xjm_2_0 = xjm_(2, 0);
    const double xjm_2_1 = xjm_(2, 1);
    const double xjm_2_2 = xjm_(2, 2);

    {
      // gradient of scalar w.r.t. reference coordinates
      static LINALG::Matrix<nsd_, 1> refgradphi;
      refgradphi.Multiply(xjm_, gradphi);

      const double refgradphi_0 = refgradphi(0);
      const double refgradphi_1 = refgradphi(1);
      const double refgradphi_2 = refgradphi(2);

      const double convelint_0 = convelint(0);
      const double convelint_1 = convelint(1);
      const double convelint_2 = convelint(2);

      const double vrhs = rhsfac * densnp / J;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +convelint_1 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_2 - deriv_(1, ui) * xjm_2_2) +
                               refgradphi_1 * (deriv_(0, ui) * xjm_2_2 - deriv_(2, ui) * xjm_0_2) +
                               refgradphi_2 * (deriv_(1, ui) * xjm_0_2 - deriv_(0, ui) * xjm_1_2)) +
            convelint_2 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_1 - deriv_(2, ui) * xjm_1_1) +
                              refgradphi_1 * (deriv_(2, ui) * xjm_0_1 - deriv_(0, ui) * xjm_2_1) +
                              refgradphi_2 * (deriv_(0, ui) * xjm_1_1 - deriv_(1, ui) * xjm_0_1));
        const double v01 =
            +convelint_0 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_2 - deriv_(2, ui) * xjm_1_2) +
                               refgradphi_1 * (deriv_(2, ui) * xjm_0_2 - deriv_(0, ui) * xjm_2_2) +
                               refgradphi_2 * (deriv_(0, ui) * xjm_1_2 - deriv_(1, ui) * xjm_0_2)) +
            convelint_2 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_0 - deriv_(1, ui) * xjm_2_0) +
                              refgradphi_1 * (deriv_(0, ui) * xjm_2_0 - deriv_(2, ui) * xjm_0_0) +
                              refgradphi_2 * (deriv_(1, ui) * xjm_0_0 - deriv_(0, ui) * xjm_1_0));
        const double v02 =
            +convelint_0 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_1 - deriv_(1, ui) * xjm_2_1) +
                               refgradphi_1 * (deriv_(0, ui) * xjm_2_1 - deriv_(2, ui) * xjm_0_1) +
                               refgradphi_2 * (deriv_(1, ui) * xjm_0_1 - deriv_(0, ui) * xjm_1_1)) +
            convelint_1 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_0 - deriv_(2, ui) * xjm_1_0) +
                              refgradphi_1 * (deriv_(2, ui) * xjm_0_0 - deriv_(0, ui) * xjm_2_0) +
                              refgradphi_2 * (deriv_(0, ui) * xjm_1_0 - deriv_(1, ui) * xjm_0_0));

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;
          const double v = vrhs * funct_(vi);

          emat(fvi, ui * 3 + 0) += v * v00;
          emat(fvi, ui * 3 + 1) += v * v01;
          emat(fvi, ui * 3 + 2) += v * v02;
        }
      }
    }
  }
  else if (nsd_ == 2)
  {
    {
      // gradient of scalar w.r.t. reference coordinates
      static LINALG::Matrix<nsd_, 1> refgradphi;
      refgradphi.Multiply(xjm_, gradphi);

      const double refgradphi_0 = refgradphi(0);
      const double refgradphi_1 = refgradphi(1);

      const double convelint_0 = convelint(0);
      const double convelint_1 = convelint(1);

      const double vrhs = rhsfac * densnp / J;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +convelint_1 * (-refgradphi_0 * deriv_(1, ui) + refgradphi_1 * deriv_(0, ui));
        const double v01 =
            +convelint_0 * (refgradphi_0 * deriv_(1, ui) - refgradphi_1 * deriv_(0, ui));

        for (unsigned vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;
          const double v = vrhs * funct_(vi);

          emat(fvi, ui * 2 + 0) += v * v00;
          emat(fvi, ui * 2 + 1) += v * v01;
        }
      }
    }
  }
  else
    dserror("shapederivatives not implemented for 1D!");


  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term (OD mesh)   vuong 08/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcDiffODMesh(Epetra_SerialDenseMatrix& emat,
    const int k, const int ndofpernodemesh, const double diffcoeff, const double fac,
    const double rhsfac, const double J, const LINALG::Matrix<nsd_, 1>& gradphi,
    const LINALG::Matrix<nsd_, 1>& convelint, const LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  const double vrhs = -rhsfac / J * diffcoeff;

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf, gradphi, vi);
    const double val = vrhs * laplawf;

    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * ndofpernodemesh;

      for (unsigned udim = 0; udim < nsd_; ++udim)
        emat(fvi, fui + udim) += val * dJ_dmesh(fui + udim);
    }
  }

  //----------------------------------------------------------------
  // standard Galerkin terms  -- "shapederivatives" diffusive term
  //----------------------------------------------------------------

  if (nsd_ == 3)
  {
    const double xjm_0_0 = xjm_(0, 0);
    const double xjm_0_1 = xjm_(0, 1);
    const double xjm_0_2 = xjm_(0, 2);
    const double xjm_1_0 = xjm_(1, 0);
    const double xjm_1_1 = xjm_(1, 1);
    const double xjm_1_2 = xjm_(1, 2);
    const double xjm_2_0 = xjm_(2, 0);
    const double xjm_2_1 = xjm_(2, 1);
    const double xjm_2_2 = xjm_(2, 2);
    {
      const double v = diffcoeff * rhsfac / J;

      const double gradphi_0 = gradphi(0);
      const double gradphi_1 = gradphi(1);
      const double gradphi_2 = gradphi(2);

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const double deriv_vi_0 = deriv_(0, vi);
        const double deriv_vi_1 = deriv_(1, vi);
        const double deriv_vi_2 = deriv_(2, vi);

        const int fvi = vi * numdofpernode_ + k;

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +gradphi_1 * (deriv_vi_0 * (deriv_(2, ui) * xjm_1_2 - deriv_(1, ui) * xjm_2_2) +
                               deriv_vi_1 * (deriv_(0, ui) * xjm_2_2 - deriv_(2, ui) * xjm_0_2) +
                               deriv_vi_2 * (deriv_(1, ui) * xjm_0_2 - deriv_(0, ui) * xjm_1_2)) +
              gradphi_2 * (deriv_vi_0 * (deriv_(1, ui) * xjm_2_1 - deriv_(2, ui) * xjm_1_1) +
                              deriv_vi_1 * (deriv_(2, ui) * xjm_0_1 - deriv_(0, ui) * xjm_2_1) +
                              deriv_vi_2 * (deriv_(0, ui) * xjm_1_1 - deriv_(1, ui) * xjm_0_1));
          const double v01 =
              +gradphi_0 * (deriv_vi_0 * (deriv_(1, ui) * xjm_2_2 - deriv_(2, ui) * xjm_1_2) +
                               deriv_vi_1 * (deriv_(2, ui) * xjm_0_2 - deriv_(0, ui) * xjm_2_2) +
                               deriv_vi_2 * (deriv_(0, ui) * xjm_1_2 - deriv_(1, ui) * xjm_0_2)) +
              gradphi_2 * (deriv_vi_0 * (deriv_(2, ui) * xjm_1_0 - deriv_(1, ui) * xjm_2_0) +
                              deriv_vi_1 * (deriv_(0, ui) * xjm_2_0 - deriv_(2, ui) * xjm_0_0) +
                              deriv_vi_2 * (deriv_(1, ui) * xjm_0_0 - deriv_(0, ui) * xjm_1_0));
          const double v02 =
              +gradphi_0 * (deriv_vi_0 * (deriv_(2, ui) * xjm_1_1 - deriv_(1, ui) * xjm_2_1) +
                               deriv_vi_1 * (deriv_(0, ui) * xjm_2_1 - deriv_(2, ui) * xjm_0_1) +
                               deriv_vi_2 * (deriv_(1, ui) * xjm_0_1 - deriv_(0, ui) * xjm_1_1)) +
              gradphi_1 * (deriv_vi_0 * (deriv_(1, ui) * xjm_2_0 - deriv_(2, ui) * xjm_1_0) +
                              deriv_vi_1 * (deriv_(2, ui) * xjm_0_0 - deriv_(0, ui) * xjm_2_0) +
                              deriv_vi_2 * (deriv_(0, ui) * xjm_1_0 - deriv_(1, ui) * xjm_0_0));

          emat(fvi, ui * ndofpernodemesh + 0) += v * v00;
          emat(fvi, ui * ndofpernodemesh + 1) += v * v01;
          emat(fvi, ui * ndofpernodemesh + 2) += v * v02;
        }
      }
    }

    const double v = diffcoeff * rhsfac / J;

    // gradient of scalar w.r.t. reference coordinates
    static LINALG::Matrix<nsd_, 1> refgradphi;
    refgradphi.Multiply(xjm_, gradphi);

    const double refgradphi_0 = refgradphi(0);
    const double refgradphi_1 = refgradphi(1);
    const double refgradphi_2 = refgradphi(2);

    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double derxy_vi_0 = derxy_(0, vi);
      const double derxy_vi_1 = derxy_(1, vi);
      const double derxy_vi_2 = derxy_(2, vi);

      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +derxy_vi_1 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_2 - deriv_(1, ui) * xjm_2_2) +
                              refgradphi_1 * (deriv_(0, ui) * xjm_2_2 - deriv_(2, ui) * xjm_0_2) +
                              refgradphi_2 * (deriv_(1, ui) * xjm_0_2 - deriv_(0, ui) * xjm_1_2)) +
            derxy_vi_2 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_1 - deriv_(2, ui) * xjm_1_1) +
                             refgradphi_1 * (deriv_(2, ui) * xjm_0_1 - deriv_(0, ui) * xjm_2_1) +
                             refgradphi_2 * (deriv_(0, ui) * xjm_1_1 - deriv_(1, ui) * xjm_0_1));
        const double v01 =
            +derxy_vi_0 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_2 - deriv_(2, ui) * xjm_1_2) +
                              refgradphi_1 * (deriv_(2, ui) * xjm_0_2 - deriv_(0, ui) * xjm_2_2) +
                              refgradphi_2 * (deriv_(0, ui) * xjm_1_2 - deriv_(1, ui) * xjm_0_2)) +
            derxy_vi_2 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_0 - deriv_(1, ui) * xjm_2_0) +
                             refgradphi_1 * (deriv_(0, ui) * xjm_2_0 - deriv_(2, ui) * xjm_0_0) +
                             refgradphi_2 * (deriv_(1, ui) * xjm_0_0 - deriv_(0, ui) * xjm_1_0));
        const double v02 =
            +derxy_vi_0 * (refgradphi_0 * (deriv_(2, ui) * xjm_1_1 - deriv_(1, ui) * xjm_2_1) +
                              refgradphi_1 * (deriv_(0, ui) * xjm_2_1 - deriv_(2, ui) * xjm_0_1) +
                              refgradphi_2 * (deriv_(1, ui) * xjm_0_1 - deriv_(0, ui) * xjm_1_1)) +
            derxy_vi_1 * (refgradphi_0 * (deriv_(1, ui) * xjm_2_0 - deriv_(2, ui) * xjm_1_0) +
                             refgradphi_1 * (deriv_(2, ui) * xjm_0_0 - deriv_(0, ui) * xjm_2_0) +
                             refgradphi_2 * (deriv_(0, ui) * xjm_1_0 - deriv_(1, ui) * xjm_0_0));

        emat(fvi, ui * ndofpernodemesh + 0) += v * v00;
        emat(fvi, ui * ndofpernodemesh + 1) += v * v01;
        emat(fvi, ui * ndofpernodemesh + 2) += v * v02;
      }
    }
  }
  else if (nsd_ == 2)
  {
    {
      const double v = diffcoeff * rhsfac / J;

      const double gradphi_0 = gradphi(0);
      const double gradphi_1 = gradphi(1);

      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const double deriv_vi_0 = deriv_(0, vi);
        const double deriv_vi_1 = deriv_(1, vi);

        const int fvi = vi * numdofpernode_ + k;

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const double v00 =
              +gradphi_1 * (-deriv_vi_0 * deriv_(1, ui) + deriv_vi_1 * deriv_(0, ui));
          const double v01 = +gradphi_0 * (deriv_vi_0 * deriv_(1, ui) - deriv_vi_1 * deriv_(0, ui));

          emat(fvi, ui * ndofpernodemesh + 0) += v * v00;
          emat(fvi, ui * ndofpernodemesh + 1) += v * v01;
        }
      }
    }

    const double v = diffcoeff * rhsfac / J;

    // gradient of scalar w.r.t. reference coordinates
    static LINALG::Matrix<nsd_, 1> refgradphi;
    refgradphi.Multiply(xjm_, gradphi);

    const double refgradphi_0 = refgradphi(0);
    const double refgradphi_1 = refgradphi(1);

    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double derxy_vi_0 = derxy_(0, vi);
      const double derxy_vi_1 = derxy_(1, vi);

      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const double v00 =
            +derxy_vi_1 * (-refgradphi_0 * deriv_(1, ui) + refgradphi_1 * deriv_(0, ui));
        const double v01 =
            +derxy_vi_0 * (refgradphi_0 * deriv_(1, ui) - refgradphi_1 * deriv_(0, ui));

        emat(fvi, ui * ndofpernodemesh + 0) += v * v00;
        emat(fvi, ui * ndofpernodemesh + 1) += v * v01;
      }
    }
  }
  else
    dserror("shapederivatives not implemented for 1D!");
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin reactive term (OD mesh)              vuong 08/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcReactODMesh(Epetra_SerialDenseMatrix& emat,
    const int k, const int ndofpernodemesh, const double rhsfac, const double rea_phi,
    const double J, const LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  if (reamanager_->Active())
  {
    // standard Galerkin term
    double vrhs = rhsfac * rea_phi / J;

    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;

      const double val = vrhs * funct_(vi);
      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * ndofpernodemesh;

        for (unsigned udim = 0; udim < nsd_; ++udim)
          emat(fvi, fui + udim) += val * dJ_dmesh(fui + udim);
      }
    }

    //        // reactive rhs stabilization
    //        if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
    //        {
    //          vrhs =
    //          scatrapara_->USFEMGLSFac()*rhstaufac*densnp*reamanager_->GetReaCoeff(k)*scatrares;
    //          for (unsigned vi=0; vi<nen_; ++vi)
    //          {
    //            const int fvi = vi*numdofpernode_+k;
    //
    //            erhs[fvi] -= vrhs*funct_(vi);
    //          }
    //        }
  }
}

/*------------------------------------------------------------------- *
 |  calculation of linearized mass (OD fluid)        kremheller 07/17 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcLinMassODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double rhsfac,
    const double fac, const double densam, const double densnp, const double phinp,
    const double hist)
{
  // do nothing
  return;
}

/*------------------------------------------------------------------------------------------ *
 |  standard Galerkin transient, old part of rhs and source term (OD fluid) kremheller 07/17 |
 *-------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcHistAndSourceODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double fac,
    const double rhsint, const double densnp)
{
  // do nothing
  return;
}

/*------------------------------------------------------------------ *
 |  standard Galerkin reactive term (OD fluid)  kremheller 07/17     |
 *----------------------------------------------------------------   */
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcReactODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodemesh, const double rhsfac,
    const double rea_phi)
{
  // do nothing
  return;
}

/*------------------------------------------------------------------ *
 |  standard Galerkin diffusive term (OD fluid)  kremheller 07/17    |
 *----------------------------------------------------------------   */
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcDiffODFluid(
    Epetra_SerialDenseMatrix& emat,  //!< element current to be filled
    const int k,                     //!< index of current scalar
    const int ndofpernodemesh,       //!< number of dofs per node of ale element
    const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
    const LINALG::Matrix<nsd_, 1>& gradphi  //!< scalar gradient at Gauss point
)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

#include "scatra_ele_calc_fwd.hpp"
