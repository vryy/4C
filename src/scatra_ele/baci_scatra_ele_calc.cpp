/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of scatra element

\level 1


*----------------------------------------------------------------------*/

#include "baci_scatra_ele_calc.hpp"

#include "baci_discretization_fem_general_utils_boundary_integration.hpp"
#include "baci_discretization_fem_general_utils_gder2.hpp"
#include "baci_fluid_rotsym_periodicbc.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_condition_utils.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_mat_electrode.hpp"
#include "baci_mat_list.hpp"
#include "baci_mat_newtonianfluid.hpp"
#include "baci_mat_scatra.hpp"
#include "baci_mat_scatra_multiscale.hpp"
#include "baci_nurbs_discret_nurbs_utils.hpp"
#include "baci_scatra_ele.hpp"
#include "baci_scatra_ele_parameter_std.hpp"
#include "baci_scatra_ele_parameter_timint.hpp"
#include "baci_scatra_ele_parameter_turbulence.hpp"
#include "baci_utils_function.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
    const int numdofpernode, const int numscal, const std::string& disname)
    : numdofpernode_(numdofpernode),
      numscal_(numscal),
      scatrapara_(DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname)),
      turbparams_(DRT::ELEMENTS::ScaTraEleParameterTurbulence::Instance(disname)),
      scatraparatimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname)),
      diffmanager_(Teuchos::rcp(new ScaTraEleDiffManager(numscal_))),
      reamanager_(Teuchos::rcp(new ScaTraEleReaManager(numscal_))),
      ephin_(numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true)),
      ephinp_(numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true)),
      ehist_(numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true)),
      fsphinp_(numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true)),
      rotsymmpbc_(Teuchos::rcp(new FLD::RotationallySymmetricPeriodicBC<distype, nsd_ + 1,
          DRT::ELEMENTS::Fluid::none>())),
      evelnp_(true),
      econvelnp_(true),
      efsvel_(true),
      eaccnp_(true),
      edispnp_(true),
      eprenp_(true),
      tpn_(0.0),
      xsi_(true),
      xyze_(true),
      funct_(true),
      deriv_(true),
      deriv2_(true),
      derxy_(true),
      derxy2_(true),
      xjm_(true),
      xij_(true),
      xder2_(true),
      bodyforce_(numdofpernode_),
      weights_(true),
      myknots_(nsd_),
      eid_(0),
      ele_(nullptr),
      scatravarmanager_(Teuchos::rcp(new ScaTraEleInternalVariableManager<nsd_, nen_>(numscal_)))
{
  dsassert(
      nsd_ >= nsd_ele_, "problem dimension has to be equal or larger than the element dimension!");

  // safety checks related with turbulence
  if (scatrapara_->ASSGD() and turbparams_->FSSGD())
  {
    dserror(
        "No combination of all-scale and fine-scale subgrid-diffusivity approach currently "
        "possible!");
  }
  if (turbparams_->BD_Gp() and not scatrapara_->MatGP())
  {
    dserror(
        "Evaluation of B and D at Gauss point should always be combined with material evaluation "
        "at Gauss point!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::SetupCalc(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // get element coordinates
  ReadElementCoordinates(ele);

  // Now do the nurbs specific stuff (for isogeometric elements)
  if (DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization, ele, myknots_, weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if (zero_size) return -1;
  }  // Nurbs specific stuff

  // set element id
  eid_ = ele->Id();
  // set element
  ele_ = ele;

  // rotationally symmetric periodic bc's: do setup for current element
  rotsymmpbc_->Setup(ele);

  Teuchos::RCP<MAT::Material> material = ele->Material();
  if (material->MaterialType() == INPAR::MAT::m_matlist or
      material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    const Teuchos::RCP<const MAT::MatList>& material_list =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (material_list->NumMat() < numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0; k < numscal_; ++k)
    {
      const int material_id = material_list->MatID(k);

      if (material_list->MaterialById(material_id)->MaterialType() == INPAR::MAT::m_scatra)
      {
        Teuchos::RCP<const MAT::ScatraMat> single_material =
            Teuchos::rcp_static_cast<const MAT::ScatraMat>(
                material_list->MaterialById(material_id));
        scatravarmanager_->SetReactsToForce(single_material->ReactsToExternalForce(), k);
      }
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
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
  // prepare turbulence models
  //--------------------------------------------------------------------------------

  int nlayer = 0;
  ExtractTurbulenceApproach(ele, params, discretization, la, nlayer);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(ele, elemat1_epetra, elevec1_epetra, elevec2_epetra);

  // perform finite difference check on element level
  if (scatrapara_->FDCheck() == INPAR::SCATRA::fdcheck_local and
      ele->Owner() == discretization.Comm().MyPID())
    FDCheck(ele, elemat1_epetra, elevec1_epetra, elevec2_epetra);

  // ---------------------------------------------------------------------
  // output values of Prt, diffeff and Cs_delta_sq_Prt (channel flow only)
  // ---------------------------------------------------------------------

  if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky and turbparams_->CsAv())
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    StoreModelParametersForOutput(
        ele, ele->Owner() == discretization.Comm().MyPID(), turbulencelist, nlayer);
  }

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ExtractElementAndNodeValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // get number of dofset associated with velocity related dofs
  const int ndsvel = scatrapara_->NdsVel();

  // get convective (velocity - mesh displacement) velocity at nodes
  auto convel = discretization.GetState(ndsvel, "convective velocity field");
  if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

  // determine number of velocity related dofs per node
  const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

  // construct location vector for velocity related dofs
  std::vector<int> lmvel(nsd_ * nen_, -1);
  for (unsigned inode = 0; inode < nen_; ++inode)
    for (unsigned idim = 0; idim < nsd_; ++idim)
      lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

  if (scatrapara_->HasExternalForce())
  {
    auto force_velocity = discretization.GetState(ndsvel, "force_velocity");
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(
        *force_velocity, eforcevelocity_, lmvel);
  }

  // extract local values of convective velocity field from global state vector
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*convel, econvelnp_, lmvel);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  rotsymmpbc_->RotateMyValuesIfNecessary(econvelnp_);

  // get additional state vector for ALE case: grid displacement
  if (scatrapara_->IsAle())
  {
    // get velocity at nodes
    Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(ndsvel, "velocity field");
    if (vel == Teuchos::null) dserror("Cannot get state vector velocity");

    // extract local values of velocity field from global state vector
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*vel, evelnp_, lmvel);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    rotsymmpbc_->RotateMyValuesIfNecessary(evelnp_);

    // get number of dofset associated with displacement related dofs
    const int ndsdisp = scatrapara_->NdsDisp();

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(nsd_ * nen_, -1);
    for (unsigned inode = 0; inode < nen_; ++inode)
      for (unsigned idim = 0; idim < nsd_; ++idim)
        lmdisp[inode * nsd_ + idim] = la[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*dispnp, edispnp_, lmdisp);

    // add nodal displacements to point coordinates
    UpdateNodeCoordinates();
  }
  else
  {
    edispnp_.Clear();

    // velocity = convective velocity for the non-ale case
    evelnp_ = econvelnp_;
  }

  // get data required for subgrid-scale velocity: acceleration and pressure
  if (scatrapara_->RBSubGrVel())
  {
    // get acceleration values at nodes
    const Teuchos::RCP<const Epetra_Vector> acc =
        discretization.GetState(ndsvel, "acceleration field");
    if (acc == Teuchos::null) dserror("Cannot get state vector acceleration field");

    // extract local values of acceleration field from global state vector
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*acc, eaccnp_, lmvel);

    // rotate the vector field in the case of rotationally symmetric boundary conditions
    rotsymmpbc_->RotateMyValuesIfNecessary(eaccnp_);

    // construct location vector for pressure dofs
    std::vector<int> lmpre(nen_, -1);
    for (unsigned inode = 0; inode < nen_; ++inode)
      lmpre[inode] = la[ndsvel].lm_[inode * numveldofpernode + nsd_];

    // extract local values of pressure field from global state vector
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*convel, eprenp_, lmpre);
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist == Teuchos::null || phinp == Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");

  // values of scatra field are always in first dofset
  const std::vector<int>& lm = la[0].lm_;
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*hist, ehist_, lm);
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phinp, ephinp_, lm);

  if (scatraparatimint_->IsGenAlpha() and not scatraparatimint_->IsIncremental())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    if (phin == Teuchos::null) dserror("Cannot get state vector 'phin'");
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phin, ephin_, lm);
  }

  // set reaction coefficient
  if (params.isParameter("rea_coeff")) reamanager_->SetReaCoeff(params.get<double>("rea_coeff"), 0);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  OtherNodeBasedSourceTerms(lm, discretization, params);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ExtractTurbulenceApproach(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, int& nlayer)
{
  if (turbparams_->TurbModel() != INPAR::FLUID::no_model or
      (scatraparatimint_->IsIncremental() and turbparams_->FSSGD()))
  {
    // do some checks first
    if (numscal_ != 1 or numdofpernode_ != 1)
      dserror("For the time being, turbulence approaches only support one scalar field!");
  }

  // set turbulent Prandt number to value given in parameterlist
  tpn_ = turbparams_->TPN();

  // if we have a dynamic model,we overwrite this value by a local element-based one here
  if (turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    Teuchos::ParameterList& turbulencelist = params.sublist("TURBULENCE MODEL");
    // remark: for dynamic estimation, this returns (Cs*h)^2 / Pr_t
    Teuchos::RCP<Epetra_Vector> ele_prt =
        turbulencelist.get<Teuchos::RCP<Epetra_Vector>>("col_ele_Prt");
    const int id = ele->LID();
    tpn_ = (*ele_prt)[id];

    // when no averaging was done, we just keep the calculated (clipped) value
    if (turbparams_->CsAv())
      GetMeanPrtOfHomogenousDirection(params.sublist("TURBULENCE MODEL"), nlayer);
  }

  // get fine-scale values
  if ((scatraparatimint_->IsIncremental() and
          (turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_all or
              turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small)) or
      turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    // get fine scale scalar field
    Teuchos::RCP<const Epetra_Vector> gfsphinp = discretization.GetState("fsphinp");
    if (gfsphinp == Teuchos::null) dserror("Cannot get state vector 'fsphinp'");

    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*gfsphinp, fsphinp_, la[0].lm_);

    if (turbparams_->WhichFssgd() == INPAR::SCATRA::fssugrdiff_smagorinsky_small or
        turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      // get number of dofset associated with velocity-related dofs
      const int ndsvel = scatrapara_->NdsVel();

      // get fine-scale velocity at nodes
      const Teuchos::RCP<const Epetra_Vector> fsvelocity =
          discretization.GetState(ndsvel, "fine-scale velocity field");
      if (fsvelocity == Teuchos::null)
        dserror("Cannot get fine-scale velocity field from scatra discretization!");

      // determine number of velocity related dofs per node
      const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

      // construct location vector for velocity related dofs
      std::vector<int> lmvel(nsd_ * nen_, -1);
      for (unsigned inode = 0; inode < nen_; ++inode)
        for (unsigned idim = 0; idim < nsd_; ++idim)
          lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

      // extract local values of fine-scale velocity field from global state vector
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nsd_, nen_>>(*fsvelocity, efsvel_, lmvel);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::Sysmat(DRT::Element* ele,
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs,
    CORE::LINALG::SerialDenseVector& subgrdiff)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  const double vol = EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(numscal_, 1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center are also necessary
  // even if the stabilization parameter is evaluated at the element center
  if (not scatrapara_->MatGP())
  {
    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    GetMaterialParams(ele, densn, densnp, densam, visc);
  }

  //----------------------------------------------------------------------
  // calculation of subgrid diffusivity and stabilization parameter(s)
  // at element center
  //----------------------------------------------------------------------

  // the stabilization parameters (one per transported scalar)
  std::vector<double> tau(numscal_, 0.0);
  // stabilization parameters for the external force term (one per transported scalar)
  // std::vector<double> tau_force(numscal_, 0.0);
  // subgrid-scale diffusion coefficient
  double sgdiff(0.0);

  if (not scatrapara_->TauGP())
  {
    for (int k = 0; k < numscal_; ++k)  // loop of each transported scalar
    {
      // get velocity at element center
      CORE::LINALG::Matrix<nsd_, 1> convelint = scatravarmanager_->ConVel(k);

      // calculation of all-scale subgrid diffusivity (by, e.g.,
      // Smagorinsky model) at element center
      if (turbparams_->TurbModel() == INPAR::FLUID::smagorinsky or
          turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky or
          turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
      {
        CalcSubgrDiff(visc, vol, k, densnp[k]);
      }

      // calculation of fine-scale artificial subgrid diffusivity at element center
      if (turbparams_->FSSGD())
        CalcFineScaleSubgrDiff(sgdiff, subgrdiff, ele, vol, k, densnp[k],
            diffmanager_->GetIsotropicDiff(k), convelint);

      // calculation of stabilization parameter at element center
      CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
          reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k], convelint,
          vol);
    }
  }

  // prepare multifractal subgrid-scale modeling
  // calculation of model coefficients B (velocity) and D (scalar)
  // at element center
  // coefficient B of fine-scale velocity
  CORE::LINALG::Matrix<nsd_, 1> B_mfs(true);
  // coefficient D of fine-scale scalar
  double D_mfs = 0.0;
  if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not turbparams_->BD_Gp())
    {
      // make sure to get material parameters at element center
      // hence, determine them if not yet available
      if (scatrapara_->MatGP())
      {
        SetInternalVariablesForMatAndRHS();

        GetMaterialParams(ele, densn, densnp, densam, visc);
      }

      // provide necessary velocities and gradients at element center
      // get velocity at element center
      CORE::LINALG::Matrix<nsd_, 1> fsvelint(true);
      fsvelint.Multiply(efsvel_, funct_);

      // calculate model coefficients
      for (int k = 0; k < numscal_; ++k)  // loop of each transported scalar
        CalcBAndDForMultifracSubgridScales(B_mfs, D_mfs, vol, k, densnp[k],
            diffmanager_->GetIsotropicDiff(k), visc, scatravarmanager_->ConVel(k), fsvelint);
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // velocity divergence required for conservative form
    double vdiv(0.0);
    if (scatrapara_->IsConservative()) GetDivergence(vdiv, evelnp_);

    // get fine-scale velocity and its derivatives at integration point
    CORE::LINALG::Matrix<nsd_, 1> fsvelint(true);
    if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
      fsvelint.Multiply(efsvel_, funct_);

    // loop all scalars
    for (int k = 0; k < numscal_; ++k)  // deal with a system of transported scalars
    {
      // reactive part of the form: (reaction coefficient)*phi
      double rea_phi(0.0);
      rea_phi = densnp[k] * scatravarmanager_->Phinp(k) * reamanager_->GetReaCoeff(k);

      // compute gradient of fine-scale part of scalar value
      CORE::LINALG::Matrix<nsd_, 1> fsgradphi(true);
      if (turbparams_->FSSGD()) fsgradphi.Multiply(derxy_, fsphinp_[k]);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      GetRhsInt(rhsint, densnp[k], k);

      //--------------------------------------------------------------------
      // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
      // velocity and stabilization parameter(s) at integration point
      //--------------------------------------------------------------------

      // subgrid-scale convective term
      CORE::LINALG::Matrix<nen_, 1> sgconv(true);
      // subgrid-scale velocity vector in gausspoint
      CORE::LINALG::Matrix<nsd_, 1> sgvelint(true);

      double scatrares(0.0);
      // calculate strong residual
      CalcStrongResidual(k, scatrares, densam[k], densnp[k], rea_phi, rhsint, tau[k]);

      if (scatrapara_->TauGP())
      {
        // artificial diffusion / shock capturing: adaption of diffusion coefficient
        if (scatrapara_->ASSGD())
        {
          // pre-calculation of stabilization parameter at integration point need for some forms of
          // artificial diffusion
          CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
              reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k],
              scatravarmanager_->ConVel(k), vol);

          // compute artificial diffusion
          CalcArtificialDiff(vol, k, densnp[k], scatravarmanager_->ConVel(k),
              scatravarmanager_->GradPhi(k), scatravarmanager_->ConvPhi(k), scatrares, tau[k]);

          // recompute strong residual since now diffus_new = diffus_old + artdiff
          CalcStrongResidual(k, scatrares, densam[k], densnp[k], rea_phi, rhsint,
              tau[k]);  // TODO:(Thon) do we really have to do this??
        }

        // calculation of all-scale subgrid diffusivity (by, e.g.,
        // Smagorinsky model) at element center
        if (turbparams_->TurbModel() == INPAR::FLUID::smagorinsky or
            turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky or
            turbparams_->TurbModel() == INPAR::FLUID::dynamic_vreman)
        {
          CalcSubgrDiff(visc, vol, k, densnp[k]);

          // recompute strong residual since now diffus_new = diffus_old + sgdiff
          CalcStrongResidual(k, scatrares, densam[k], densnp[k], rea_phi, rhsint, tau[k]);
        }

        // calculation of fine-scale artificial subgrid diffusivity at element center
        if (turbparams_->FSSGD())
          CalcFineScaleSubgrDiff(sgdiff, subgrdiff, ele, vol, k, densnp[k],
              diffmanager_->GetIsotropicDiff(k), scatravarmanager_->ConVel(k));

        // calculation of subgrid-scale velocity at integration point if required
        if (scatrapara_->RBSubGrVel())
        {
          // calculation of stabilization parameter related to fluid momentum
          // equation at integration point
          CalcTau(tau[k], visc, 0.0, densnp[k], scatravarmanager_->ConVel(k), vol);
          // calculation of residual-based subgrid-scale velocity
          CalcSubgrVelocity(
              ele, sgvelint, densam[k], densnp[k], visc, scatravarmanager_->ConVel(k), tau[k]);

          // calculation of subgrid-scale convective part
          sgconv.MultiplyTN(derxy_, sgvelint);
        }

        // (re)compute stabilization parameter at integration point, since diffusion may have
        // changed
        CalcTau(tau[k], diffmanager_->GetIsotropicDiff(k),
            reamanager_->GetStabilizationCoeff(k, scatravarmanager_->Phinp(k)), densnp[k],
            scatravarmanager_->ConVel(k), vol);  // TODO:(Thon) do we really have to do this??
      }

      CORE::LINALG::Matrix<nen_, 1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
      }

      // prepare multifractal subgrid-scale modeling
      // calculation of model coefficients B (velocity) and D (scalar)
      // at Gauss point as well as calculation
      // of multifractal subgrid-scale quantities
      CORE::LINALG::Matrix<nsd_, 1> mfsgvelint(true);
      double mfsvdiv(0.0);
      double mfssgphi(0.0);
      CORE::LINALG::Matrix<nsd_, 1> mfsggradphi(true);
      if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
      {
        if (turbparams_->BD_Gp())
        {
          // calculate model coefficients
          CalcBAndDForMultifracSubgridScales(B_mfs, D_mfs, vol, k, densnp[k],
              diffmanager_->GetIsotropicDiff(k), visc, scatravarmanager_->ConVel(k), fsvelint);
        }

        // calculate fine-scale velocity, its derivative and divergence for multifractal
        // subgrid-scale modeling
        for (unsigned idim = 0; idim < nsd_; idim++)
          mfsgvelint(idim, 0) = fsvelint(idim, 0) * B_mfs(idim, 0);
        // required for conservative formulation in the context of passive scalar transport
        if (turbparams_->MfsConservative() or scatrapara_->IsConservative())
        {
          // get divergence of subgrid-scale velocity
          CORE::LINALG::Matrix<nsd_, nsd_> mfsvderxy;
          mfsvderxy.MultiplyNT(efsvel_, derxy_);
          for (unsigned idim = 0; idim < nsd_; idim++)
            mfsvdiv += mfsvderxy(idim, idim) * B_mfs(idim, 0);
        }

        // calculate fine-scale scalar and its derivative for multifractal subgrid-scale modeling
        mfssgphi = D_mfs * funct_.Dot(fsphinp_[k]);
        mfsggradphi.Multiply(derxy_, fsphinp_[k]);
        mfsggradphi.Scale(D_mfs);
      }

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      // stabilization parameter and integration factors
      const double taufac = tau[k] * fac;
      const double timefacfac = scatraparatimint_->TimeFac() * fac;
      const double timetaufac = scatraparatimint_->TimeFac() * taufac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of convective element matrix in convective form
      CalcMatConv(emat, k, timefacfac, densnp[k], sgconv);

      // add conservative contributions
      if (scatrapara_->IsConservative()) CalcMatConvAddCons(emat, k, timefacfac, vdiv, densnp[k]);

      // calculation of diffusive element matrix
      CalcMatDiff(emat, k, timefacfac);

      //----------------------------------------------------------------
      // convective stabilization term
      //----------------------------------------------------------------

      // convective stabilization of convective term (in convective form)
      // transient stabilization of convective term (in convective form)
      if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
        CalcMatTransConvDiffStab(emat, k, timetaufac, densnp[k], sgconv, diff);

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------

      if (not scatraparatimint_->IsStationary())
      {
        CalcMatMass(emat, k, fac, densam[k]);

        if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
          CalcMatMassStab(emat, k, taufac, densam[k], densnp[k], sgconv, diff);
      }

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      // including stabilization
      if (reamanager_->Active())
        CalcMatReact(emat, k, timefacfac, timetaufac, taufac, densnp[k], sgconv, diff);

      //----------------------------------------------------------------
      // 4) element matrix: chemotactic term
      //----------------------------------------------------------------

      // including stabilization
      CalcMatChemo(emat, k, timefacfac, timetaufac, densnp[k], scatrares, sgconv, diff);

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac = scatraparatimint_->TimeFacRhs() * fac;
      double rhstaufac = scatraparatimint_->TimeFacRhsTau() * taufac;

      if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsStationary())
        CalcRHSLinMass(erhs, k, rhsfac, fac, densam[k], densnp[k]);

      // the order of the following three functions is important
      // and must not be changed
      ComputeRhsInt(rhsint, densam[k], densnp[k], scatravarmanager_->Hist(k));

      RecomputeScatraResForRhs(scatrares, k, diff, densn[k], densnp[k], rea_phi, rhsint);

      RecomputeConvPhiForRhs(k, sgvelint, densnp[k], densn[k], vdiv);

      //----------------------------------------------------------------
      // standard Galerkin transient, old part of rhs and bodyforce term
      //----------------------------------------------------------------
      CalcRHSHistAndSource(erhs, k, fac, rhsint);

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // convective term
      CalcRHSConv(erhs, k, rhsfac);

      // diffusive term
      CalcRHSDiff(erhs, k, rhsfac);

      //----------------------------------------------------------------
      // stabilization terms
      //----------------------------------------------------------------
      if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
        CalcRHSTransConvDiffStab(erhs, k, rhstaufac, densnp[k], scatrares, sgconv, diff);

      //----------------------------------------------------------------
      // reactive terms (standard Galerkin and stabilization) on rhs
      //----------------------------------------------------------------

      if (reamanager_->Active())
        CalcRHSReact(erhs, k, rhsfac, rhstaufac, rea_phi, densnp[k], scatrares);

      //----------------------------------------------------------------
      // chemotactic terms (standard Galerkin and stabilization) on rhs
      //----------------------------------------------------------------

      CalcRHSChemo(erhs, k, rhsfac, rhstaufac, scatrares, densnp[k]);

      //----------------------------------------------------------------
      // 6) advanced turbulence models
      //----------------------------------------------------------------

      //----------------------------------------------------------------
      // fine-scale subgrid-diffusivity term on right hand side
      //----------------------------------------------------------------
      if (scatraparatimint_->IsIncremental() and turbparams_->FSSGD())
        CalcRHSFSSGD(erhs, k, rhsfac, sgdiff, fsgradphi);

      //---------------------------------------------------------------
      // multifractal subgrid-scale modeling on right hand side only
      //---------------------------------------------------------------
      if (scatraparatimint_->IsIncremental() and
          turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales)
        CalcRHSMFS(erhs, k, rhsfac, densnp[k], mfsggradphi, mfsgvelint, mfssgphi, mfsvdiv);

      //----------------------------------------------------------------
      // 7) macro-scale matrix and vector contributions arising from
      //    macro-micro coupling in multi-scale simulations
      //----------------------------------------------------------------
      if (ele->Material()->MaterialType() == INPAR::MAT::m_scatra_multiscale)
        CalcMatAndRhsMultiScale(ele, emat, erhs, k, iquad, timefacfac, rhsfac);

      //----------------------------------------------------------------
      // 8) Compute Rhs for ElectroMagnetic Diffusion equation
      // the term includes the divergence og the electric current
      //----------------------------------------------------------------
      if (scatrapara_->IsEMD()) CalcRHSEMD(ele, erhs, rhsfac);
    }  // end loop all scalars
  }    // end loop Gauss points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::BodyForce(const DRT::Element* ele)
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  switch (nsd_ele_)
  {
    case 3:
      DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
      break;
    case 2:
      DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
      break;
    case 1:
      DRT::UTILS::FindElementConditions(ele, "LineNeumann", myneumcond);
      break;
    default:
      dserror("Illegal number of spatial dimensions: %d", nsd_ele_);
      break;
  }

  if (myneumcond.size() > 1) dserror("More than one Neumann condition on one node!");

  if (myneumcond.size() == 1)
  {
    // (SPATIAL) FUNCTION BUSINESS
    const auto* funct = myneumcond[0]->Get<std::vector<int>>("funct");

    // get values and switches from the condition
    const auto* onoff = myneumcond[0]->Get<std::vector<int>>("onoff");
    const auto* val = myneumcond[0]->Get<std::vector<double>>("val");

    // set this condition to the bodyforce array
    for (int idof = 0; idof < numdofpernode_; idof++)
    {
      // function evaluation
      const int functnum = (funct) ? (*funct)[idof] : -1;
      for (unsigned jnode = 0; jnode < nen_; jnode++)
      {
        const double functfac =
            (functnum > 0)
                ? GLOBAL::Problem::Instance()
                      ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                      .Evaluate((ele->Nodes()[jnode])->X().data(), scatraparatimint_->Time(), idof)
                : 1.0;
        (bodyforce_[idof])(jnode) = (*onoff)[idof] * (*val)[idof] * functfac;
      }
    }
  }
  else
  {
    for (int idof = 0; idof < numdofpernode_; idof++)
    {
      // no bodyforce
      bodyforce_[idof].Clear();
    }
  }
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::OtherNodeBasedSourceTerms(
    const std::vector<int>& lm, DRT::Discretization& discretization, Teuchos::ParameterList& params)
{
  // set externally calculated source term instead of body force by volume
  // Neumann boundary condition of input file
  if (turbparams_->ScalarForcing() == INPAR::FLUID::scalarforcing_isotropic)
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> source = discretization.GetState("forcing");
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*source, bodyforce_, lm);
  }
  // special forcing mean scalar gradient
  else if (turbparams_->ScalarForcing() == INPAR::FLUID::scalarforcing_mean_scalar_gradient)
  {
    // get mean-scalar gradient
    const double grad_phi = params.sublist("TURBULENCE MODEL").get<double>("MEAN_SCALAR_GRADIENT");

    // fill element array
    for (unsigned i = 0; i < nen_; ++i)
    {
      for (int k = 0; k < numdofpernode_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        bodyforce_[k](i, 0) = -grad_phi * evelnp_(2, i);
      }
    }  // for i
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ReadElementCoordinates(const DRT::Element* ele)
{
  // Directly copy the coordinates since in 3D the transformation is just the identity
  CORE::GEO::fillInitialPositionArray<distype, nsd_, CORE::LINALG::Matrix<nsd_, nen_>>(ele, xyze_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvalShapeFuncAndDerivsAtEleCenter()
{
  // use one-point Gauss rule to do calculations at the element center
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

  return vol;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvalShapeFuncAndDerivsAtIntPoint(
    const CORE::FE::IntPointsAndWeights<nsd_ele_>& intpoints, const int iquad)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (unsigned idim = 0; idim < nsd_ele_; idim++) xsi_(idim) = gpcoord[idim];

  const double det = EvalShapeFuncAndDerivsInParameterSpace();

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

  // compute global spatial derivatives
  derxy_.Multiply(xij_, deriv_);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad] * det;

  // compute second global derivatives (if needed)
  if (use2ndderiv_)
  {
    // get global second derivatives
    CORE::FE::gder2<distype, nen_, probdim>(xjm_, derxy_, deriv2_, xyze_, derxy2_);
  }
  else
    derxy2_.Clear();

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::EvalShapeFuncAndDerivsInParameterSpace()
{
  double det = 0.0;

  if (nsd_ == nsd_ele_)  // standard case
  {
    if (not DRT::NURBS::IsNurbs(distype))
    {
      // shape functions and their first derivatives
      CORE::FE::shape_function<distype>(xsi_, funct_);
      CORE::FE::shape_function_deriv1<distype>(xsi_, deriv_);
      if (use2ndderiv_)
      {
        // get the second derivatives of standard element at current GP
        CORE::FE::shape_function_deriv2<distype>(xsi_, deriv2_);
      }
    }
    else  // nurbs elements are always somewhat special...
    {
      if (use2ndderiv_)
      {
        CORE::FE::NURBS::nurbs_get_funct_deriv_deriv2(
            funct_, deriv_, deriv2_, xsi_, myknots_, weights_, distype);
      }
      else
      {
        CORE::FE::NURBS::nurbs_get_funct_deriv(funct_, deriv_, xsi_, myknots_, weights_, distype);
      }
    }  // IsNurbs()

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

    xjm_.MultiplyNT(deriv_, xyze_);
    det = xij_.Invert(xjm_);
  }
  else  // element dimension is smaller than problem dimension -> manifold
  {
    static CORE::LINALG::Matrix<nsd_ele_, nen_> deriv_red;

    if (not DRT::NURBS::IsNurbs(distype))
    {
      // shape functions and their first derivatives
      CORE::FE::shape_function<distype>(xsi_, funct_);
      CORE::FE::shape_function_deriv1<distype>(xsi_, deriv_red);
      if (use2ndderiv_)
      {
        // get the second derivatives of standard element at current GP
        CORE::FE::shape_function_deriv2<distype>(xsi_, deriv2_);
      }
    }
    else  // nurbs elements are always somewhat special...
    {
      if (use2ndderiv_)
      {
        CORE::FE::NURBS::nurbs_get_funct_deriv_deriv2(
            funct_, deriv_red, deriv2_, xsi_, myknots_, weights_, distype);
      }
      else
      {
        CORE::FE::NURBS::nurbs_get_funct_deriv(
            funct_, deriv_red, xsi_, myknots_, weights_, distype);
      }
    }  // IsNurbs()

    //! metric tensor at integration point
    static CORE::LINALG::Matrix<nsd_ele_, nsd_ele_> metrictensor;
    static CORE::LINALG::Matrix<nsd_, 1> normalvec;

    // the metric tensor and the area of an infinitesimal surface/line element
    // optional: get unit normal at integration point as well
    const bool throw_error_if_negative_determinant(true);
    CORE::FE::ComputeMetricTensorForBoundaryEle<distype, nsd_>(
        xyze_, deriv_red, metrictensor, det, throw_error_if_negative_determinant, &normalvec);

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

    // transform the derivatives and Jacobians to the higher dimensional coordinates(problem
    // dimension)
    static CORE::LINALG::Matrix<nsd_ele_, nsd_> xjm_red;
    xjm_red.MultiplyNT(deriv_red, xyze_);

    for (unsigned i = 0; i < nsd_; ++i)
    {
      for (unsigned j = 0; j < nsd_ele_; ++j) xjm_(j, i) = xjm_red(j, i);
      xjm_(nsd_ele_, i) = normalvec(i, 0);
    }

    for (unsigned i = 0; i < nen_; ++i)
    {
      for (unsigned j = 0; j < nsd_ele_; ++j) deriv_(j, i) = deriv_red(j, i);
      deriv_(nsd_ele_, i) = 0.0;
    }

    // special case: 1D element embedded in 3D problem
    if (nsd_ele_ == 1 and nsd_ == 3)
    {
      // compute second unit normal
      const double normalvec2_0 = xjm_red(0, 1) * normalvec(2, 0) - normalvec(1, 0) * xjm_red(0, 2);
      const double normalvec2_1 = xjm_red(0, 2) * normalvec(0, 0) - normalvec(2, 0) * xjm_red(0, 0);
      const double normalvec2_2 = xjm_red(0, 0) * normalvec(1, 0) - normalvec(0, 0) * xjm_red(0, 1);

      // norm
      const double norm2 = std::sqrt(
          normalvec2_0 * normalvec2_0 + normalvec2_1 * normalvec2_1 + normalvec2_2 * normalvec2_2);

      xjm_(2, 0) = normalvec2_0 / norm2;
      xjm_(2, 1) = normalvec2_1 / norm2;
      xjm_(2, 2) = normalvec2_2 / norm2;

      for (unsigned i = 0; i < nen_; i++) deriv_(2, i) = 0.0;
    }

    xij_.Invert(xjm_);
  }

  // modify Jacobian determinant in case of spherical coordinates
  if (scatrapara_->SphericalCoords())
  {
    static CORE::LINALG::Matrix<nsd_, 1> xyzint;

    // evaluate radial coordinate
    xyzint.Multiply(xyze_, funct_);

    // multiply standard Jacobian determinant by square of radial coordinate and 4 pi
    constexpr double four_pi = 4.0 * M_PI;
    det *= xyzint(0) * xyzint(0) * four_pi;
  }

  return det;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::GetMaterialParams(const DRT::Element* ele,
    std::vector<double>& densn, std::vector<double>& densnp, std::vector<double>& densam,
    double& visc, const int iquad)
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0; k < numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat, k, densn[k], densnp[k], densam[k], visc, iquad);
    }
  }
  else
    Materials(material, 0, densn[0], densnp[0], densam[0], visc, iquad);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::Materials(
    const Teuchos::RCP<const MAT::Material> material, const int k, double& densn, double& densnp,
    double& densam, double& visc, const int iquad)
{
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_electrode:
    {
      // safety check
      if (k != 0) dserror("Invalid species ID!");

      MatElectrode(material);
      break;
    }

    case INPAR::MAT::m_scatra:
    {
      MatScaTra(material, k, densn, densnp, densam, visc, iquad);
      break;
    }

    case INPAR::MAT::m_scatra_multiscale:
    {
      MatScaTraMultiScale(material, densn, densnp, densam);
      break;
    }

    default:
    {
      dserror("Material type %i is not supported!", material->MaterialType());
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::MatScaTra(
    const Teuchos::RCP<const MAT::Material> material, const int k, double& densn, double& densnp,
    double& densam, double& visc, const int iquad)
{
  const Teuchos::RCP<const MAT::ScatraMat>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  // get constant diffusivity
  diffmanager_->SetIsotropicDiff(actmat->Diffusivity(), k);

  // get reaction coefficient
  reamanager_->SetReaCoeff(actmat->ReaCoeff(), k);

  // in case of multifractal subgrid-scales, read Schmidt number
  if (turbparams_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales or
      scatrapara_->RBSubGrVel() or turbparams_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
  {
    // access fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
    fluiddis = GLOBAL::Problem::Instance()->GetDis("fluid");
    // get corresponding fluid element (it has the same global ID as the scatra element)
    DRT::Element* fluidele = fluiddis->gElement(eid_);
    if (fluidele == nullptr)
      dserror("Fluid element %i not on local processor", eid_);
    else
    {
      // get fluid material
      Teuchos::RCP<MAT::Material> fluidmat = fluidele->Material();
      if (fluidmat->MaterialType() != INPAR::MAT::m_fluid)
        dserror("Invalid fluid material for passive scalar transport in turbulent flow!");

      const Teuchos::RCP<const MAT::NewtonianFluid>& actfluidmat =
          Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat);

      // get constant dynamic viscosity
      visc = actfluidmat->Viscosity();
      densn = actfluidmat->Density();
      densnp = actfluidmat->Density();
      densam = actfluidmat->Density();

      if (densam != 1.0 or densnp != 1.0 or densn != 1.0)
        dserror("Check your parameters! Read comment!");
      // For all implementations, dens=1.0 is assumed, in particular for
      // multifractal_subgrid_scales. Hence, visc and diffus are kinematic quantities. Using
      // dens!=1.0 should basically work, but you should check it before application.
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::MatScaTraMultiScale(
    const Teuchos::RCP<const MAT::Material> material, double& densn, double& densnp,
    double& densam) const
{
  // safety check
  if (numscal_ > 1)
    dserror("Multi-scale scalar transport only implemented for one transported scalar!");

  // extract multi-scale scalar transport material
  const auto* matmultiscale = static_cast<const MAT::ScatraMultiScale*>(material.get());

  // set densities equal to porosity
  densn = densnp = densam = matmultiscale->Porosity();

  // set effective diffusion coefficient in diffusion manager
  // effective diffusion coefficient = intrinsic diffusion coefficient * porosity / tortuosity
  diffmanager_->SetIsotropicDiff(
      matmultiscale->Diffusivity() * matmultiscale->Porosity() / matmultiscale->Tortuosity(), 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::MatElectrode(
    const Teuchos::RCP<const MAT::Material> material)
{
  // set constant diffusivity
  diffmanager_->SetIsotropicDiff(
      Teuchos::rcp_static_cast<const MAT::Electrode>(material)
          ->ComputeDiffusionCoefficientConcentrationDependent(scatravarmanager_->Phinp(0)),
      0);
}

/*---------------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::GetLaplacianStrongForm(
    CORE::LINALG::Matrix<nen_, 1>& diff)
{
  diff.Clear();
  // compute N,xx  +  N,yy +  N,zz for each shape function at integration point
  for (unsigned i = 0; i < nen_; ++i)
  {
    for (unsigned j = 0; j < nsd_; ++j)
    {
      diff(i) += derxy2_(j, i);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::GetDivergence(
    double& vdiv, const CORE::LINALG::Matrix<nsd_, nen_>& evel)
{
  CORE::LINALG::Matrix<nsd_, nsd_> vderxy;
  vderxy.MultiplyNT(evel, derxy_);

  vdiv = 0.0;
  // compute vel x,x  + vel y,y +  vel z,z at integration point
  for (unsigned j = 0; j < nsd_; ++j)
  {
    vdiv += vderxy(j, j);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::GetRhsInt(
    double& rhsint, const double densnp, const int k)
{
  // compute rhs containing bodyforce (divided by specific heat capacity) and,
  // for temperature equation, the time derivative of thermodynamic pressure,
  // if not constant, and for temperature equation of a reactive
  // equation system, the reaction-rate term
  rhsint = bodyforce_[k].Dot(funct_);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatConv(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double densnp, const CORE::LINALG::Matrix<nen_, 1>& sgconv)
{
  const CORE::LINALG::Matrix<nen_, 1>& conv = scatravarmanager_->Conv(k);

  // convective term in convective form
  const double densfac = timefacfac * densnp;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = densfac * funct_(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * (conv(ui) + sgconv(ui));
    }
  }
}

/*------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatConvAddCons(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac, const double vdiv,
    const double densnp)
{
  const double consfac = timefacfac * densnp * vdiv;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = consfac * funct_(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * funct_(ui);
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatDiff(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac)
{
  // diffusive term
  const double fac_diffus = timefacfac * diffmanager_->GetIsotropicDiff(k);
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf, ui, vi);
      emat(fvi, fui) += fac_diffus * laplawf;
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatTransConvDiffStab(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timetaufac,
    const double densnp, const CORE::LINALG::Matrix<nen_, 1>& sgconv,
    const CORE::LINALG::Matrix<nen_, 1>& diff)
{
  const CORE::LINALG::Matrix<nen_, 1>& conv = scatravarmanager_->Conv(k);

  const double dens2taufac = timetaufac * densnp * densnp;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = dens2taufac * (conv(vi) + sgconv(vi) +
                                       scatrapara_->USFEMGLSFac() * 1.0 /
                                           scatraparatimint_->TimeFac() * funct_(vi));
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * conv(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization terms for higher-order elements
  //----------------------------------------------------------------
  if (use2ndderiv_)
  {
    const double denstaufac = timetaufac * densnp;
    // convective stabilization of diffusive term (in convective form)
    // transient stabilization of diffusive term (in convective form)
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = denstaufac * (conv(vi) + sgconv(vi) +
                                        scatrapara_->USFEMGLSFac() * 1.0 /
                                            scatraparatimint_->TimeFac() * funct_(vi));
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) -= v * diff(ui);
      }
    }

    const double densdifftaufac = scatrapara_->USFEMGLSFac() * denstaufac;
    // diffusive stabilization of convective term (in convective form)
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = densdifftaufac * diff(vi);
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) -= v * conv(ui);
      }
    }

    const double difftaufac = scatrapara_->USFEMGLSFac() * timetaufac;
    // diffusive stabilization of diffusive term
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = difftaufac * diff(vi);
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) += v * diff(ui);
      }
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatMass(
    CORE::LINALG::SerialDenseMatrix& emat, const int& k, const double& fac, const double& densam)
{
  CalcMatMass(emat, k, fac, densam, funct_, funct_);
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatMass(
    CORE::LINALG::SerialDenseMatrix& emat, const int& k, const double& fac, const double& densam,
    const CORE::LINALG::Matrix<nen_, 1>& sfunct, const CORE::LINALG::Matrix<nen_, 1>& tfunct) const
{
  const double densamfac = fac * densam;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = densamfac * tfunct(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * sfunct(ui);
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatMassStab(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double taufac, const double densam,
    const double densnp, const CORE::LINALG::Matrix<nen_, 1>& sgconv,
    const CORE::LINALG::Matrix<nen_, 1>& diff)
{
  const CORE::LINALG::Matrix<nen_, 1>& conv = scatravarmanager_->Conv(k);
  const double densamnptaufac = taufac * densam * densnp;
  //----------------------------------------------------------------
  // stabilization of transient term
  //----------------------------------------------------------------
  // convective stabilization of transient term (in convective form)
  // transient stabilization of transient term
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = densamnptaufac * (conv(vi) + sgconv(vi) +
                                          scatrapara_->USFEMGLSFac() * 1.0 /
                                              scatraparatimint_->TimeFac() * funct_(vi));
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * funct_(ui);
    }
  }

  if (use2ndderiv_)
  {
    const double densamreataufac = scatrapara_->USFEMGLSFac() * taufac * densam;
    // diffusive stabilization of transient term
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = densamreataufac * diff(vi);
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) -= v * funct_(ui);
      }
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatReact(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const double timetaufac, const double taufac, const double densnp,
    const CORE::LINALG::Matrix<nen_, 1>& sgconv, const CORE::LINALG::Matrix<nen_, 1>& diff)
{
  const CORE::LINALG::Matrix<nen_, 1>& conv = scatravarmanager_->Conv(k);

  // NOTE: it is important that the reaction coefficient reamanager_->GetReaCoeff(k) does not depend
  // on ANY concentrations.
  const double fac_reac = timefacfac * densnp * reamanager_->GetReaCoeff(k);
  const double timetaufac_reac = timetaufac * densnp * reamanager_->GetReaCoeff(k);

  //----------------------------------------------------------------
  // standard Galerkin reactive term
  //----------------------------------------------------------------
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = fac_reac * funct_(vi);
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * numdofpernode_ + k;

      emat(fvi, fui) += v * funct_(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization of reactive term
  //----------------------------------------------------------------
  if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
  {
    double densreataufac = timetaufac_reac * densnp;
    // convective stabilization of reactive term (in convective form)
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = densreataufac * (conv(vi) + sgconv(vi) +
                                           scatrapara_->USFEMGLSFac() * 1.0 /
                                               scatraparatimint_->TimeFac() * funct_(vi));
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) += v * funct_(ui);
      }
    }

    if (use2ndderiv_)
    {
      // diffusive stabilization of reactive term
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac() * timetaufac_reac * diff(vi);
        const int fvi = vi * numdofpernode_ + k;

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const int fui = ui * numdofpernode_ + k;

          emat(fvi, fui) -= v * funct_(ui);
        }
      }
    }

    //----------------------------------------------------------------
    // reactive stabilization
    //----------------------------------------------------------------
    densreataufac = scatrapara_->USFEMGLSFac() * timetaufac_reac * densnp;

    // reactive stabilization of convective (in convective form) and reactive term
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double v = densreataufac * funct_(vi);
      const int fvi = vi * numdofpernode_ + k;

      for (unsigned ui = 0; ui < nen_; ++ui)
      {
        const int fui = ui * numdofpernode_ + k;

        emat(fvi, fui) += v * (conv(ui) + reamanager_->GetReaCoeff(k) * funct_(ui));
      }
    }

    if (use2ndderiv_)
    {
      // reactive stabilization of diffusive term
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac() * timetaufac_reac * funct_(vi);
        const int fvi = vi * numdofpernode_ + k;

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const int fui = ui * numdofpernode_ + k;

          emat(fvi, fui) -= v * diff(ui);
        }
      }
    }


    if (not scatraparatimint_->IsStationary())
    {
      // reactive stabilization of transient term
      for (unsigned vi = 0; vi < nen_; ++vi)
      {
        const double v = scatrapara_->USFEMGLSFac() * taufac * densnp *
                         reamanager_->GetReaCoeff(k) * densnp * funct_(vi);
        const int fvi = vi * numdofpernode_ + k;

        for (unsigned ui = 0; ui < nen_; ++ui)
        {
          const int fui = ui * numdofpernode_ + k;

          emat(fvi, fui) += v * funct_(ui);
        }
      }

      if (use2ndderiv_ and reamanager_->GetReaCoeff(k) != 0.0)
        dserror("Second order reactive stabilization is not fully implemented!! ");
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSLinMass(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double fac,
    const double densam, const double densnp)
{
  const double& phinp = scatravarmanager_->Phinp(k);
  const double& hist = scatravarmanager_->Hist(k);

  double vtrans = 0.0;

  if (scatraparatimint_->IsGenAlpha())
    vtrans = rhsfac * densam * hist;
  else
  {
    // compute scalar at integration point
    vtrans = fac * densnp * phinp;
  }

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    erhs[fvi] -= vtrans * funct_(vi);
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ComputeRhsInt(
    double& rhsint, const double densam, const double densnp, const double hist)
{
  if (scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsIncremental())
      rhsint += densam * hist * (scatraparatimint_->AlphaF() / scatraparatimint_->TimeFac());

    rhsint *= (scatraparatimint_->TimeFac() / scatraparatimint_->AlphaF());
  }
  else  // OST, BDF2, stationary
  {
    if (not scatraparatimint_->IsStationary())  // OST, BDF2
    {
      rhsint *= scatraparatimint_->TimeFac();
      rhsint += densnp * hist;
    }
  }
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::RecomputeScatraResForRhs(double& scatrares,
    const int k, const CORE::LINALG::Matrix<nen_, 1>& diff, const double densn, const double densnp,
    double& rea_phi, const double rhsint)
{
  const CORE::LINALG::Matrix<nsd_, 1>& convelint = scatravarmanager_->ConVel(k);
  const double& phin = scatravarmanager_->Phin(k);

  if (scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsIncremental())
    {
      // for this case, gradphi_ (i.e. the gradient
      // at time n+1) is overwritten by the gradient at time n
      // analogously, conv_phi_ at time n+1 is replace by its
      // value at time n
      // gradient of scalar value at n
      CORE::LINALG::Matrix<nsd_, 1> gradphi;
      gradphi.Multiply(derxy_, ephin_[k]);

      // convective term using scalar value at n
      double conv_phi = convelint.Dot(gradphi);

      // diffusive term using current scalar value for higher-order elements
      double diff_phin = 0.0;
      if (use2ndderiv_) diff_phin = diff.Dot(ephin_[k]);

      // reactive term using scalar value at n
      // if no reaction is chosen, GetReaCoeff(k) returns 0.0
      rea_phi = densnp * reamanager_->GetReaCoeff(k) * phin;
      // reacterm_[k] must be evaluated at t^n to be used in the line above!

      scatrares = (1.0 - scatraparatimint_->AlphaF()) * (densn * conv_phi - diff_phin + rea_phi) -
                  rhsint * scatraparatimint_->AlphaF() / scatraparatimint_->TimeFac();

      scatravarmanager_->SetGradPhi(k, gradphi);
      scatravarmanager_->SetConvPhi(k, conv_phi);
    }
  }
  else if (scatraparatimint_->IsIncremental() and not scatraparatimint_->IsGenAlpha())
  {
    if (not scatraparatimint_->IsStationary()) scatrares *= scatraparatimint_->Dt();
  }
  else
    scatrares = -rhsint;
}

/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::RecomputeConvPhiForRhs(const int k,
    const CORE::LINALG::Matrix<nsd_, 1>& sgvelint, const double densnp, const double densn,
    const double vdiv)
{
  double conv_phi = 0.0;
  const double& phinp = scatravarmanager_->Phinp(k);
  const double& phin = scatravarmanager_->Phin(k);
  const CORE::LINALG::Matrix<nsd_, 1>& gradphi = scatravarmanager_->GradPhi(k);


  if (scatraparatimint_->IsIncremental())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (scatrapara_->IsConservative())
    {
      // convective term in conservative form
      conv_phi += phinp * vdiv;
    }

    // multiply convective term by density
    conv_phi *= densnp;

    // multiply convective term by density
    scatravarmanager_->ScaleConvPhi(k, densnp);
    // addition to convective term due to subgrid-scale velocity
    scatravarmanager_->AddToConvPhi(k, conv_phi);
  }
  else if (not scatraparatimint_->IsIncremental() and scatraparatimint_->IsGenAlpha())
  {
    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint.Dot(gradphi);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (scatrapara_->IsConservative())
    {
      // convective term in conservative form
      // caution: velocity divergence is for n+1 and not for n!
      // -> hopefully, this inconsistency is of small amount
      conv_phi += phin * vdiv;
    }

    // multiply convective term by density
    conv_phi *= densn;

    // multiply convective term by density
    scatravarmanager_->ScaleConvPhi(k, densn);
    // addition to convective term due to subgrid-scale velocity
    scatravarmanager_->AddToConvPhi(k, conv_phi);
  }
}

/*-------------------------------------------------------------------------------------- *
 *---------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSHistAndSource(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double fac, const double rhsint)
{
  double vrhs = fac * rhsint;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    erhs[fvi] += vrhs * funct_(vi);
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSConv(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac)
{
  const double& conv_phi = scatravarmanager_->ConvPhi(k);

  double vrhs = rhsfac * conv_phi;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    erhs[fvi] -= vrhs * funct_(vi);
  }
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSDiff(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac)
{
  const CORE::LINALG::Matrix<nsd_, 1>& gradphi = scatravarmanager_->GradPhi(k);

  double vrhs = rhsfac * diffmanager_->GetIsotropicDiff(k);

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf, gradphi, vi);
    erhs[fvi] -= vrhs * laplawf;
  }
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSTransConvDiffStab(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhstaufac, const double densnp,
    const double scatrares, const CORE::LINALG::Matrix<nen_, 1>& sgconv,
    const CORE::LINALG::Matrix<nen_, 1>& diff)
{
  const CORE::LINALG::Matrix<nen_, 1>& conv = scatravarmanager_->Conv(k);

  // convective rhs stabilization (in convective form)
  double vrhs = rhstaufac * scatrares * densnp;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    erhs[fvi] -=
        vrhs * (conv(vi) + sgconv(vi) +
                   scatrapara_->USFEMGLSFac() * 1.0 / scatraparatimint_->TimeFac() * funct_(vi));
  }

  // diffusive rhs stabilization
  if (use2ndderiv_)
  {
    vrhs = rhstaufac * scatrares;
    // diffusive stabilization of convective temporal rhs term (in convective form)
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;

      erhs[fvi] += scatrapara_->USFEMGLSFac() * vrhs * diff(vi);
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSReact(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double rhstaufac,
    const double rea_phi, const double densnp, const double scatrares)
{
  // standard Galerkin term
  double vrhs = rhsfac * rea_phi;

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    erhs[fvi] -= vrhs * funct_(vi);
  }

  // reactive rhs stabilization
  if (scatrapara_->StabType() != INPAR::SCATRA::stabtype_no_stabilization)
  {
    vrhs =
        scatrapara_->USFEMGLSFac() * rhstaufac * densnp * reamanager_->GetReaCoeff(k) * scatrares;
    // TODO: this is not totally correct since GetReaCoeff(k) can depend on phinp(k)...
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;

      erhs[fvi] -= vrhs * funct_(vi);
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSFSSGD(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double sgdiff,
    const CORE::LINALG::Matrix<nsd_, 1> fsgradphi)
{
  const double vrhs = rhsfac * sgdiff;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * numdofpernode_ + k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf, fsgradphi, vi);
    erhs[fvi] -= (vrhs * laplawf);
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSMFS(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac, const double densnp,
    const CORE::LINALG::Matrix<nsd_, 1> mfsggradphi, const CORE::LINALG::Matrix<nsd_, 1> mfsgvelint,
    const double mfssgphi, const double mfsvdiv)
{
  const double& phinp = scatravarmanager_->Phinp(k);
  const CORE::LINALG::Matrix<nsd_, 1>& gradphi = scatravarmanager_->GradPhi(k);
  const CORE::LINALG::Matrix<nsd_, 1>& convelint = scatravarmanager_->ConVel(k);

  if (nsd_ < 3) dserror("Turbulence is 3D!");
  // fixed-point iteration only (i.e. beta=0.0 assumed), cf
  // turbulence part in Evaluate()
  {
    double cross = convelint.Dot(mfsggradphi) + mfsgvelint.Dot(gradphi);
    double reynolds = mfsgvelint.Dot(mfsggradphi);

    // conservative formulation
    double conserv = 0.0;
    if (turbparams_->MfsConservative() or scatrapara_->IsConservative())
    {
      double convdiv = 0.0;
      GetDivergence(convdiv, econvelnp_);

      conserv = mfssgphi * convdiv + phinp * mfsvdiv + mfssgphi * mfsvdiv;
    }

    const double vrhs = rhsfac * densnp * (cross + reynolds + conserv);

    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const int fvi = vi * numdofpernode_ + k;
      erhs[fvi] -= vrhs * funct_(vi);
    }
  }
}

/*-----------------------------------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcMatAndRhsMultiScale(
    const DRT::Element* const ele, CORE::LINALG::SerialDenseMatrix& emat,
    CORE::LINALG::SerialDenseVector& erhs, const int k, const int iquad, const double timefacfac,
    const double rhsfac)
{
  // extract multi-scale scalar transport material
  const MAT::ScatraMultiScale* matmultiscale =
      static_cast<const MAT::ScatraMultiScale*>(ele->Material().get());

  // initialize variables for micro-scale coupling flux and derivative of micro-scale coupling flux
  // w.r.t. macro-scale state variable
  double q_micro(0.0);
  std::vector<double> dq_dphi_micro(1, 0.0);

  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const double detF = EvalDetFAtIntPoint(ele, intpoints, iquad);

  // evaluate multi-scale scalar transport material
  matmultiscale->Evaluate(
      iquad, std::vector<double>(1, scatravarmanager_->Phinp(k)), q_micro, dq_dphi_micro, detF);

  // macro-scale matrix contribution
  const double matrixterm =
      timefacfac * dq_dphi_micro[0] * matmultiscale->SpecificMicroScaleSurfaceArea(detF);
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const double v = funct_(vi) * matrixterm;
    const int fvi = vi * numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui) emat(fvi, ui * numdofpernode_ + k) += v * funct_(ui);
  }

  // macro-scale vector contribution
  const double rhsterm = rhsfac * q_micro * matmultiscale->SpecificMicroScaleSurfaceArea(detF);
  for (unsigned vi = 0; vi < nen_; ++vi) erhs[vi * numdofpernode_ + k] -= funct_(vi) * rhsterm;
}

/*-------------------------------------------------------------------- *
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::CalcRHSEMD(
    const DRT::Element* const ele, CORE::LINALG::SerialDenseVector& erhs, const double rhsfac)
{
  // (SPATIAL) FUNCTION BUSINESS
  const int functno = scatrapara_->EMDSource();
  if (functno <= 0)
  {
    dserror(
        "For electromagnetic diffusion simulations a current density source function has to be "
        "given.");
  }

  std::vector<double> current(nsd_, 0);
  for (unsigned jnode = 0; jnode < nen_; jnode++)
  {
    for (int d = 0; d < static_cast<int>(nsd_); ++d)
    {
      current[d] += funct_(jnode) *
                    GLOBAL::Problem::Instance()
                        ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functno - 1)
                        .Evaluate((ele->Nodes()[jnode])->X().data(), scatraparatimint_->Time(), d);
    }
  }

  for (int vi = 0; vi < static_cast<int>(nen_); ++vi)
  {
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      erhs[vi] += derxy_(d, vi) * current[d] * rhsfac;
    }
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::SetInternalVariablesForMatAndRHS()
{
  scatravarmanager_->SetInternalVariables(
      funct_, derxy_, ephinp_, ephin_, econvelnp_, ehist_, eforcevelocity_);
}

BACI_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
#include "baci_scatra_ele_calc_fwd.hpp"
