/*-----------------------------------------------------------*/
/*! \file

\brief Setting of specific XFEM based fluid parameter for element evaluation


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_ele_parameter_xfem.hpp"

#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidEleParameterXFEM* Discret::ELEMENTS::FluidEleParameterXFEM::Instance(
    Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleParameterXFEM>(
            new Discret::ELEMENTS::FluidEleParameterXFEM());
      });

  return singleton_owner.Instance(action);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleParameterXFEM::FluidEleParameterXFEM()
    : Discret::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd(),
      vcellgausspts_(Inpar::Cut::VCellGaussPts_DirectDivergence),
      bcellgausspts_(Inpar::Cut::BCellGaussPts_Tessellation),
      coupling_method_(Inpar::XFEM::Nitsche),
      hybrid_lm_l2_proj_(Inpar::XFEM::Hybrid_LM_L2_Proj_part),
      visc_stab_trace_estimate_(Inpar::XFEM::ViscStab_TraceEstimate_CT_div_by_hk),
      visc_stab_hk_(Inpar::XFEM::ViscStab_hk_vol_equivalent),
      nit_stab_gamma_(0.0),
      nit_stab_gamma_tang_(0.0),
      visc_adjoint_scaling_(Inpar::XFEM::adj_sym),
      is_pseudo_2_d_(false),
      xff_conv_stab_scaling_(Inpar::XFEM::XFF_ConvStabScaling_none),
      conv_stab_scaling_(Inpar::XFEM::ConvStabScaling_none),
      mass_conservation_combo_(Inpar::XFEM::MassConservationCombination_max),
      mass_conservation_scaling_(Inpar::XFEM::MassConservationScaling_only_visc)
{
}

//----------------------------------------------------------------------*/
//    check parameter combination for consistency
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterXFEM::check_parameter_consistency(int myrank) const
{
  if (visc_adjoint_scaling_ == Inpar::XFEM::adj_sym &&
      coupling_method_ == Inpar::XFEM::Hybrid_LM_viscous_stress)
  {
    if (myrank == 0)
      Core::IO::cout
          << "Be warned: the symmetric hybrid/viscous stress-based LM approach is known for "
             "unstable behaviour in xfluid-fluid problems."
          << Core::IO::endl;
  }

#ifdef FOUR_C_ENABLE_ASSERTIONS
  switch (intterms_prev_state_)
  {
    case Inpar::XFEM::PreviousState_only_consistency:
    {
      if (myrank == 0)
        Core::IO::cout << "Treatment of interface terms in case of new OST: ONLY CONSISTENCY \n"
                       << "only standard consistency terms at t_n!\n"
                       << "Be careful in cases of non-stationary XFEM interfaces."
                       << Core::IO::endl;
      break;
    }
    case Inpar::XFEM::PreviousState_full:
    {
      if (myrank == 0)
        Core::IO::cout << "Treatment of interface terms in case of new OST: FULL \n"
                       << "all interface terms at t_n (standard + adjoint + penalty)!\n"
                       << "Be careful in cases of non-stationary XFEM interfaces."
                       << Core::IO::endl;
      break;
    }
    default:
      FOUR_C_THROW("Treatment of interface terms for new OST not specified.");
      break;
  }
#endif
}

//----------------------------------------------------------------------*/
//    check parameter combination for consistency
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterXFEM::check_parameter_consistency_for_averaging_strategy(
    int myrank, Inpar::XFEM::AveragingStrategy averaging_strategy) const
{
  // Determine, whether this is an embedded-sided Nitsche-approach
  const bool isEmbNitsche = (coupling_method_ == Inpar::XFEM::Nitsche &&
                             averaging_strategy == Inpar::XFEM::Embedded_Sided);

  if (visc_stab_trace_estimate_ == Inpar::XFEM::ViscStab_TraceEstimate_eigenvalue && !isEmbNitsche)
    FOUR_C_THROW(
        "Solution of eigenvalue problem to estimate parameter from trace inequality is only "
        "reasonable for embedded-sided Nitsche coupling.");

  // Consistency Checks for characteristic element length definitions
  if (isEmbNitsche)
  {
    // check element l
    if (visc_stab_trac_estimate() != Inpar::XFEM::ViscStab_TraceEstimate_eigenvalue)
    {
      if (ViscStabHK() == Inpar::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf or
          ViscStabHK() == Inpar::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf)
        FOUR_C_THROW(
            "chosen characteristic element length definition ViscStabHK is not supported for "
            "embedded-sided Nitsche method");

      if (ViscStabHK() == Inpar::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf)
        FOUR_C_THROW(
            "chosen characteristic element length definition "
            "ViscStab_hk_ele_vol_div_by_max_ele_surf is supported for embedded-sided Nitsche "
            "method,"
            "however not a reasonable choice, as it can lead to an ill-conditioning due to a to "
            "large penalty scaling for anisotropic embedded meshes!"
            "If you want, you can try it!");
    }
  }

  // Consistency Check for EVP and Xfluid-sided Nitsche
  if (!isEmbNitsche and visc_stab_trac_estimate() == Inpar::XFEM::ViscStab_TraceEstimate_eigenvalue)
    FOUR_C_THROW(
        "estimating trace inequality scaling via solving eigenvalue problems not supported for "
        "xfluid-sided Nitsche!");

  if (!isEmbNitsche and ViscStabHK() == Inpar::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf)
    FOUR_C_THROW(
        "chosen characteristic element length definition ViscStabHK is not supported for "
        "xfluid-sided Nitsche method as the element surface cannot be specified for cut elements");

  if (!isEmbNitsche and (ViscStabHK() == Inpar::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf or
                            ViscStabHK() == Inpar::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf))
  {
    if (myrank == 0)
      Core::IO::cout
          << "Be warned: the chosen characteristic element length definition ViscStabHK can "
             "become critical for xfluid-sided Nitsche method as the current definition "
             "either can not guarantee a sufficient estimate of the inverse inequality or can "
             "lead to cut position dependent error behaviour!"
          << Core::IO::endl;
  }
}

//----------------------------------------------------------------------*/
//    set parameters
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterXFEM::set_element_xfem_parameter(
    Teuchos::ParameterList& params,  ///< parameter list
    int myrank                       ///< pid (for output purpose)
)
{
  //----------------------------------------------------------------
  Teuchos::ParameterList& params_xfem = params.sublist("XFEM");

  vcellgausspts_ =
      Core::UTILS::IntegralValue<Inpar::Cut::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  bcellgausspts_ = Core::UTILS::IntegralValue<Inpar::Cut::BCellGaussPts>(
      params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  //----------------------------------------------------------------
  Teuchos::ParameterList& params_xf_gen = params.sublist("XFLUID DYNAMIC/GENERAL");

  //----------------------------------------------------------------
  Teuchos::ParameterList& params_xf_stab = params.sublist("XFLUID DYNAMIC/STABILIZATION");

  //--------------------------------------------
  // parameters describing the coupling approach
  //---------------------------------------------
  coupling_method_ =
      Core::UTILS::IntegralValue<Inpar::XFEM::CouplingMethod>(params_xf_stab, "COUPLING_METHOD");

  hybrid_lm_l2_proj_ =
      Core::UTILS::IntegralValue<Inpar::XFEM::HybridLmL2Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  //--------------------------------------------
  // parameters for the viscous stabilization in Nitsche's method and MixedHybrid_LM methods
  //---------------------------------------------

  visc_stab_trace_estimate_ = Core::UTILS::IntegralValue<Inpar::XFEM::ViscStabTraceEstimate>(
      params_xf_stab, "VISC_STAB_TRACE_ESTIMATE");

  visc_stab_hk_ =
      Core::UTILS::IntegralValue<Inpar::XFEM::ViscStabHk>(params_xf_stab, "VISC_STAB_HK");

  nit_stab_gamma_ = params_xf_stab.get<double>("NIT_STAB_FAC");

  nit_stab_gamma_tang_ = params_xf_stab.get<double>("NIT_STAB_FAC_TANG");

  visc_adjoint_scaling_ = Core::UTILS::IntegralValue<Inpar::XFEM::AdjointScaling>(
      params_xf_stab, "VISC_ADJOINT_SYMMETRY");

  is_pseudo_2_d_ = (bool)Core::UTILS::IntegralValue<int>(params_xf_stab, "IS_PSEUDO_2D");


  // TODO: add a comment how to define the visc-stab-fac for eigenvalue problem
  // TODO add or adapt a comment like this when using the eigenvalue problem
  // be at least 2 to get a stable formulation. To reach the value corresponding to alpha
  // equal to 35 it should be higher than 2.

  //--------------------------------------------
  // parameters for the convective interface stabilizations
  //---------------------------------------------

  xff_conv_stab_scaling_ = Core::UTILS::IntegralValue<Inpar::XFEM::XffConvStabScaling>(
      params_xf_stab, "XFF_CONV_STAB_SCALING");
  conv_stab_scaling_ =
      Core::UTILS::IntegralValue<Inpar::XFEM::ConvStabScaling>(params_xf_stab, "CONV_STAB_SCALING");

  //--------------------------------------------
  // settings for computation of Nitsche's penalty parameter
  //---------------------------------------------

  mass_conservation_combo_ = Core::UTILS::IntegralValue<Inpar::XFEM::MassConservationCombination>(
      params_xf_stab, "MASS_CONSERVATION_COMBO");
  mass_conservation_scaling_ = Core::UTILS::IntegralValue<Inpar::XFEM::MassConservationScaling>(
      params_xf_stab, "MASS_CONSERVATION_SCALING");

  intterms_prev_state_ = Core::UTILS::IntegralValue<Inpar::XFEM::InterfaceTermsPreviousState>(
      params_xf_gen, "INTERFACE_TERMS_PREVIOUS_STATE");

  //--------------------------------------------
  // CONSISTENCY CHECKS for PARAMETERS
  //--------------------------------------------

  // finally, perform the consistency check
  check_parameter_consistency(myrank);

  return;
}

FOUR_C_NAMESPACE_CLOSE
