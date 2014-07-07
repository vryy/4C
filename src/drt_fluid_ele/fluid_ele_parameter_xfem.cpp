/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_xfem.cpp

\brief Setting of specific XFEM based fluid parameter for element evaluation

<pre>
Maintainers: Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#include "fluid_ele_parameter_xfem.H"
#include "../drt_io/io_pstream.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterXFEM* DRT::ELEMENTS::FluidEleParameterXFEM::Instance( bool create )
{
  static FluidEleParameterXFEM* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleParameterXFEM();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterXFEM::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterXFEM::FluidEleParameterXFEM()
  : DRT::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd(),
    vcellgausspts_(INPAR::CUT::VCellGaussPts_DirectDivergence),
    bcellgausspts_(INPAR::CUT::BCellGaussPts_Tessellation),
    coupling_method_(INPAR::XFEM::Nitsche),
    coupling_strategy_(INPAR::XFEM::Xfluid_Sided_Coupling),
    hybrid_lm_l2_proj_(INPAR::XFEM::Hybrid_LM_L2_Proj_part),
    visc_stab_trace_estimate_(INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk),
    visc_stab_hk_(INPAR::XFEM::ViscStab_hk_vol_equivalent),
    visc_stab_gamma_(0.0),
    is_visc_adjoint_symmetric_(true),
    is_pseudo_2D_(false),
    xff_conv_stab_scaling_(INPAR::XFEM::XFF_ConvStabScaling_none),
    conv_stab_scaling_(INPAR::XFEM::ConvStabScaling_none),
    is_velgrad_interface_stab_(false),
    velgrad_interface_stab_fac_(0.0),
    is_pressure_gradient_interface_stab_(false),
    pressure_interface_stab_fac_(0.0)
{
}

//----------------------------------------------------------------------*/
//    check parameter combination for consistency
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterXFEM::CheckParameterConsistency(int myrank) const
{
  // short consistency check
  const bool isEmbNitsche = (coupling_method_ == INPAR::XFEM::Nitsche && coupling_strategy_ == INPAR::XFEM::Embedded_Sided_Coupling);

  if (visc_stab_trace_estimate_ == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue && !isEmbNitsche)
    dserror("Solution of eigenvalue problem to estimate parameter from trace inequality is only reasonable for embedded-sided Nitsche coupling.");

  if (is_visc_adjoint_symmetric_ && coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress && coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_Coupling)
  {
    if (myrank == 0)
      IO::cout << "Be warned: the symmetric hybrid/viscous stress-based LM approach is known for unstable behaviour in xfluid-fluid problems." << IO::endl;
  }


  // Consistency Checks for characteristic element length definitions
  if(isEmbNitsche)
  {
    // check element l
    if(ViscStabTracEstimate() != INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
    {
      if(ViscStabHK() == INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf or
         ViscStabHK() == INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf)
        dserror("chosen characteristic element length definition ViscStabHK is not supported for embedded-sided Nitsche method");

      if(ViscStabHK() == INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf)
        dserror("chosen characteristic element length definition ViscStab_hk_ele_vol_div_by_max_ele_surf is supported for embedded-sided Nitsche method,"
            "however not a reasonable choice, as it can lead to an ill-conditioning due to a to large penalty scaling for anisotropic embedded meshes!"
            "If you want, you can try it!");
    }
  }

  // Consistency Check for EVP and Xfluid-sided Nitsche
  if(!isEmbNitsche and ViscStabTracEstimate() == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
    dserror("estimating trace inequality scaling via solving eigenvalue problems not supported for xfluid-sided Nitsche!");

  if(!isEmbNitsche and ViscStabHK() == INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf)
    dserror("chosen characteristic element length definition ViscStabHK is not supported for xfluid-sided Nitsche method as the element surface cannot be specified for cut elements");

  if(!isEmbNitsche and
      (ViscStabHK() == INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf or
       ViscStabHK() == INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf ))
    IO::cout << "Be warned: the chosen characteristic element length definition ViscStabHK can become critical for xfluid-sided Nitsche method as the current definition either can not guarantee a sufficient estimate of the inverse inequality or can lead to cut position dependent error behaviour!" << IO::endl;

}


//----------------------------------------------------------------------*/
//    set parameters
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterXFEM::SetElementXFEMParameter(
    Teuchos::ParameterList& params,  ///< parameter list
    int myrank                       ///< pid (for output purpose)
    )
{

  //----------------------------------------------------------------
  Teuchos::ParameterList& params_xfem = params.sublist("XFEM");

  vcellgausspts_  = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  bcellgausspts_  = DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  //----------------------------------------------------------------
  //Teuchos::ParameterList&   params_xf_gen  = params.sublist("XFLUID DYNAMIC/GENERAL");

  //----------------------------------------------------------------
  Teuchos::ParameterList&   params_xf_stab = params.sublist("XFLUID DYNAMIC/STABILIZATION");

  //--------------------------------------------
  // parameters describing the coupling approach
  //---------------------------------------------
  coupling_method_ = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingMethod>(params_xf_stab,"COUPLING_METHOD");

  coupling_strategy_  = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(params_xf_stab,"COUPLING_STRATEGY");

  hybrid_lm_l2_proj_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab,"HYBRID_LM_L2_PROJ");

  //--------------------------------------------
  // parameters for the viscous stabilization in Nitsche's method and MixedHybrid_LM methods
  //---------------------------------------------

  visc_stab_trace_estimate_   = DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE");

  visc_stab_hk_               = DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_hk>(params_xf_stab, "VISC_STAB_HK");

  visc_stab_gamma_            = params_xf_stab.get<double>("VISC_STAB_FAC");

  is_visc_adjoint_symmetric_  = DRT::INPUT::IntegralValue<bool>(params_xf_stab,"VISC_ADJOINT_SYMMETRY");

  is_pseudo_2D_               = (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab, "IS_PSEUDO_2D");


  // TODO: add a comment how to define the visc-stab-fac for eigenvalue problem
  // TODO add or adapt a comment like this when using the eigenvalue problem
      // be at least 2 to get a stable formulation. To reach the value corresponding to alpha
      // equal to 35 it should be higher than 2.

  //--------------------------------------------
  // parameters for the convective interface stabilizations
  //---------------------------------------------

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");
  conv_stab_scaling_     = DRT::INPUT::IntegralValue<INPAR::XFEM::ConvStabScaling>(params_xf_stab,"CONV_STAB_SCALING");


  //--------------------------------------------
  // special interface stabilizations for fluidfluid problems
  //---------------------------------------------

  is_velgrad_interface_stab_   =  DRT::INPUT::IntegralValue<bool>(params_xf_stab,"VELGRAD_INTERFACE_STAB");

  velgrad_interface_stab_fac_  = params_xf_stab.get<double>("GHOST_PENALTY_FAC",0.0);

  is_pressure_gradient_interface_stab_ = DRT::INPUT::IntegralValue<bool>(params_xf_stab,"PRESSCOUPLING_INTERFACE_STAB");

  pressure_interface_stab_fac_ = params_xf_stab.get<double>("PRESSCOUPLING_INTERFACE_FAC",0.0);


  //--------------------------------------------
  // CONSISTENCY CHECKS for PARAMETERS
  //--------------------------------------------

  // finally, perform the consistency check
  CheckParameterConsistency(myrank);

  return;
}

