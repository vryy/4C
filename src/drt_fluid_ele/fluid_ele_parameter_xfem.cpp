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
    visc_stab_trace_estimate_(INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk),
    visc_stab_hk_(INPAR::XFEM::ViscStab_hk_vol_equivalent),
    visc_stab_gamma_(0.0)
{
}

void DRT::ELEMENTS::FluidEleParameterXFEM::SetElementXFEMParameter(
    Teuchos::ParameterList& params   ///< parameter list
    )
{

  //----------------------------------------------------------------
  Teuchos::ParameterList& params_xfem = params.sublist("XFEM");

  vcellgausspts_  = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  bcellgausspts_  = DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  //----------------------------------------------------------------
  //Teuchos::ParameterList&   params_xf_gen  = params.sublist("XFLUID DYNAMIC/GENERAL");
  // TODO: fill

  //----------------------------------------------------------------
  Teuchos::ParameterList&   params_xf_stab = params.sublist("XFLUID DYNAMIC/STABILIZATION");

  visc_stab_trace_estimate_ = DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE");
  //TODO: add safety checks for a reasonable usage of the eigenvalue problem

  //TODO: exclude cut-based stab_hk for embedded-sided formulations

  visc_stab_hk_             =  DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_hk>(params_xf_stab, "VISC_STAB_HK");


  // TODO: add a comment how to define the visc-stab-fac for eigenvalue problem

  // TODO add or adapt a comment like this when using the eigenvalue problem
      // be at least 2 to get a stable formulation. To reach the value corresponding to alpha
      // equal to 35 it should be higher than 2.

  visc_stab_gamma_          = params_xf_stab.get<double>("VISC_STAB_FAC");

  return;
}

