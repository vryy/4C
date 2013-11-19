/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_loma.H

\brief Evaluation of general fluid parameter for loma

       As it is currently not feasible to split std and loma terms
       this class holds the loma specific paramters. Though this
       class provides the construction methods of the stdfluid it
       is nice to keep the loma parameters separated, since they are
       definitely not needed for poro / fpsi problems.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/
#include "fluid_ele_parameter_std.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd* DRT::ELEMENTS::FluidEleParameterStd::Instance( bool create )
{
  static FluidEleParameterStd* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleParameterStd();
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
void DRT::ELEMENTS::FluidEleParameterStd::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterStd::FluidEleParameterStd()
  : DRT::ELEMENTS::FluidEleParameter::FluidEleParameter(),
    update_mat_(false),
    conti_supg_(true),
    conti_cross_(INPAR::FLUID::cross_stress_stab_none),
    conti_reynolds_(INPAR::FLUID::reynolds_stress_stab_none),
    multifrac_loma_conti_(false)
{
}

//----------------------------------------------------------------------*
//  set loma parameters                                  rasthofer 03/12|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterStd::SetElementLomaParameter( Teuchos::ParameterList& params )
{
  // get parameter lists
  Teuchos::ParameterList& lomaparams = params.sublist("LOMA");
  Teuchos::ParameterList& stabparams = params.sublist("RESIDUAL-BASED STABILIZATION");
  Teuchos::ParameterList& turbmodelparamsmfs = params.sublist("MULTIFRACTAL SUBGRID SCALES");

  //---------------------------------------------------------------------------------
  // material update with subgrid-scale temperature
  //---------------------------------------------------------------------------------

  update_mat_ = lomaparams.get<bool>("update material",false);

  //---------------------------------------------------------------------------------
  // parameter for additional rbvmm terms in continuity equation
  //---------------------------------------------------------------------------------

  conti_supg_     = DRT::INPUT::IntegralValue<int>(stabparams,"LOMA_CONTI_SUPG");
  conti_cross_    = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stabparams,"LOMA_CONTI_CROSS_STRESS");
  conti_reynolds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stabparams,"LOMA_CONTI_REYNOLDS_STRESS");

  //---------------------------------------------------------------------------------
  // parameter for additional multifractal subgrid-scale terms
  //---------------------------------------------------------------------------------

  if (turb_mod_action_ == INPAR::FLUID::multifractal_subgrid_scales)
   multifrac_loma_conti_ = DRT::INPUT::IntegralValue<int>(turbmodelparamsmfs,"LOMA_CONTI");

  return;
}
