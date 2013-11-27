/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_poro.cpp

\brief Evaluation of general fluid parameter for fluid in poroelast problem

FluidEleParameter::SetElementPoroParameter(Teuchos::ParameterList& params)
set all general porofluid parameter once for all elements.

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/
#include "fluid_ele_parameter_poro.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_globalproblem.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterPoro* DRT::ELEMENTS::FluidEleParameterPoro::Instance( bool create )
{
  static FluidEleParameterPoro* instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleParameterPoro();
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
void DRT::ELEMENTS::FluidEleParameterPoro::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidEleParameterPoro::FluidEleParameterPoro()
  : DRT::ELEMENTS::FluidEleParameter::FluidEleParameter(),
    set_fluid_parameter_poro_(false),
    poro_conti_partint_(false)
{
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
}

//----------------------------------------------------------------------*
//  set poro parameters                                      vuong 11/12|
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterPoro::SetElementPoroParameter( Teuchos::ParameterList& params, int myrank)
{
  SetElementGeneralFluidParameter(params,myrank);

  set_fluid_parameter_poro_ = true;
  poro_conti_partint_ = params.get<bool>("conti partial integration",false);
  reaction_= true;
  reaction_topopt_= false;
  darcy_= true;
  graddiv_=false;

  if (DRT::Problem::Instance()->ProblemType()==prb_fpsi)
  {
    Teuchos::ParameterList stablist = params.sublist("POROUS-FLOW STABILIZATION");

    // no safety check necessary since all options are used
    tds_      = DRT::INPUT::IntegralValue<INPAR::FLUID::SubscalesTD>(stablist,"TDS");
    transient_= DRT::INPUT::IntegralValue<INPAR::FLUID::Transient>(stablist,"TRANSIENT");
    pspg_     = DRT::INPUT::IntegralValue<int>(stablist,"PSPG");
    supg_     = DRT::INPUT::IntegralValue<int>(stablist,"SUPG");
    vstab_    = DRT::INPUT::IntegralValue<INPAR::FLUID::VStab>(stablist,"VSTAB");
    rstab_    = DRT::INPUT::IntegralValue<INPAR::FLUID::RStab>(stablist,"RSTAB");
    cross_    = DRT::INPUT::IntegralValue<INPAR::FLUID::CrossStress>(stablist,"CROSS-STRESS");
    reynolds_ = DRT::INPUT::IntegralValue<INPAR::FLUID::ReynoldsStress>(stablist,"REYNOLDS-STRESS");

    // overrule higher_order_ele if input-parameter is set
    // this might be interesting for fast (but slightly
    // less accurate) computations
    is_inconsistent_ = DRT::INPUT::IntegralValue<int>(stablist,"INCONSISTENT");
    //-------------------------------
    // get tau definition
    //-------------------------------

    whichtau_ =  DRT::INPUT::IntegralValue<INPAR::FLUID::TauType>(stablist,"DEFINITION_TAU");
    // check if tau can be handled
    if (not(whichtau_ ==
        INPAR::FLUID::tau_franca_madureira_valentin_badia_codina or
        INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt))
      dserror("Definition of Tau cannot be handled by the element");

    // set correct stationary definition of stabilization parameter automatically
    if (fldparatimint_->IsStationary())
    {
      if (whichtau_ == INPAR::FLUID::tau_franca_madureira_valentin_badia_codina)
        whichtau_ = INPAR::FLUID::tau_franca_madureira_valentin_badia_codina_wo_dt;
    }

    //---------------------------------
    // set flags for potential evaluation of tau and material law at int. point
    // default value: evaluation at element center
    const std::string tauloc = stablist.get<std::string>("EVALUATION_TAU");
    if (tauloc == "integration_point") tau_gp_ = true;
    else                               tau_gp_ = false;
    const std::string matloc = stablist.get<std::string>("EVALUATION_MAT");
    if (matloc == "integration_point") mat_gp_ = true;
    else                               mat_gp_ = false;
  }

}

//----------------------------------------------------------------------*/
// print fluid parameter to screen                          rauch 11/13 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidEleParameterPoro::PrintFluidParameterPoro()
{
  std::cout << std::endl << "|-----------------------------------------------------------------------------" << std::endl;
  std::cout << "|  Poro Fluid parameter: " << std::endl;
  std::cout << "|-----------------------------------------------------------------------------" << std::endl;
  //! flag SetGeneralParameter was called
  std::cout << "|    method SetElementParameterPoro was called:    " << set_fluid_parameter_poro_ << std::endl;
  //! flag to (de)activate stationary formulation
  std::cout << "|    Partial integration of conti equation:    " << poro_conti_partint_ << std::endl;
  //! flag to (de)activate Newton linearization
  std::cout << "|    Type of stabilization:    " << stabtype_ << std::endl;

  std::cout << "|---------------------------------------------------------------------------" << std::endl;

}
