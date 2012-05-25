/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_parameter.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer_ele_parameter.H"


using namespace std;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::TopOptParam> DRT::ELEMENTS::TopOptParam::instance_;


//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::TopOptParam> DRT::ELEMENTS::TopOptParam::Instance()
{
  if (instance_==Teuchos::null)
    instance_ = Teuchos::rcp(new TopOptParam());
  return instance_;
}


//----------------------------------------------------------------------*/
// private constructor of FluidImplParameter
//----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptParam::TopOptParam()
: dissipation_(false),
  pressure_drop_(false),
  dissipation_fac_(0.0),
  pressure_drop_fac_(0.0),
  is_stationary_(false),
  timealgo_(INPAR::FLUID::timeint_one_step_theta),
  dt_(-1.0),
  max_timesteps_(-1)
{
}


//----------------------------------------------------------------------*
//  set general parameters                                   ehrl 04/10 |
//---------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptParam::SetGeneralOptimizationParameter( Teuchos::ParameterList& params )
{
//----------------------------------------------------------------------
// get flags to switch on/off different fluid formulations
//----------------------------------------------------------------------

  // set flag, time integration scheme
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

  // set time integration scheme-specific element parameters
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    is_stationary_ = true;
  else
    is_stationary_ = false;

  // set if objective contains dissipation and according factor
  dissipation_ = params.get<bool>("dissipation");
  if (dissipation_) dissipation_fac_ = params.get<double>("dissipation_fac");
  else              dissipation_fac_ = 0.0;

  // set if objective contains pressure drop and according factor
  pressure_drop_ = params.get<bool>("pres_drop");
  if (pressure_drop_) pressure_drop_fac_ = params.get<double>("pres_drop_fac");
  else                pressure_drop_fac_ = 0.0;

  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  if (not is_stationary_)
  {
    dt_ = params.get<double>("dt");
    max_timesteps_ = params.get<int>("maxtimesteps");
  }
  else
  {
    dt_ = 1.0;
    max_timesteps_ = 1;
  }
}



//----------------------------------------------------------------------*/
// print fluid parameter to screen (AE 01-11)
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::TopOptParam::PrintAdjointParameter() const
{
  cout << endl << "|-----------------------------------------------------------------------------" << endl;
  cout << "|  General optimization parameter: " << endl;
  cout << "|-----------------------------------------------------------------------------" << endl;
  // boolean is true if objective contains dissipation
  cout << "|    objective dissipation on?    " << dissipation_ << endl;
  // boolean is true if objective contains pressure drop
  cout << "|    objective pressure drop on?    " << pressure_drop_ << endl;
  // objective's dissipation factor
  cout << "|    objective dissipation factor:    " << dissipation_fac_ << endl;
  // objective's pressure drop factor
  cout << "|    objective pressure drop factor:    " << pressure_drop_fac_ << endl;

  cout << endl << "|---------------------------------------------------------------------------" << endl;
  cout << "|  Flow parameter: " << endl;
  cout << "|---------------------------------------------------------------------------" << endl;
  //! flag to (de)activate stationary formulation
  cout << "|    steady state:    " << is_stationary_ << endl;
  //! time algorithm
  cout << "|    time algorithm:    " << timealgo_ << endl;
  //! time-step length
  cout << "|    time step:    " << dt_ << endl;
  /// maximal number of time steps
  cout << "|    maximal number of time steps:     " << max_timesteps_ << endl;
  cout << "|---------------------------------------------------------------------------" << endl;
}  /// @name objective parameters

