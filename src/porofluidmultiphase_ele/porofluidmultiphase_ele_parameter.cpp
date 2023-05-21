/*----------------------------------------------------------------------*/
/*! \file
 \brief container class holding parameters for element evaluation (singleton)

   \level 3

 *----------------------------------------------------------------------*/


#include "porofluidmultiphase_ele_parameter.H"

#include "utils_exceptions.H"
#include "utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 | singleton access method                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter*
DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>(
      [](const std::string& disname)
      {
        return std::unique_ptr<PoroFluidMultiPhaseEleParameter>(
            new PoroFluidMultiPhaseEleParameter(disname));
      });

  return singleton_map[disname].Instance(::UTILS::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                       vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::PoroFluidMultiPhaseEleParameter(
    const std::string& disname  //!< name of discretization
    )
    : time_(-1.0),
      dt_(0.0),
      timefac_(0.0),
      timefacrhs_(0.0),
      timefacrhstau_(0.0),
      alphaF_(0.0),
      is_genalpha_(false),
      is_stationary_(false),
      is_ale_(false),
      stab_biot_(false),
      nds_disp_(-1),
      nds_vel_(-1),
      nds_solidpressure_(-1),
      nds_scalar_(-1),
      isset_generalparams_(false),
      domainint_funct_(0)
{
  return;
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     vuong 08/16 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::SetTimeStepParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  if (not isset_generalparams_)
    dserror("General parameters have to be set before time step parameters!");

  // get current time and time-step length
  time_ = parameters.get<double>("total time");
  dt_ = parameters.get<double>("time-step length");

  // get time factor and alpha_F if required
  // one-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = alphaF * (gamma/alpha_M) * dt

  //-----------------------------------------------------
  //       |          timefac         |    timefacrhs   |
  // ----------------------------------------------------
  // OST   |                  dt*theta                  |
  //-----------------------------------------------------
  // BDF2  |               2/3 * dt                     |
  //-----------------------------------------------------
  // Af GA | alphaF*gamma*dt/alphaM   | gamma*dt/alphaM |
  //-----------------------------------------------------

  timefac_ = 1.0;
  alphaF_ = 1.0;
  timefacrhs_ = 1.0;
  timefacrhstau_ = 1.0;

  if (not is_stationary_)
  {
    timefac_ = parameters.get<double>("time factor");

    if (is_genalpha_)
    {
      alphaF_ = parameters.get<double>("alpha_F");
      timefac_ *= alphaF_;
    }
    if (timefac_ < 0.0) dserror("time factor is negative.");
  }

  if (not is_stationary_)
  {
    if (is_genalpha_)
    {
      timefacrhs_ = timefac_ / alphaF_;
      timefacrhstau_ = timefacrhs_;
    }
    else
    {
      timefacrhs_ = timefac_;
    }
  }
  else
  {
    timefacrhs_ = 0.0;
  }

  return;
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     vuong 08/16 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::SetGeneralParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get control parameters
  is_stationary_ = parameters.get<bool>("using stationary formulation");
  is_genalpha_ = parameters.get<bool>("using generalized-alpha time integration");

  // set ale case
  is_ale_ = parameters.get<bool>("isale", false);
  // set biot stabilization
  stab_biot_ = parameters.get<bool>("stab_biot", false);

  // set number of dof set related to mesh displacements
  nds_disp_ = parameters.get<int>("nds_disp", false);
  // set number of dof set related to mesh velocities
  nds_vel_ = parameters.get<int>("nds_vel", false);
  // set number of dof set related to solid pressure
  nds_solidpressure_ = parameters.get<int>("nds_solidpressure", false);
  // set number of dof set related to solid pressure
  nds_scalar_ = parameters.get<int>("nds_scalar", false);
  // get number of domain integral functions and resize vector
  const int num_domainint_funct = parameters.get<int>("num_domainint_funct", false);
  domainint_funct_.resize(num_domainint_funct);

  // set functions into vector
  for (int ifunct = 0; ifunct < num_domainint_funct; ifunct++)
  {
    domainint_funct_[ifunct] =
        parameters.get<int>("domainint_funct_" + std::to_string(ifunct), false);
  }

  // done
  isset_generalparams_ = true;
}
