/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static time integration parameters required for scalar transport
element evaluation

This singleton class holds all static time integration parameters required for scalar transport
element evaluation. All parameters are usually set only once at the beginning of a simulation,
namely during initialization of the global time integrator, and then never touched again throughout
the simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

\maintainer Anh-Tu Vuong

\level 1
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"
#include "scatra_ele_parameter_timint.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterTimInt* DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(
    const std::string& disname,                //!< name of discretization
    const ScaTraEleParameterTimInt* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterTimInt*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterTimInt();
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterTimInt*>::iterator i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  // return existing or newly created instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 08/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterTimInt::Done()
{
  // delete singleton
  Instance("", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 08/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterTimInt::ScaTraEleParameterTimInt()
    : is_genalpha_(false),
      is_stationary_(false),
      is_incremental_(false),
      time_(-1.0),
      dt_(0.0),
      timefac_(0.0),
      timefacrhs_(0.0),
      timefacrhstau_(0.0),
      alphaF_(0.0)
{
  return;
}


//----------------------------------------------------------------------*/
// set time parameters which are equal for every fluid  rasthofer 11/13 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterTimInt::SetParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get control parameters
  is_stationary_ = parameters.get<bool>("using stationary formulation");
  is_genalpha_ = parameters.get<bool>("using generalized-alpha time integration");
  is_incremental_ = parameters.get<bool>("incremental solver");

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
      if (not is_incremental_) timefacrhs_ *= (1.0 - alphaF_);
    }
    else
    {
      if (not is_incremental_)
        timefacrhs_ = 0.0;
      else
        timefacrhs_ = timefac_;
    }
  }
  else
  {
    if (not is_incremental_) timefacrhs_ = 0.0;
  }

  return;
}


//----------------------------------------------------------------------*/
// print fluid time parameter to screen                 rasthofer 11/13 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterTimInt::PrintScaTraTimeParameter()
{
  //  std::cout << std::endl << "|-----------------------------------------------------------" <<
  //  std::endl; std::cout << "|  Fluid Element Time parameter: " << std::endl; std::cout <<
  //  "|-------------------------------------------------------------------" << std::endl;
  //  //! time parameters set?
  //  std::cout << "|    time algorithm:         " << timealgo_ << std::endl;
  //  //! is stationary?
  //  std::cout << "|    is stationary?:         " << is_stationary_ << std::endl;
  //  //! time algorithm
  //  std::cout << "|    time parameters set:    " << set_general_fluid_timeparameter_ << std::endl;
  //  //! actual time to evaluate the body BC
  //  std::cout << "|    time:                   " << time_ << std::endl;
  //  //! time-step length
  //  std::cout << "|    time step:              " << dt_ << std::endl;
  //  //! timefac = dt_ * ("pseudo"-)theta_
  //  std::cout << "|    time factor:            " << timefac_ << std::endl;
  //  //! factor for left-hand side due to one-step-theta time-integration scheme
  //  std::cout << "|    theta:                  " << theta_ << std::endl;
  //  //! factor for right-hand side due to one-step-theta time-integration scheme
  //  std::cout << "|    (1-theta):              " << omtheta_ << std::endl;
  //  //! generalised-alpha parameter (connecting velocity and acceleration)
  //  std::cout << "|    gamma:                  " << gamma_ << std::endl;
  //  //! generalised-alpha parameter (velocity)
  //  std::cout << "|    alpha_F:                " << alphaF_ << std::endl;
  //  //! generalised-alpha parameter (acceleration)
  //  std::cout << "|    alpha_M:                " << alphaM_ << std::endl;
  //  //! generalised-alpha parameter, alphaF_*gamma_*dt_
  //  std::cout << "|    time factor mat_u:      " << afgdt_ << std::endl;
  //  //! time integration factor for the right hand side (boundary elements)
  //  std::cout << "|    time factor rhs:        " << timefacrhs_ << std::endl;
  //  //! time integration factor for the left hand side (pressure)
  //  std::cout << "|    time factor mat_p:      " << timefacpre_ << std::endl;
  //  std::cout << std::endl << "|-----------------------------------------------------------" <<
  //  std::endl;

  return;
}
