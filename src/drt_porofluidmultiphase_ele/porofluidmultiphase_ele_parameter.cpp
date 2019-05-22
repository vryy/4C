/*----------------------------------------------------------------------*/
/*!
 \brief container class holding parameters for element evaluation (singleton)

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/


#include "porofluidmultiphase_ele_parameter.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | singleton access method                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter*
DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(
    const std::string& disname,                       //!< name of discretization
    const PoroFluidMultiPhaseEleParameter* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, PoroFluidMultiPhaseEleParameter*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new PoroFluidMultiPhaseEleParameter(disname);
  }

  // destruct instance given to the destructor
  else
  {
    for (std::map<std::string, PoroFluidMultiPhaseEleParameter*>::iterator i = instances.begin();
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
 | singleton destruction                                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Done()
{
  // delete singleton
  Instance("", this);

  return;
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
      isset_generalparams_(false)
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

  // done
  isset_generalparams_ = true;
}
