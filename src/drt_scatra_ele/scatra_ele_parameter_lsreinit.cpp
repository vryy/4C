/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_lsreinit.cpp

\brief singleton class holding all static levelset reinitialization parameters required for element
evaluation

This singleton class holds all static levelset reinitialization parameters required for element
evaluation. All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237

\level 2
</pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_parameter_lsreinit.H"
#include "../drt_lib/drt_dserror.H"

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterLsReinit* DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance(
    const std::string& disname,                  //!< name of discretization
    const ScaTraEleParameterLsReinit* delete_me  //!< creation/destruction indication
)
{
  // each discretization is associated with exactly one instance of this class according to a static
  // map
  static std::map<std::string, ScaTraEleParameterLsReinit*> instances;

  // check whether instance already exists for current discretization, and perform instantiation if
  // not
  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleParameterLsReinit(disname);
  }

  // destruct instance
  else
  {
    for (std::map<std::string, ScaTraEleParameterLsReinit*>::iterator i = instances.begin();
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

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterLsReinit::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance("", this);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterLsReinit::ScaTraEleParameterLsReinit(
    const std::string& disname  //!< name of discretization
    )
    : reinittype_(INPAR::SCATRA::reinitaction_none),
      signtype_(INPAR::SCATRA::signtype_nonsmoothed),
      charelelengthreinit_(INPAR::SCATRA::root_of_volume_reinit),
      interfacethicknessfac_(1.0),
      useprojectedreinitvel_(false),
      linform_(INPAR::SCATRA::fixed_point),
      artdiff_(INPAR::SCATRA::artdiff_none),
      alphapen_(0.0),
      project_(true),
      projectdiff_(0.0),
      lumping_(false),
      difffct_(INPAR::SCATRA::hyperbolic)
{
}


//----------------------------------------------------------------------*
//  set parameters                                      rasthofer 12/13 |
//----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterLsReinit::SetParameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get reinitialization parameters list
  Teuchos::ParameterList& reinitlist = parameters.sublist("REINITIALIZATION");

  // reinitialization strategy
  reinittype_ =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialAction>(reinitlist, "REINITIALIZATION");

  // get signum function
  signtype_ =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::SmoothedSignType>(reinitlist, "SMOOTHED_SIGN_TYPE");

  // characteristic element length for signum function
  charelelengthreinit_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLengthReinit>(
      reinitlist, "CHARELELENGTHREINIT");

  // interface thickness for signum function
  interfacethicknessfac_ = reinitlist.get<double>("INTERFACE_THICKNESS");

  // form of linearization for nonlinear terms
  linform_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::LinReinit>(reinitlist, "LINEARIZATIONREINIT");

  // set form of velocity evaluation
  INPAR::SCATRA::VelReinit velreinit =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::VelReinit>(reinitlist, "VELREINIT");
  if (velreinit == INPAR::SCATRA::vel_reinit_node_based) useprojectedreinitvel_ = true;

  // set flag for artificial diffusion term
  artdiff_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::ArtDiff>(reinitlist, "ARTDIFFREINIT");

  // set penalty parameter for elliptic reinitialization
  alphapen_ = reinitlist.get<double>("PENALTY_PARA");

  // get diffusivity function
  difffct_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::DiffFunc>(reinitlist, "DIFF_FUNC");

  // L2-projection
  project_ = DRT::INPUT::IntegralValue<bool>(reinitlist, "PROJECTION");

  // diffusion for L2-projection
  projectdiff_ = reinitlist.get<double>("PROJECTION_DIFF");
  if (projectdiff_ < 0.0) dserror("Diffusivity has to be positive!");

  // lumping for L2-projection
  lumping_ = DRT::INPUT::IntegralValue<bool>(reinitlist, "LUMPING");

  // check for illegal combination
  if (projectdiff_ > 0.0 and lumping_ == true) dserror("Illegal combination!");
  if (projectdiff_ > 0.0 and reinittype_ == INPAR::SCATRA::reinitaction_sussman)
    dserror("Illegal combination!");
  // The second dserror is added here for safety reasons. I think that using a diffusive term for
  // the reconstruction of the velocity for reinitialization is possible, but I have not yet further
  // investigated this option. Therefore, you should test it first.

  return;
}
