/*----------------------------------------------------------------------*/
/*! \file

\brief singleton class holding all static electrochemistry parameters required for element
evaluation

This singleton class holds all static electrochemistry parameters required for element evaluation.
All parameters are usually set only once at the beginning of a simulation, namely during
initialization of the global time integrator, and then never touched again throughout the
simulation. This parameter class needs to coexist with the general parameter class holding all
general static parameters required for scalar transport element evaluation.

\level 2

*/
/*----------------------------------------------------------------------*/
#include "lib_dserror.H"

#include "scatra_ele_parameter_elch.H"
#include "utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch* DRT::ELEMENTS::ScaTraEleParameterElch::Instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map = ::UTILS::MakeSingletonMap<std::string>([](const std::string& disname)
      { return std::unique_ptr<ScaTraEleParameterElch>(new ScaTraEleParameterElch(disname)); });

  return singleton_map[disname].Instance(::UTILS::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch(
    const std::string& disname  //!< name of discretization
    )
    : boundaryfluxcoupling_(true),
      equpot_(INPAR::ELCH::equpot_undefined),
      faraday_(0.0),
      gas_constant_(0.0),
      epsilon_(INPAR::ELCH::epsilon_const),
      frt_(0.0),
      temperature_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ScaTraEleParameterElch::SetParameters(Teuchos::ParameterList& parameters)
{
  // coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann
  // boundaries
  boundaryfluxcoupling_ = parameters.get<bool>("boundaryfluxcoupling");

  // type of closing equation for electric potential
  equpot_ = DRT::INPUT::get<INPAR::ELCH::EquPot>(parameters, "equpot");
  if (equpot_ == INPAR::ELCH::equpot_undefined)
    dserror("Invalid type of closing equation for electric potential!");

  // get parameters
  faraday_ = parameters.get<double>("faraday", -1.0);
  gas_constant_ = parameters.get<double>("gas_constant", -1.0);
  frt_ = parameters.get<double>("frt", -1.0);
  temperature_ = parameters.get<double>("temperature", -1.0);

  // safety checks
  if (frt_ <= 0.0) dserror("Factor F/RT is non-positive!");
  if (faraday_ <= 0.0) dserror("Faraday constant is non-positive!");
  if (gas_constant_ <= 0.0) dserror("(universal) gas constant is non-positive!");
  if (temperature_ < 0.0) dserror("temperature is non-positive!");
}
