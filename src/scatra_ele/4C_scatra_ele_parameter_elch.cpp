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
#include "4C_scatra_ele_parameter_elch.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterElch* Discret::ELEMENTS::ScaTraEleParameterElch::instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map =
      Core::UTILS::MakeSingletonMap<std::string>([](const std::string& disname)
          { return std::unique_ptr<ScaTraEleParameterElch>(new ScaTraEleParameterElch(disname)); });

  return singleton_map[disname].instance(Core::UTILS::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterElch::ScaTraEleParameterElch(
    const std::string& disname  //!< name of discretization
    )
    : boundaryfluxcoupling_(true),
      equpot_(Inpar::ElCh::equpot_undefined),
      faraday_(0.0),
      gas_constant_(0.0),
      epsilon_(Inpar::ElCh::epsilon_const),
      frt_(0.0),
      temperature_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterElch::set_parameters(Teuchos::ParameterList& parameters)
{
  // coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann
  // boundaries
  boundaryfluxcoupling_ = parameters.get<bool>("boundaryfluxcoupling");

  // type of closing equation for electric potential
  equpot_ = Core::UTILS::GetAsEnum<Inpar::ElCh::EquPot>(parameters, "equpot");
  if (equpot_ == Inpar::ElCh::equpot_undefined)
    FOUR_C_THROW("Invalid type of closing equation for electric potential!");

  // get parameters
  faraday_ = parameters.get<double>("faraday", -1.0);
  gas_constant_ = parameters.get<double>("gas_constant", -1.0);
  frt_ = parameters.get<double>("frt", -1.0);
  temperature_ = parameters.get<double>("temperature", -1.0);

  // safety checks
  if (frt_ <= 0.0) FOUR_C_THROW("Factor F/RT is non-positive!");
  if (faraday_ <= 0.0) FOUR_C_THROW("Faraday constant is non-positive!");
  if (gas_constant_ <= 0.0) FOUR_C_THROW("(universal) gas constant is non-positive!");
  if (temperature_ < 0.0) FOUR_C_THROW("temperature is non-positive!");
}

FOUR_C_NAMESPACE_CLOSE
