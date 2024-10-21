// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_parameter_lsreinit.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterLsReinit*
Discret::ELEMENTS::ScaTraEleParameterLsReinit::instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [](const std::string& disname) {
        return std::unique_ptr<ScaTraEleParameterLsReinit>(new ScaTraEleParameterLsReinit(disname));
      });

  return singleton_map[disname].instance(Core::Utils::SingletonAction::create, disname);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraEleParameterLsReinit::ScaTraEleParameterLsReinit(
    const std::string& disname  //!< name of discretization
    )
    : reinittype_(Inpar::ScaTra::reinitaction_none),
      signtype_(Inpar::ScaTra::signtype_nonsmoothed),
      charelelengthreinit_(Inpar::ScaTra::root_of_volume_reinit),
      interfacethicknessfac_(1.0),
      useprojectedreinitvel_(false),
      linform_(Inpar::ScaTra::fixed_point),
      artdiff_(Inpar::ScaTra::artdiff_none),
      alphapen_(0.0),
      project_(true),
      projectdiff_(0.0),
      lumping_(false),
      difffct_(Inpar::ScaTra::hyperbolic)
{
}


//----------------------------------------------------------------------*
//  set parameters                                      rasthofer 12/13 |
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::ScaTraEleParameterLsReinit::set_parameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get reinitialization parameters list
  Teuchos::ParameterList& reinitlist = parameters.sublist("REINITIALIZATION");

  // reinitialization strategy
  reinittype_ =
      Teuchos::getIntegralValue<Inpar::ScaTra::ReInitialAction>(reinitlist, "REINITIALIZATION");

  // get signum function
  signtype_ =
      Teuchos::getIntegralValue<Inpar::ScaTra::SmoothedSignType>(reinitlist, "SMOOTHED_SIGN_TYPE");

  // characteristic element length for signum function
  charelelengthreinit_ = Teuchos::getIntegralValue<Inpar::ScaTra::CharEleLengthReinit>(
      reinitlist, "CHARELELENGTHREINIT");

  // interface thickness for signum function
  interfacethicknessfac_ = reinitlist.get<double>("INTERFACE_THICKNESS");

  // form of linearization for nonlinear terms
  linform_ = Teuchos::getIntegralValue<Inpar::ScaTra::LinReinit>(reinitlist, "LINEARIZATIONREINIT");

  // set form of velocity evaluation
  Inpar::ScaTra::VelReinit velreinit =
      Teuchos::getIntegralValue<Inpar::ScaTra::VelReinit>(reinitlist, "VELREINIT");
  if (velreinit == Inpar::ScaTra::vel_reinit_node_based) useprojectedreinitvel_ = true;

  // set flag for artificial diffusion term
  artdiff_ = Teuchos::getIntegralValue<Inpar::ScaTra::ArtDiff>(reinitlist, "ARTDIFFREINIT");

  // set penalty parameter for elliptic reinitialization
  alphapen_ = reinitlist.get<double>("PENALTY_PARA");

  // get diffusivity function
  difffct_ = Teuchos::getIntegralValue<Inpar::ScaTra::DiffFunc>(reinitlist, "DIFF_FUNC");

  // L2-projection
  project_ = reinitlist.get<bool>("PROJECTION");

  // diffusion for L2-projection
  projectdiff_ = reinitlist.get<double>("PROJECTION_DIFF");
  if (projectdiff_ < 0.0) FOUR_C_THROW("Diffusivity has to be positive!");

  // lumping for L2-projection
  lumping_ = reinitlist.get<bool>("LUMPING");

  // check for illegal combination
  if (projectdiff_ > 0.0 and lumping_ == true) FOUR_C_THROW("Illegal combination!");
  if (projectdiff_ > 0.0 and reinittype_ == Inpar::ScaTra::reinitaction_sussman)
    FOUR_C_THROW("Illegal combination!");
  // The second FOUR_C_THROW is added here for safety reasons. I think that using a diffusive term
  // for the reconstruction of the velocity for reinitialization is possible, but I have not yet
  // further investigated this option. Therefore, you should test it first.

  return;
}

FOUR_C_NAMESPACE_CLOSE
