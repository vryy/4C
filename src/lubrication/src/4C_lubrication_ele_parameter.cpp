// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_lubrication_ele_parameter.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | singleton access method                                  wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::LubricationEleParameter* Discret::ELEMENTS::LubricationEleParameter::instance(
    const std::string& disname  //!< name of discretization
)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [](const std::string& disname)
      { return std::unique_ptr<LubricationEleParameter>(new LubricationEleParameter(disname)); });

  return singleton_map[disname].instance(Core::Utils::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                       wirtz 10/15 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::LubricationEleParameter::LubricationEleParameter(
    const std::string& disname  //!< name of discretization
    )

    : time_(-1.0),
      modified_reynolds_(true),
      addsqz_(true),
      purelub_(true),
      roughness_deviation_(0.0)
{
  return;
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     wirtz 10/15 |
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationEleParameter::set_time_parameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  // get current time and time-step length
  time_ = parameters.get<double>("total time");
}

//----------------------------------------------------------------------*/
// set parameters which are equal for every lubrication     wirtz 10/15 |
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::LubricationEleParameter::set_general_parameters(
    Teuchos::ParameterList& parameters  //!< parameter list
)
{
  modified_reynolds_ = parameters.get<bool>("ismodifiedrey");
  addsqz_ = parameters.get<bool>("addsqz");
  purelub_ = parameters.get<bool>("purelub");
  roughness_deviation_ = parameters.get<double>("roughnessdeviation");
}

FOUR_C_NAMESPACE_CLOSE
