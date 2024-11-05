// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FPSIStructureWrapper : public FSIStructureWrapper
  {
   public:
    /// constructor
    explicit FPSIStructureWrapper(std::shared_ptr<Structure> structure);

    /*!
    \brief extract interface displacements at \f$t_{n}\f$

    \param FPSI (in) : if true perform some pre-stress and fpsi specific stuff

    \note if param FPSI = false the base class version of this method is called

    */
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_dispn(
        bool FPSI = false);

    /*!
    \brief  extract interface displacements at \f$t_{n+1}\f$

    \param FPSI (in) : if true perform some pre-stress and fpsi specific stuff

    \note if param FPSI = false the base class version of this method is called

    */
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_dispnp(
        bool FPSI = false);
  };
}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
