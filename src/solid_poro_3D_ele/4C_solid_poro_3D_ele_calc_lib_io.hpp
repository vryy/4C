// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_IO_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_IO_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

#include <Teuchos_ParameterList.hpp>


FOUR_C_NAMESPACE_OPEN

namespace Discret::Elements
{

  template <typename T>
  inline std::vector<char>& get_couplstress_data(const T& ele, const Teuchos::ParameterList& params)
  {
    if (ele.is_solid_params_interface())
    {
      return *ele.get_solid_params_interface().coupling_stress_data_ptr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("couplstress");
    }
  }

  template <typename T>
  inline Inpar::Solid::StressType get_io_couplstress_type(
      const T& ele, Teuchos::ParameterList& params)
  {
    if (ele.is_solid_params_interface())
    {
      return ele.get_solid_params_interface().get_coupling_stress_output_type();
    }
    else
    {
      return params.get<Inpar::Solid::StressType>(
          "iocouplstress", Inpar::Solid::StressType::stress_none);
    }
  }

}  // namespace Discret::Elements
FOUR_C_NAMESPACE_CLOSE

#endif