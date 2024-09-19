/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_IO_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_IO_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

#include <Teuchos_ParameterList.hpp>


FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
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

}  // namespace Discret::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif