// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_FACTORY_HPP
#define FOUR_C_SCATRA_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_scatra_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    class ScaTraFactory
    {
     public:
      //! ctor
      ScaTraFactory() = default;

      //! dtor
      virtual ~ScaTraFactory() = default;

      //! ProvideImpl
      static ScaTraEleInterface* provide_impl(Core::FE::CellType distype,
          Inpar::ScaTra::ImplType problem, const int numdofpernode, const int numscal,
          const std::string& disname);

      //! ProvideImplHDG
      static ScaTraEleInterface* provide_impl_hdg(Core::FE::CellType distype,
          Inpar::ScaTra::ImplType problem, const int numdofpernode, const int numscal,
          const std::string& disname);

     private:
      //! define ScatraEle instances dependent on problem
      template <Core::FE::CellType distype, int probdim>
      static ScaTraEleInterface* define_problem_type(Inpar::ScaTra::ImplType problem,
          const int numdofpernode, const int numscal, const std::string& disname);

      //! define ScatraEle instances dependent on problem
      template <Core::FE::CellType distype, int probdim>
      static ScaTraEleInterface* define_problem_type_hdg(Inpar::ScaTra::ImplType problem,
          const int numdofpernode, const int numscal, const std::string& disname);

    };  // end class ScaTraFactory

  }  // namespace Elements

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
