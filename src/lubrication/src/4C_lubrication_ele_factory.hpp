// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LUBRICATION_ELE_FACTORY_HPP
#define FOUR_C_LUBRICATION_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    // forward declaration
    class LubricationEleInterface;

    class LubricationFactory
    {
     public:
      //! ctor
      LubricationFactory() { return; }

      //! dtor
      virtual ~LubricationFactory() = default;
      //! ProvideImpl
      static LubricationEleInterface* provide_impl(
          Core::FE::CellType distype, const std::string& disname);

     private:
      //! define LubricationEle instances dependent on problem
      template <Core::FE::CellType distype, int probdim>
      static LubricationEleInterface* define_problem_type(const std::string& disname);


    };  // end class LubricationFactory

  }  // namespace Elements

}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
