// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ELEMAG_ELE_FACTORY_HPP
#define FOUR_C_ELEMAG_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class ElemagEleInterface;

    class ElemagFactory
    {
     public:
      //! ctor
      ElemagFactory() { return; }

      //! dtor
      virtual ~ElemagFactory() = default;
      //! ProvideImpl
      static ElemagEleInterface* provide_impl(Core::FE::CellType distype, std::string problem);

     private:
      //! define ElemagEle instances dependent on problem
      template <Core::FE::CellType distype>
      static ElemagEleInterface* define_problem_type(std::string problem);
    };

  }  // namespace ELEMENTS

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
