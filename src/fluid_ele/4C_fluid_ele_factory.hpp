// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_FACTORY_HPP
#define FOUR_C_FLUID_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class FluidEleInterface;

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief to do
     *
     *
     * \date March, 2012
     */
    /*--------------------------------------------------------------------------*/
    class FluidFactory
    {
     public:
      //! ctor
      FluidFactory() { return; }

      //! dtor
      virtual ~FluidFactory() = default;
      //! ProvideImpl
      static FluidEleInterface* provide_impl(Core::FE::CellType distype, std::string problem);

      //! special ProvideImpl for XFEM problems to reduce created template combinations
      static FluidEleInterface* provide_impl_xfem(Core::FE::CellType distype, std::string problem);

     private:
      //! define FluidEle instances dependent on problem
      template <Core::FE::CellType distype>
      static FluidEleInterface* define_problem_type(std::string problem);

      //! special define_problem_type_xfem for XFEM problems
      template <Core::FE::CellType distype>
      static FluidEleInterface* define_problem_type_xfem(std::string problem);

    };  // end class FluidFactory

  }  // namespace ELEMENTS

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
