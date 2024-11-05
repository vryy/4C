// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_SCATRA_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_POROELAST_SCATRA_UTILS_CLONESTRATEGY_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace PoroElastScaTra
{
  namespace Utils
  {

    class PoroelastCloneStrategyforScatraElements : public PoroElast::Utils::PoroelastCloneStrategy
    {
     public:
      //! constructor
      explicit PoroelastCloneStrategyforScatraElements() = default;
      //! determine element type string and whether element is copied or not
      bool determine_ele_type(
          Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);
    };

    //! \brief implementation of special clone strategy for automatic generation
    //!        of scatra discretization from a given structure discretization for porous media
    class PoroScatraCloneStrategy
    {
     public:
      //! constructor
      explicit PoroScatraCloneStrategy() = default;
      //! destructor
      virtual ~PoroScatraCloneStrategy() = default;

      //! return ScaTra::ImplType of the element
      Inpar::ScaTra::ImplType get_impl_type(
          Core::Elements::Element* ele  //! element whose ScaTra::ImplType shall be determined
      );

     protected:
      //! determine element type string and whether element is copied or not
      bool determine_ele_type(
          Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      //! check for correct material
      void check_material_type(const int matid);

      //! provide cloned element with element specific data (material etc.)
      void set_element_data(std::shared_ptr<Core::Elements::Element>
                                newele,     //! current cloned element on target discretization
          Core::Elements::Element* oldele,  //! current element on source discretization
          const int matid,                  //! material of cloned element
          const bool isnurbs                //! nurbs flag
      );

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> conditions_to_copy() const;
    };
  }  // namespace Utils
}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
