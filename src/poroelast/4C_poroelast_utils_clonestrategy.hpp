// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROELAST_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_POROELAST_UTILS_CLONESTRATEGY_HPP



#include "4C_config.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace PoroElast
{
  namespace Utils
  {
    //! \brief implementation of special clone strategy for automatic generation
    //!        of fluid discretization from a given structure discretization
    class PoroelastCloneStrategy
    {
     public:
      //! constructor
      explicit PoroelastCloneStrategy() = default;
      //! destructor
      virtual ~PoroelastCloneStrategy() = default;

      //! returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> conditions_to_copy() const;


      //! determine element type string and whether element is copied or not
      bool determine_ele_type(
          Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

     protected:
      //! set element-specific data (material etc.)
      virtual void set_element_data(std::shared_ptr<Core::Elements::Element> newele,
          Core::Elements::Element* oldele, const int matid, const bool isnurbs);

      //! check for correct material
      void check_material_type(const int matid);

      //! set anisotropic permeability directions onto fluid element
      void set_anisotropic_permeability_directions_onto_fluid(
          std::shared_ptr<Core::Elements::Element> newele, Core::Elements::Element* oldele);

      //! set anisotropic permeability nodal coefficients onto fluid element
      void set_anisotropic_permeability_nodal_coeffs_onto_fluid(
          std::shared_ptr<Core::Elements::Element> newele, Core::Elements::Element* oldele);
    };

  }  // namespace Utils
}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
