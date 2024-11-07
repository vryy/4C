// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROMULTIPHASE_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_POROMULTIPHASE_UTILS_CLONESTRATEGY_HPP

#include "4C_config.hpp"

#include <map>
#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace POROMULTIPHASE
{
  namespace Utils
  {
    /*!
    \brief implementation of special clone strategy for automatic generation
           of scatra from a given fluid discretization

     */
    class PoroFluidMultiPhaseCloneStrategy
    {
     public:
      /// constructor
      explicit PoroFluidMultiPhaseCloneStrategy() {}
      /// destructor
      virtual ~PoroFluidMultiPhaseCloneStrategy() = default;
      /// returns conditions names to be copied (source and target name)
      virtual std::map<std::string, std::string> conditions_to_copy() const;

     protected:
      /// determine element type std::string and whether element is copied or not
      virtual bool determine_ele_type(
          Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      /// set element-specific data (material etc.)
      void set_element_data(std::shared_ptr<Core::Elements::Element> newele,
          Core::Elements::Element* oldele, const int matid, const bool isnurbs);

      /// check for correct material
      void check_material_type(const int matid);

     private:
    };  // class PoroFluidMultiPhaseCloneStrategy

  }  // namespace Utils
}  // namespace POROMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif
