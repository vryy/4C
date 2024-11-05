// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ALE_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_ALE_UTILS_CLONESTRATEGY_HPP

/*----------------------------------------------------------------------------*/
/*header inclusions */
#include "4C_config.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace Core::Elements
{
  class Element;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ALE
{
  namespace Utils
  {
    /*!
    \brief Implementation of special clone strategy for automatic generation
           of ale from a given fluid discretization

    */
    class AleCloneStrategy
    {
     public:
      /// constructor
      explicit AleCloneStrategy() {}
      /// destructor
      virtual ~AleCloneStrategy() = default;

     protected:
      /// determine element type string and whether element is copied or not
      bool determine_ele_type(Core::Elements::Element* actele,  ///< current element
          const bool ismyele,                ///< true if element belongs to my proc
          std::vector<std::string>& eletype  ///< element type
      );

      /*! \brief Set element-specific data (material etc.)
       *
       *  We need to set material and possibly other things to complete element
       *  setup. This is again really ugly as we have to extract the actual
       *  element type in order to access the material property.
       */
      void set_element_data(std::shared_ptr<Core::Elements::Element>
                                newele,     ///< newly created element where data has to be set
          Core::Elements::Element* oldele,  ///< existing element, that has been cloned
          const int matid,                  ///< ID of material law
          const bool nurbsdis               ///< Is this a Nurbs-based discretization?
      );

      /// returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> conditions_to_copy() const;

      /// check for correct material
      void check_material_type(const int matid);

     private:
    };  // class AleCloneStrategy
  }     // namespace Utils
}  // namespace ALE

FOUR_C_NAMESPACE_CLOSE

#endif
