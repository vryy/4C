/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone scatra discretization from elemag discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_ELEMAG_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_ELEMAG_UTILS_CLONESTRATEGY_HPP

/*----------------------------------------------------------------------------*/
/*header inclusions */
#include "4C_config.hpp"

#include "4C_discretization_fem_general_shape_function_type.hpp"

#include <Teuchos_RCP.hpp>

#include <map>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{
  class Element;
}  // namespace DRT

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ELEMAG
{
  namespace UTILS
  {
    /*!
    \brief Implementation of special clone strategy for automatic generation
           of scatra from a given elemag discretization

    */
    template <CORE::FE::ShapeFunctionType sft>
    class ScatraCloneStrategy
    {
     public:
      /// constructor
      explicit ScatraCloneStrategy() {}
      /// destructor
      virtual ~ScatraCloneStrategy() = default;

     protected:
      /// determine element type string and whether element is copied or not
      bool determine_ele_type(DRT::Element* actele,  ///< current element
          const bool ismyele,                        ///< true if element belongs to my proc
          std::vector<std::string>& eletype          ///< element type
      );

      /*! \brief Set element-specific data (material etc.)
       *
       *  We need to set material and possibly other things to complete element
       *  setup. This is again really ugly as we have to extract the actual
       *  element type in order to access the material property.
       */
      void set_element_data(
          Teuchos::RCP<DRT::Element> newele,  ///< newly created element where data has to be set
          DRT::Element* oldele,               ///< existing element, that has been cloned
          const int matid,                    ///< ID of material law
          const bool nurbsdis                 ///< Is this a Nurbs-based discretization?
      );

      /// returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> conditions_to_copy() const;

      /// check for correct material
      void check_material_type(const int matid);
    };  // class ScatraCloneStrategy
  }     // namespace UTILS
}  // namespace ELEMAG

FOUR_C_NAMESPACE_CLOSE

#endif
