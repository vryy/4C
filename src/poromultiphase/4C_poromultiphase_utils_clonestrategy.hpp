/*----------------------------------------------------------------------*/
/*! \file
 \brief utils methods for cloning the porofluid discretization


   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_UTILS_CLONESTRATEGY_HPP
#define FOUR_C_POROMULTIPHASE_UTILS_CLONESTRATEGY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::Elements
{
  class Element;
}

namespace POROMULTIPHASE
{
  namespace UTILS
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
          CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype);

      /// set element-specific data (material etc.)
      void set_element_data(Teuchos::RCP<CORE::Elements::Element> newele,
          CORE::Elements::Element* oldele, const int matid, const bool isnurbs);

      /// check for correct material
      void check_material_type(const int matid);

     private:
    };  // class PoroFluidMultiPhaseCloneStrategy

  }  // namespace UTILS
}  // namespace POROMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif
