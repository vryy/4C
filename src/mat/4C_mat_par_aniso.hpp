/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a base class for anisotropic material parameters

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef FOUR_C_MAT_PAR_ANISO_HPP
#define FOUR_C_MAT_PAR_ANISO_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_material_parameter_base.hpp"


FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    class StructuralTensorStrategyBase;
  }
}  // namespace MAT

/*----------------------------------------------------------------------*/
/* declarations */
namespace MAT
{
  namespace PAR
  {
    /// Extension to hold 'quick' access material parameters for anisotropy
    class ParameterAniso : public CORE::MAT::PAR::Parameter
    {
     public:
      /// construct the material object given material parameters
      ParameterAniso(Teuchos::RCP<const CORE::MAT::PAR::Material> matdata);

      /// return pointer to strategy
      const Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase>& StructuralTensorStrategy()
      {
        return structural_tensor_strategy_;
      };

     private:
      /// structural tensor strategy
      Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase> structural_tensor_strategy_;

    };  // class ParameterAniso

  }  // namespace PAR

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif
