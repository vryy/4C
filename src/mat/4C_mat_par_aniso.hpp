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

namespace Mat
{
  namespace Elastic
  {
    class StructuralTensorStrategyBase;
  }
}  // namespace Mat

/*----------------------------------------------------------------------*/
/* declarations */
namespace Mat
{
  namespace PAR
  {
    /// Extension to hold 'quick' access material parameters for anisotropy
    class ParameterAniso : public Core::Mat::PAR::Parameter
    {
     public:
      /// construct the material object given material parameters
      ParameterAniso(Teuchos::RCP<const Core::Mat::PAR::Material> matdata);

      /// return pointer to strategy
      const Teuchos::RCP<Mat::Elastic::StructuralTensorStrategyBase>& structural_tensor_strategy()
      {
        return structural_tensor_strategy_;
      };

     private:
      /// structural tensor strategy
      Teuchos::RCP<Mat::Elastic::StructuralTensorStrategyBase> structural_tensor_strategy_;

    };  // class ParameterAniso

  }  // namespace PAR

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif
