/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of an exponential strain energy function with a simple active
contribution
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_ACTIVE_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_ACTIVE_HPP

#include "baci_config.hpp"

#include "baci_mixture_constituent_remodelfiber_material_exponential.hpp"

#include <Teuchos_RCPDecl.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace MAT::PAR
{
  class Material;
}  // namespace MAT::PAR

namespace MIXTURE
{
  // forward declarations
  template <typename T>
  class RemodelFiberMaterialExponentialActive;

  namespace PAR
  {
    template <typename T>
    class RemodelFiberMaterialExponentialActive : public RemodelFiberMaterial<T>
    {
      friend class MIXTURE::RemodelFiberMaterialExponentialActive<T>;

     public:
      explicit RemodelFiberMaterialExponentialActive(
          const Teuchos::RCP<MAT::PAR::Material>& matdata);

      [[nodiscard]] std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>> CreateRemodelFiberMaterial()
          const override;

      MIXTURE::PAR::ExponentialFiberParameters passive_params_;

      /// @name parameters of the exponential strain energy function
      /// @{
      const double initial_reference_density_;
      const double sigma_act_max_;
      const double lambda_act_max_;
      const double lambda_act_0_;
      const double lambda_act_;
      const double dPsiAct_;
      /// @}
    };
  }  // namespace PAR

  template <typename T>
  class RemodelFiberMaterialExponentialActive : public RemodelFiberMaterial<T>
  {
   public:
    RemodelFiberMaterialExponentialActive(
        const PAR::RemodelFiberMaterialExponentialActive<T>* matdata);

    [[nodiscard]] T GetCauchyStress(T I4) const final;

    [[nodiscard]] T GetDCauchyStressDI4(T I4) const final;

    [[nodiscard]] T GetDCauchyStressDI4DI4(T I4) const final;

   private:
    const PAR::RemodelFiberMaterialExponentialActive<T>* params_;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
