/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_HPP

#include "baci_config.hpp"

#include "baci_mixture_constituent_remodelfiber_lib.hpp"
#include "baci_mixture_constituent_remodelfiber_material.hpp"

#include <Teuchos_RCP.hpp>

#include <memory>

BACI_NAMESPACE_OPEN

// forward declarations
namespace MAT::PAR
{
  class Material;
}  // namespace MAT::PAR

namespace MIXTURE
{
  // forward declarations
  template <typename T>
  class RemodelFiberMaterialExponential;

  namespace PAR
  {
    template <typename T>
    class RemodelFiberMaterialExponential : public RemodelFiberMaterial<T>
    {
      friend class MIXTURE::RemodelFiberMaterialExponential<T>;

     public:
      explicit RemodelFiberMaterialExponential(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      [[nodiscard]] std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>> CreateRemodelFiberMaterial()
          const override;

      /// @name parameters of the exponential strain energy function
      /// @{
      MIXTURE::PAR::ExponentialFiberParameters params_;
      /// @}
    };
  }  // namespace PAR

  template <typename T>
  class RemodelFiberMaterialExponential : public RemodelFiberMaterial<T>
  {
   public:
    RemodelFiberMaterialExponential(const PAR::RemodelFiberMaterialExponential<T>* matdata);

    [[nodiscard]] T GetCauchyStress(T I4) const final;

    [[nodiscard]] T GetDCauchyStressDI4(T I4) const final;

    [[nodiscard]] T GetDCauchyStressDI4DI4(T I4) const final;

   private:
    const PAR::RemodelFiberMaterialExponential<T>* params_;
  };

}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif
