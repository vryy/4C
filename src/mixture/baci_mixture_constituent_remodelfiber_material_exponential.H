/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H

#include "baci_mixture_constituent_remodelfiber_lib.H"
#include "baci_mixture_constituent_remodelfiber_material.H"

#include <Teuchos_RCP.hpp>

#include <memory>

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

#endif  // MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H
