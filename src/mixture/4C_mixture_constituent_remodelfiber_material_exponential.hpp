/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_HPP

#include "4C_config.hpp"

#include "4C_mixture_constituent_remodelfiber_lib.hpp"
#include "4C_mixture_constituent_remodelfiber_material.hpp"

#include <Teuchos_RCP.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

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
      explicit RemodelFiberMaterialExponential(const Core::Mat::PAR::Parameter::Data& matdata);

      [[nodiscard]] std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
      create_remodel_fiber_material() const override;

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

    [[nodiscard]] T get_d_cauchy_stress_d_i4_d_i4(T I4) const final;

   private:
    const PAR::RemodelFiberMaterialExponential<T>* params_;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
