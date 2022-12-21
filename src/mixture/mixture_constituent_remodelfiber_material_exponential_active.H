/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of an exponential strain energy function with a simple active
contribution
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_ACTIVE_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_ACTIVE_H

#include <Teuchos_RCPDecl.hpp>
#include <memory>
#include "mixture_constituent_remodelfiber_material_exponential.H"
#include "mixture_constituent_remodelfiber_material_active.H"

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
    class RemodelFiberMaterialExponentialActive : public RemodelFiberMaterialExponential<T>
    {
      friend class MIXTURE::RemodelFiberMaterialExponentialActive<T>;

     public:
      explicit RemodelFiberMaterialExponentialActive(
          const Teuchos::RCP<MAT::PAR::Material>& matdata);

      [[nodiscard]] std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>> CreateRemodelFiberMaterial()
          const override;

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
  class RemodelFiberMaterialExponentialActive : public RemodelFiberMaterialExponential<T>,
                                                public RemodelFiberMaterialActive<T>
  {
   public:
    RemodelFiberMaterialExponentialActive(
        const PAR::RemodelFiberMaterialExponentialActive<T>* matdata);


    [[nodiscard]] T GetPsi(T I4) const override;


    [[nodiscard]] T GetFirstDerivativeActive() const override;

   private:
    const PAR::RemodelFiberMaterialExponentialActive<T>* params_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_ACTIVE_H
