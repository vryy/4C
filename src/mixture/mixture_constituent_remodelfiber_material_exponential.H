/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel fiber material with exponential strain energy function.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H
#define MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H

#include <Teuchos_RCPDecl.hpp>
#include <memory>
#include "mixture_constituent_remodelfiber_material.H"

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

      const double k1_;
      const double k2_;
      const bool supports_compression_;

      /// @}
    };
  }  // namespace PAR

  template <typename T>
  class RemodelFiberMaterialExponential : public RemodelFiberMaterial<T>
  {
   public:
    RemodelFiberMaterialExponential(const PAR::RemodelFiberMaterialExponential<T>* matdata);


    [[nodiscard]] T GetPsi(T I4) const override;


    [[nodiscard]] T GetFirstDerivativeI4(T I4) const override;


    [[nodiscard]] T GetSecondDerivativeI4(T I4) const override;


    [[nodiscard]] T GetThirdDerivativeI4(T I4) const override;

   private:
    const PAR::RemodelFiberMaterialExponential<T>* params_;
  };

}  // namespace MIXTURE

#endif  // MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_EXPONENTIAL_H
