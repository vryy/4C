/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a general material for a remodel fiber constituent.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_H
#define MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_H


#include <memory>
#include "mat_par_parameter.H"

namespace MIXTURE
{
  // forward declaration
  template <typename T>
  class RemodelFiberMaterial;

  namespace PAR
  {
    template <typename T>
    class RemodelFiberMaterial : public MAT::PAR::Parameter
    {
     public:
      RemodelFiberMaterial(const Teuchos::RCP<MAT::PAR::Material>& matdata);
      Teuchos::RCP<MAT::Material> CreateMaterial() override
      {
        dserror("This type of material is not created with CreateMaterial()");
        std::exit(1);
      }

      [[nodiscard]] virtual std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
      CreateRemodelFiberMaterial() const = 0;
    };
  }  // namespace PAR

  template <typename T>
  class RemodelFiberMaterial
  {
   public:
    virtual ~RemodelFiberMaterial() = default;

    /*!
     * @brief Evaluates the strain energy function
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetPsi(T I4) const = 0;

    /*!
     * @brief Evaluates the first derivative of the strain energy function w.r.t. I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetFirstDerivativeI4(T I4) const = 0;

    /*!
     * @brief Evaluates the second derivative of the strain energy function w.r.t. I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetSecondDerivativeI4(T I4) const = 0;

    /*!
     * @brief Evaluates the third derivative of the strain energy function w.r.t. I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetThirdDerivativeI4(T I4) const = 0;
  };

}  // namespace MIXTURE

#endif  // MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_H
