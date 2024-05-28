/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a general material for a remodel fiber constituent.
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_MATERIAL_HPP


#include "4C_config.hpp"

#include "4C_material_parameter_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  // forward declaration
  template <typename T>
  class RemodelFiberMaterial;

  namespace PAR
  {
    template <typename T>
    class RemodelFiberMaterial : public CORE::MAT::PAR::Parameter
    {
     public:
      RemodelFiberMaterial(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);
      Teuchos::RCP<CORE::MAT::Material> create_material() override
      {
        FOUR_C_THROW("This type of material is not created with create_material()");
        std::exit(1);
      }

      [[nodiscard]] virtual std::unique_ptr<MIXTURE::RemodelFiberMaterial<T>>
      create_remodel_fiber_material() const = 0;
    };
  }  // namespace PAR

  template <typename T>
  class RemodelFiberMaterial
  {
   public:
    virtual ~RemodelFiberMaterial() = default;

    /*!
     * @brief Evaluates the Cauchy stress as a function of I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetCauchyStress(T I4) const = 0;

    /*!
     * @brief Evaluates the first derivative of the Cauchy stress w.r.t. I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T GetDCauchyStressDI4(T I4) const = 0;

    /*!
     * @brief Evaluates the second derivative of the Cauchy stress w.r.t. I4
     *
     * @param I4 Fourth invariant of the Cauchy-Green tensor
     * @return T
     */
    [[nodiscard]] virtual T get_d_cauchy_stress_d_i4_d_i4(T I4) const = 0;
  };

}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
