/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a mixture growth strategy interface

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_GROWTH_STRATEGY_HPP
#define FOUR_C_MIXTURE_GROWTH_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mixture_rule.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  class MixtureGrowthStrategy;
  namespace PAR
  {
    class MixtureGrowthStrategy : public CORE::MAT::PAR::Parameter
    {
      friend class MIXTURE::MixtureGrowthStrategy;

     public:
      /// constructor
      explicit MixtureGrowthStrategy(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

      /// Override this method and throw error, as only the create_growth_strategy() should be used.
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() final
      {
        FOUR_C_THROW(
            "Cannot create mixture growth strategy from this method. Use create_growth_strategy() "
            "instead.");
        return Teuchos::null;
      }

      /// create material instance of matching type with my parameters
      virtual std::unique_ptr<MIXTURE::MixtureGrowthStrategy> create_growth_strategy() = 0;

      /*!
       * \brief Factory of the mixture growth strategy parameters
       *
       * This static method generates the specific class of the mixture growth strategy defined in
       * the datfile at the corresponding material id
       *
       * @param matid Material id of the mixturerule
       * @return Parameters of the referenced mixture rule
       */
      static MIXTURE::PAR::MixtureGrowthStrategy* Factory(int matid);
    };
  }  // namespace PAR

  class MixtureGrowthStrategy
  {
   public:
    virtual ~MixtureGrowthStrategy() = default;

    MixtureGrowthStrategy() = default;
    MixtureGrowthStrategy(const MixtureGrowthStrategy& copyFrom) = default;
    MixtureGrowthStrategy& operator=(const MixtureGrowthStrategy& copyFrom) = default;
    MixtureGrowthStrategy(MixtureGrowthStrategy&&) noexcept = default;
    MixtureGrowthStrategy& operator=(MixtureGrowthStrategy&&) noexcept = default;

    virtual void pack_mixture_growth_strategy(CORE::COMM::PackBuffer& data) const {}

    virtual void unpack_mixture_growth_strategy(
        std::vector<char>::size_type& position, const std::vector<char>& data)
    {
    }

    virtual void register_anisotropy_extensions(MAT::Anisotropy& anisotropy)
    {
      // do nothing in the default case
    }

    [[nodiscard]] virtual bool has_inelastic_growth_deformation_gradient() const = 0;

    /*!
     * @brief Evaluates the inverse growth deformation gradient at the Gausspoint #gp
     *
     * The growth deformation gradient describes the deformation of the solid by addition/removal
     * of materials.
     *
     * @param iFgM (out) : Inverse of the growth deformation gradient
     * @param mixtureRule (in) : mixture rule
     * @param currentReferenceGrowthScalar (in) : current reference growth scalar
     * @param gp (in) : Gauss point
     */
    virtual void evaluate_inverse_growth_deformation_gradient(CORE::LINALG::Matrix<3, 3>& iFgM,
        const MIXTURE::MixtureRule& mixtureRule, double currentReferenceGrowthScalar,
        int gp) const = 0;

    /*!
     * @brief Evaluates the contribution of the growth strategy to the stress tensor and the
     * linearization.
     *
     * This is meant for growth strategies that use some kind of penalty formulation to ensure
     * growth
     *
     * @param mixtureRule (in) : mixture rule
     * @param currentReferenceGrowthScalar (in) : current reference growth scalar (volume change in
     * percent)
     * @param dCurrentReferenceGrowthScalarDC (in) : Derivative of the current reference growth
     * scalar w.r.t. Cauchy green deformation tensor
     * @param F (in) : deformation gradient
     * @param E_strain (in) : Green-Langrange strain tensor
     * @param params (in) : Container for additional information
     * @param S_stress (out) : 2nd Piola-Kirchhoff stress tensor in stress like Voigt notation
     * @param cmat (out) : linearization of the 2nd Piola-Kirchhoff stress tensor
     * @param gp (in) : Gauss point
     * @param eleGID (in) : global element id
     */
    virtual void evaluate_growth_stress_cmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar,
        const CORE::LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC,
        const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, const int gp, const int eleGID) const = 0;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif