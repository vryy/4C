/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a mixture growth strategy interface

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_GROWTH_STRATEGY_H
#define MIXTURE_GROWTH_STRATEGY_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "mixture_rule.H"
#include <memory>

namespace MIXTURE
{
  class MixtureGrowthStrategy;
  namespace PAR
  {
    class MixtureGrowthStrategy : public MAT::PAR::Parameter
    {
      friend class MIXTURE::MixtureGrowthStrategy;

     public:
      /// constructor
      explicit MixtureGrowthStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// Override this method and throw error, as only the CreateGrowthStrategy() should be used.
      Teuchos::RCP<MAT::Material> CreateMaterial() final
      {
        dserror(
            "Cannot create mixture growth strategy from this method. Use CreateGrowthStrategy() "
            "instead.");
        return Teuchos::null;
      }

      /// create material instance of matching type with my parameters
      virtual std::unique_ptr<MIXTURE::MixtureGrowthStrategy> CreateGrowthStrategy() = 0;

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
    virtual void PackMixtureGrowthStrategy(DRT::PackBuffer& data) const {}

    virtual void UnpackMixtureGrowthStrategy(
        std::vector<char>::size_type& position, const std::vector<char>& data)
    {
    }

    virtual void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
    {
      // do nothing in the default case
    }

    [[nodiscard]] virtual bool HasInelasticGrowthDeformationGradient() const = 0;

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
    virtual void EvaluateInverseGrowthDeformationGradient(LINALG::Matrix<3, 3>& iFgM,
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
    virtual void EvaluateGrowthStressCmat(const MIXTURE::MixtureRule& mixtureRule,
        double currentReferenceGrowthScalar,
        const LINALG::Matrix<1, 6>& dCurrentReferenceGrowthScalarDC, const LINALG::Matrix<3, 3>& F,
        const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, const int gp,
        const int eleGID) const = 0;
  };
}  // namespace MIXTURE

#endif  // MIXTURE_GROWTH_STRATEGY_H