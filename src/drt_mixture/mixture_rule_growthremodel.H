/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_GROWTHREMODEL_H
#define BACI_MIXTURE_RULE_GROWTHREMODEL_H

#include "mixture_rule.H"
#include "../drt_matelast/elast_summand.H"

namespace MIXTURE
{
  // forward declaration
  class GrowthRemodelMixtureRule;

  namespace PAR
  {
    class GrowthRemodelMixtureRule : public MIXTURE::PAR::MixtureRule
    {
      friend class MIXTURE::GrowthRemodelMixtureRule;

     public:
      static constexpr int GROWTH_TYPE_ISOTROPIC = 0;
      static constexpr int GROWTH_TYPE_ANISOTROPIC = 1;

      /// constructor
      explicit GrowthRemodelMixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// default destructor
      ~GrowthRemodelMixtureRule() override = default;

      /// Create mixturerule instance
      Teuchos::RCP<MIXTURE::MixtureRule> CreateRule() override;

      /// @name parameters of the mixture rule
      /// @{

      const int growth_type_;

      /// @}
    };

  }  // namespace PAR

  /*!
   * \brief This mixture rule controls the evaluation of growth and remodel simulations with
   * homogenized constrained mixture models
   */
  class GrowthRemodelMixtureRule : public MIXTURE::MixtureRule
  {
   public:
    /// Constructor for mixture rule given the input parameters
    explicit GrowthRemodelMixtureRule(MIXTURE::PAR::GrowthRemodelMixtureRule* params);

    void PackMixtureRule(DRT::PackBuffer& data) const override;

    void UnpackMixtureRule(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    void Update(
        LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params, int gp, int eleGID) override;

    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * @brief Evaluates the inverse growth deformation gradient at the Gausspoint #gp
     *
     * The growth deformation gradient describes the deformation of the solid by addition/removal
     * of materials. It depends on the growth type.
     *
     * In case of isotropic growth (#GROWTH_TYPE_ISOTROPIC), the growth deformation gradient is
     *
     * $\boldsymbol{F}_g = (\frac{\rho_0(s)}{\rho_0(0)})^(\frac 1 3)$
     *
     * @param iFgM (out) : Inverse of the growth deformation gradient
     * @param gp (in) : Gauss point
     */
    void EvaluateInverseGrowthDeformationGradient(LINALG::Matrix<3, 3>& iFgM, int gp) const;

   private:
    ///! Rule parameters as defined in the input file
    PAR::GrowthRemodelMixtureRule* params_{};
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_RULE_GROWTHREMODEL_H