/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_RULE_GROWTHREMODEL_HPP
#define FOUR_C_MIXTURE_RULE_GROWTHREMODEL_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mixture_rule.hpp"

#include <Teuchos_RCPDecl.hpp>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace LINALG
{
  class SerialDenseMatrix;
}
namespace CORE::COMM
{
  class PackBuffer;
}
namespace MAT
{
  namespace PAR
  {
    class Material;
  }
}  // namespace MAT

namespace MIXTURE
{
  // forward declaration
  class GrowthRemodelMixtureRule;
  class MixtureConstituent;
  class MixtureGrowthStrategy;

  namespace PAR
  {
    class GrowthRemodelMixtureRule : public MIXTURE::PAR::MixtureRule
    {
      friend class MIXTURE::GrowthRemodelMixtureRule;

     public:
      static constexpr int GROWTH_TYPE_ISOTROPIC = 0;
      static constexpr int GROWTH_TYPE_ANISOTROPIC = 1;

      /// constructor
      explicit GrowthRemodelMixtureRule(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

      /// Create mixturerule instance
      std::unique_ptr<MIXTURE::MixtureRule> CreateRule() override;

      /// @name parameters of the mixture rule
      /// @{

      const int growth_strategy_matid_;

      /// Initial reference density of the whole material
      const double initial_reference_density_;

      /// Initial mass fractions
      const std::vector<double> mass_fractions_;

      /// @}
    };

  }  // namespace PAR

  /*!
   * \brief This mixture rule controls the evaluation of growth and remodel simulations with
   * homogenized constrained mixture models
   */
  class GrowthRemodelMixtureRule : public MIXTURE::MixtureRule
  {
   private:
    static constexpr auto OUTPUT_CURRENT_REFERENCE_DENSITY = "current_reference_density";

   public:
    /// Constructor for mixture rule given the input parameters
    explicit GrowthRemodelMixtureRule(MIXTURE::PAR::GrowthRemodelMixtureRule* params);

    void PackMixtureRule(CORE::COMM::PackBuffer& data) const override;

    void UnpackMixtureRule(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void register_anisotropy_extensions(MAT::Anisotropy& anisotropy) override;

    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    void Update(CORE::LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    /*!
     * \brief Returns the initial reference mass denstity of the constituent
     *
     * \param constituent
     * \return double
     */
    double get_constituent_initial_reference_mass_density(
        const MIXTURE::MixtureConstituent& constituent) const;

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateOutputData(
        const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const override;

   private:
    double compute_current_reference_growth_scalar(int gp) const;

    ///! Rule parameters as defined in the input file
    PAR::GrowthRemodelMixtureRule* params_{};

    std::unique_ptr<MIXTURE::MixtureGrowthStrategy> growth_strategy_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif