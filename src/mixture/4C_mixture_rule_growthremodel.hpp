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

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace LinAlg
{
  class SerialDenseMatrix;
}
namespace Core::Communication
{
  class PackBuffer;
}
namespace Mat
{
  namespace PAR
  {
    class Material;
  }
}  // namespace Mat

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
      explicit GrowthRemodelMixtureRule(const Core::Mat::PAR::Parameter::Data& matdata);

      /// Create mixturerule instance
      std::unique_ptr<MIXTURE::MixtureRule> create_rule() override;

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

    void pack_mixture_rule(Core::Communication::PackBuffer& data) const override;

    void unpack_mixture_rule(Core::Communication::UnpackBuffer& buffer) override;

    void register_anisotropy_extensions(Mat::Anisotropy& anisotropy) override;

    void setup(Teuchos::ParameterList& params, int eleGID) override;

    void update(Core::LinAlg::Matrix<3, 3> const& F, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void evaluate(const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
        Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    /*!
     * \brief Returns the initial reference mass denstity of the constituent
     *
     * \param constituent
     * \return double
     */
    [[nodiscard]] double get_constituent_initial_reference_mass_density(
        const MIXTURE::MixtureConstituent& constituent) const;

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

   private:
    [[nodiscard]] double compute_current_reference_growth_scalar(int gp) const;

    ///! Rule parameters as defined in the input file
    PAR::GrowthRemodelMixtureRule* params_{};

    std::unique_ptr<MIXTURE::MixtureGrowthStrategy> growth_strategy_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif