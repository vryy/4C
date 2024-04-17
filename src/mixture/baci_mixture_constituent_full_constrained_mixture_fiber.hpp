/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel constituent with implicit update rule
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_FULL_CONSTRAINED_MIXTURE_FIBER_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_FULL_CONSTRAINED_MIXTURE_FIBER_HPP

#include "baci_config.hpp"

#include "baci_mat_anisotropy_extension_default.hpp"
#include "baci_mat_par_material.hpp"
#include "baci_mixture_constituent.hpp"
#include "baci_mixture_constituent_remodelfiber_material.hpp"
#include "baci_mixture_full_constrained_mixture_fiber.hpp"

#include <Teuchos_RCP.hpp>

#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  class MixtureConstituent;
  template <typename T>
  class FullConstrainedMixtureFiber;

  namespace PAR
  {
    class MixtureConstituent_FullConstrainedMixtureFiber : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_FullConstrainedMixtureFiber(
          const Teuchos::RCP<MAT::PAR::Material>& matdata);
      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

      const int fiber_id_;
      const int init_;

      const int fiber_material_id_;
      const MIXTURE::PAR::RemodelFiberMaterial<double>* fiber_material_;

      const bool growth_enabled_;
      const double poisson_decay_time_;
      const double growth_constant_;

      const double deposition_stretch_;
      const int initial_deposition_stretch_timefunc_num_;

      const HistoryAdaptionStrategy adaptive_history_strategy_;
      const double adaptive_history_tolerance_;
    };
  }  // namespace PAR

  /*!
   * \brief Full constrained mixture fiber constituent
   */
  class MixtureConstituent_FullConstrainedMixtureFiber : public MIXTURE::MixtureConstituent
  {
   public:
    MixtureConstituent_FullConstrainedMixtureFiber(
        MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber* params, int id);

    [[nodiscard]] INPAR::MAT::MaterialType MaterialType() const override;

    void PackConstituent(CORE::COMM::PackBuffer& data) const override;

    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    void ReadElement(int numgp, INPUT::LineDefinition* linedef) override;

    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    void Update(const CORE::LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    void EvaluateElasticPart(const CORE::LINALG::Matrix<3, 3>& FM,
        const CORE::LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) override;

    [[nodiscard]] double GetGrowthScalar(int gp) const override;
    [[nodiscard]] CORE::LINALG::Matrix<1, 6> GetDGrowthScalarDC(int gp, int eleGID) const override;

    void RegisterOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateOutputData(
        const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const override;

   private:
    [[nodiscard]] double EvaluateLambdaf(
        const CORE::LINALG::Matrix<3, 3>& C, int gp, int eleGID) const;
    [[nodiscard]] CORE::LINALG::Matrix<1, 6> EvaluateDLambdafsqDC(int gp, int eleGID) const;

    [[nodiscard]] CORE::LINALG::Matrix<6, 1> EvaluateCurrentPK2(int gp, int eleGID) const;
    [[nodiscard]] CORE::LINALG::Matrix<6, 6> EvaluateCurrentCmat(int gp, int eleGID) const;
    [[nodiscard]] double EvaluateInitialDepositionStretch(double time) const;

    void Initialize();

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber* params_;

    /// An instance of the full constrained mixture fiber for each Gauss point
    std::vector<FullConstrainedMixtureFiber<double>> full_constrained_mixture_fiber_;

    /// Store the last converged lambda_f for initializing history of full constrained mixture model
    std::vector<double> last_lambda_f_;

    /// Handler for anisotropic input
    MAT::DefaultAnisotropyExtension<1> anisotropy_extension_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
