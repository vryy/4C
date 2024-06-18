/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel constituent with explicit update rule
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_HPP
#define FOUR_C_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_HPP

#include "4C_config.hpp"

#include "4C_mat_anisotropy_extension_default.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_mixture_constituent_remodelfiber_material.hpp"
#include "4C_mixture_remodelfiber.hpp"

#include <Teuchos_RCPDecl.hpp>

#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  class MixtureConstituent;
  template <typename T>
  class RemodelFiberMaterial;

  namespace PAR
  {
    class MixtureConstituentRemodelFiberExpl : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituentRemodelFiberExpl(const Core::Mat::PAR::Parameter::Data& matdata);
      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> create_constituent(int id) override;

      const int fiber_id_;
      const int init_;
      const double gamma_;

      const int fiber_material_id_;
      const MIXTURE::PAR::RemodelFiberMaterial<double>* fiber_material_;

      const bool growth_enabled_;
      const double poisson_decay_time_;
      const double growth_constant_;

      const double deposition_stretch_;
      const int deposition_stretch_timefunc_num_;

      const bool inelastic_external_deformation_;
    };
  }  // namespace PAR

  /*!
   * \brief Remodel fiber constituent with an explicit update rule
   */
  class MixtureConstituentRemodelFiberExpl : public MIXTURE::MixtureConstituent
  {
   public:
    explicit MixtureConstituentRemodelFiberExpl(
        MIXTURE::PAR::MixtureConstituentRemodelFiberExpl* params, int id);

    [[nodiscard]] Core::Materials::MaterialType material_type() const override;

    void pack_constituent(Core::Communication::PackBuffer& data) const override;

    void unpack_constituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void register_anisotropy_extensions(Mat::Anisotropy& anisotropy) override;

    void read_element(int numgp, Input::LineDefinition* linedef) override;

    void setup(Teuchos::ParameterList& params, int eleGID) override;

    void update(const Core::LinAlg::Matrix<3, 3>& F, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void update_elastic_part(const Core::LinAlg::Matrix<3, 3>& F,
        const Core::LinAlg::Matrix<3, 3>& iFext, Teuchos::ParameterList& params, double dt, int gp,
        int eleGID) override;

    void evaluate(const Core::LinAlg::Matrix<3, 3>& F, const Core::LinAlg::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
        Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    void evaluate_elastic_part(const Core::LinAlg::Matrix<3, 3>& FM,
        const Core::LinAlg::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, int gp,
        int eleGID) override;

    [[nodiscard]] double get_growth_scalar(int gp) const override;

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

   private:
    [[nodiscard]] double evaluate_lambdaf(
        const Core::LinAlg::Matrix<3, 3>& C, int gp, int eleGID) const;
    [[nodiscard]] double evaluate_lambda_ext(
        const Core::LinAlg::Matrix<3, 3>& iFext, int gp, int eleGID) const;

    [[nodiscard]] Core::LinAlg::Matrix<6, 1> evaluate_current_p_k2(int gp, int eleGID) const;
    [[nodiscard]] Core::LinAlg::Matrix<6, 6> evaluate_current_cmat(int gp, int eleGID) const;

    [[nodiscard]] double evaluate_deposition_stretch(double time) const;
    void update_homeostatic_values(const Teuchos::ParameterList& params, int eleGID);

    void initialize();

    /// my material parameters
    MIXTURE::PAR::MixtureConstituentRemodelFiberExpl* params_;

    /// An instance of the remodel fiber
    std::vector<RemodelFiber<2>> remodel_fiber_;

    /// Handler for anisotropic input
    Mat::DefaultAnisotropyExtension<1> anisotropy_extension_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
