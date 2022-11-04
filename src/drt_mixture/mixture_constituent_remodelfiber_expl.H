/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel constituent with explicit update rule
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H
#define BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H

#include <Teuchos_RCPDecl.hpp>
#include <cmath>
#include "matpar_material.H"
#include "growth_evolution.H"
#include "mixture_constituent.H"
#include "mixture_constituent_remodelfiber_material.H"
#include "anisotropy_extension_default.H"
#include "remodelfiber.H"

namespace MIXTURE
{
  class MixtureConstituent;
  template <typename T>
  class RemodelFiberMaterial;

  namespace PAR
  {
    class MixtureConstituent_RemodelFiberExpl : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_RemodelFiberExpl(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      ~MixtureConstituent_RemodelFiberExpl() override = default;

      /// create material instance of matching type with my parameters
      std::unique_ptr<MIXTURE::MixtureConstituent> CreateConstituent(int id) override;

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
  class MixtureConstituent_RemodelFiberExpl : public MIXTURE::MixtureConstituent
  {
   public:
    explicit MixtureConstituent_RemodelFiberExpl(
        MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl* params, int id);

    [[nodiscard]] INPAR::MAT::MaterialType MaterialType() const override;

    void PackConstituent(DRT::PackBuffer& data) const override;

    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    void Update(
        const LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params, int gp, int eleGID) override;

    void UpdateElasticPart(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFext,
        Teuchos::ParameterList& params, double dt, int gp, int eleGID) override;

    void UpdatePrestress(
        const LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params, int gp, int eleGID) override;

    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    void EvaluateElasticPart(const LINALG::Matrix<3, 3>& FM, const LINALG::Matrix<3, 3>& iFextin,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    [[nodiscard]] double GetGrowthScalar(int gp) const override;

    void RegisterVtkOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateVtkOutputData(
        const std::string& name, Epetra_SerialDenseMatrix& data) const override;

   private:
    [[nodiscard]] double EvaluateLambdaf(const LINALG::Matrix<3, 3>& C, int gp, int eleGID) const;
    [[nodiscard]] double EvaluateLambdaExt(
        const LINALG::Matrix<3, 3>& iFext, int gp, int eleGID) const;

    [[nodiscard]] LINALG::Matrix<6, 1> EvaluateCurrentPK2(int gp, int eleGID) const;
    [[nodiscard]] LINALG::Matrix<6, 6> EvaluateCurrentCmat(int gp, int eleGID) const;

    [[nodiscard]] double EvaluateDepositionStretch(double time) const;
    void UpdateHomeostaticValues(Teuchos::ParameterList& params, int eleGID);

    void Initialize();

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl* params_;

    /// An instance of the remodel fiber
    std::vector<RemodelFiber<2>> remodel_fiber_;

    /// Handler for anisotropic input
    MAT::DefaultAnisotropyExtension<1> anisotropy_extension_;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_REMODELFIBER_EXPL_H
