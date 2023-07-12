/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a remodel constituent with implicit update rule
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_CONSTITUENT_REMODELFIBER_IMPL_H
#define MIXTURE_CONSTITUENT_REMODELFIBER_IMPL_H

#include <Teuchos_RCPDecl.hpp>
#include <cmath>
#include "mat_par_material.H"
#include "mixture_constituent.H"
#include "mixture_constituent_remodelfiber_material.H"
#include "mat_anisotropy_extension_default.H"
#include "mixture_remodelfiber.H"

namespace MIXTURE
{
  class MixtureConstituent;
  template <typename T>
  class RemodelFiberMaterial;

  namespace PAR
  {
    class MixtureConstituent_RemodelFiberImpl : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_RemodelFiberImpl(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      ~MixtureConstituent_RemodelFiberImpl() override = default;

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
      const int deposition_stretch_timefunc_num_;
    };
  }  // namespace PAR

  /*!
   * \brief Remodel fiber constituent with an implicit update rule
   */
  class MixtureConstituent_RemodelFiberImpl : public MIXTURE::MixtureConstituent
  {
   public:
    explicit MixtureConstituent_RemodelFiberImpl(
        MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl* params, int id);

    /// Returns the material type enum
    [[nodiscard]] INPAR::MAT::MaterialType MaterialType() const override;

    void PackConstituent(DRT::PackBuffer& data) const override;

    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    void Update(const CORE::LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void UpdatePrestress(const CORE::LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params,
        int gp, int eleGID) override;

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    void EvaluateElasticPart(const CORE::LINALG::Matrix<3, 3>& FM,
        const CORE::LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) override;

    [[nodiscard]] double GetGrowthScalar(int gp) const override;
    [[nodiscard]] CORE::LINALG::Matrix<1, 6> GetDGrowthScalarDC(int gp, int eleGID) const override;

    void RegisterVtkOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateVtkOutputData(
        const std::string& name, Epetra_SerialDenseMatrix& data) const override;

   private:
    void IntegrateLocalEvolutionEquations(double dt, int gp, int eleGID);
    [[nodiscard]] double EvaluateLambdaf(
        const CORE::LINALG::Matrix<3, 3>& C, int gp, int eleGID) const;
    [[nodiscard]] CORE::LINALG::Matrix<1, 6> EvaluateDLambdafsqDC(int gp, int eleGID) const;

    [[nodiscard]] CORE::LINALG::Matrix<6, 1> EvaluateCurrentPK2(int gp, int eleGID) const;
    [[nodiscard]] CORE::LINALG::Matrix<6, 6> EvaluateCurrentCmat(int gp, int eleGID) const;

    [[nodiscard]] double EvaluateDepositionStretch(double time) const;
    void UpdateHomeostaticValues(const Teuchos::ParameterList& params, int eleGID);

    void Initialize();

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl* params_;

    /// An instance of the remodel fiber
    std::vector<RemodelFiber<2>> remodel_fiber_;

    /// Handler for anisotropic input
    MAT::DefaultAnisotropyExtension<1> anisotropy_extension_;

    /// temporary variables used within one evaluation
    /// @{
    std::vector<CORE::LINALG::Matrix<1, 6>> dgrowthscalardC_;
    std::vector<CORE::LINALG::Matrix<1, 6>> dlambdardC_;
    /// @}
  };
}  // namespace MIXTURE

#endif  // MIXTURE_CONSTITUENT_REMODELFIBER_IMPL_H
