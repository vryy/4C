/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for growth and remodeling simulations with homogenized constrained mixtures

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_GROWTHREMODEL_H
#define BACI_MIXTURE_RULE_GROWTHREMODEL_H

#include "mixture_constituent.H"
#include "mixture_rule.H"
#include <Teuchos_RCPDecl.hpp>
#include <string>
#include <unordered_map>
#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
class Epetra_SerialDenseMatrix;
namespace DRT
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
namespace Teuchos
{
  class ParameterList;
}

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
     * \brief Returns the initial reference mass denstity of the constituent
     *
     * \param constituent
     * \return double
     */
    double GetConstituentInitialReferenceMassDensity(
        const MIXTURE::MixtureConstituent& constituent) const;

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

    void RegisterVtkOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool EvaluateVtkOutputData(
        const std::string& name, Epetra_SerialDenseMatrix& data) const override;

   private:
    double ComputeCurrentReferenceDensity(int gp) const;
    ///! Rule parameters as defined in the input file
    PAR::GrowthRemodelMixtureRule* params_{};
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_RULE_GROWTHREMODEL_H