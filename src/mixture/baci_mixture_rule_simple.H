/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for very simple mixture materials

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_SIMPLE_H
#define BACI_MIXTURE_RULE_SIMPLE_H

#include "baci_config.H"

#include "baci_linalg_fixedsizematrix.H"
#include "baci_mixture_rule.H"

#include <Teuchos_RCPDecl.hpp>

#include <memory>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

// forward declarations
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
  class SimpleMixtureRule;

  namespace PAR
  {
    class SimpleMixtureRule : public MIXTURE::PAR::MixtureRule
    {
      friend class MIXTURE::SimpleMixtureRule;

     public:
      /// constructor
      explicit SimpleMixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// Create mixturerule instance
      std::unique_ptr<MIXTURE::MixtureRule> CreateRule() override;

      /// @name parameters of the mixture rule
      /// @{
      const double initial_reference_density_;

      const std::vector<double> mass_fractions_;
      /// @}
    };

  }  // namespace PAR

  /*!
   * \brief This mixture rule controls the evaluation of growth and remodel simulations with
   * homogenized constrained mixture models
   */
  class SimpleMixtureRule : public MIXTURE::MixtureRule
  {
   public:
    /// Constructor for mixture rule given the input parameters
    explicit SimpleMixtureRule(MIXTURE::PAR::SimpleMixtureRule* params);

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    double ReturnMassDensity() const override { return params_->initial_reference_density_; };

   private:
    ///! Rule parameters as defined in the input file
    PAR::SimpleMixtureRule* params_{};
  };
}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif  // MIXTURE_RULE_SIMPLE_H