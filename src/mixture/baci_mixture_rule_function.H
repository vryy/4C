/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule using functions to define the massfractions of the mixture materials

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_FUNCTION_H
#define BACI_MIXTURE_RULE_FUNCTION_H

#include "baci_config.H"

#include "baci_linalg_fixedsizematrix.H"
#include "baci_mixture_rule.H"
#include "baci_utils_function.H"

#include <Teuchos_RCPDecl.hpp>

#include <memory>
#include <vector>

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

// forward declarations
namespace MAT::PAR
{
  class Material;
}

namespace MIXTURE
{
  // forward declaration
  class FunctionMixtureRule;

  namespace PAR
  {
    class FunctionMixtureRule : public MIXTURE::PAR::MixtureRule
    {
      friend class MIXTURE::FunctionMixtureRule;

     public:
      /// constructor
      explicit FunctionMixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// Create mixturerule instance
      std::unique_ptr<MIXTURE::MixtureRule> CreateRule() override;

      /// @name parameters of the mixture rule
      /// @{
      const double initial_reference_density_;

      const std::vector<int> mass_fractions_funct_ids_;
      /// @}
    };

  }  // namespace PAR

  /*!
   * \brief Mixture rule to be used in homogenized constrained mixture models. It scales the mass
   * fractions of the individual constitutents by functions of space and time.
   */
  class FunctionMixtureRule : public MIXTURE::MixtureRule
  {
   public:
    /// Constructor for mixture rule given the input parameters
    explicit FunctionMixtureRule(MIXTURE::PAR::FunctionMixtureRule* params);

    void Evaluate(const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
        CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) override;

    [[nodiscard]] double ReturnMassDensity() const override
    {
      return params_->initial_reference_density_;
    };

    void Setup(Teuchos::ParameterList& params, const int eleGID) override;

    void UnpackMixtureRule(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

   private:
    ///! Rule parameters as defined in the input file
    PAR::FunctionMixtureRule* params_{};
    std::vector<const CORE::UTILS::FunctionOfSpaceTime*> mass_fractions_functions_;
  };
}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif  // MIXTURE_RULE_FUNCTION_H