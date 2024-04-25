/*----------------------------------------------------------------------*/
/*! \file

\brief Constant prestretch strategy

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_PRESTRESS_STRATEGY_CONSTANT_HPP
#define FOUR_C_MIXTURE_PRESTRESS_STRATEGY_CONSTANT_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mixture_elastin_membrane_prestress_strategy.hpp"
#include "4C_mixture_prestress_strategy.hpp"

#include <NOX_LAPACK.H>

#include <array>

FOUR_C_NAMESPACE_OPEN

namespace MIXTURE
{
  // forward declaration
  class ConstantPrestressStrategy;

  namespace PAR
  {
    class ConstantPrestressStrategy : public MIXTURE::PAR::PrestressStrategy
    {
      friend class MIXTURE::ConstantPrestressStrategy;

     public:
      /// constructor
      explicit ConstantPrestressStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// create prestress strategy instance of matching type with my parameters
      std::unique_ptr<MIXTURE::PrestressStrategy> CreatePrestressStrategy() override;

      /// @name parameters of the prestress strategy
      /// @{
      std::array<double, 9> prestretch_;
      /// @}
    };
  }  // namespace PAR


  /*!
   * \brief Prestressing strategy for a constant predefined prestretch tensor.
   */
  class ConstantPrestressStrategy : public PrestressStrategy
  {
   public:
    /// Constructor for the material given the material parameters
    explicit ConstantPrestressStrategy(MIXTURE::PAR::ConstantPrestressStrategy* params);

    void Setup(MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void EvaluatePrestress(const MixtureRule& mixtureRule,
        const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
        MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    void Update(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
        CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) override;

   private:
    /// Holder for internal parameters
    const PAR::ConstantPrestressStrategy* params_;
  };
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif
