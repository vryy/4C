/*----------------------------------------------------------------------*/
/*! \file

\brief Constant prestretch strategy

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_PRESTRESS_STRATEGY_CONSTANT_H
#define BACI_MIXTURE_PRESTRESS_STRATEGY_CONSTANT_H

#include <NOX_LAPACK.H>
#include "../linalg/linalg_fixedsizematrix.H"
#include "mixture_prestress_strategy.H"
#include "elastin_membrane_prestress_strategy.H"

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

      /// default destructor
      ~ConstantPrestressStrategy() override = default;

      /// create prestress strategy instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::PrestressStrategy> CreatePrestressStrategy() override;

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

    void EvaluatePrestress(const MixtureRule& mixtureRule,
        const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
        MIXTURE::MixtureConstituent& constituent, LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    void UpdatePrestress(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const LINALG::Matrix<3, 3>& F,
        LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) override;

   private:
    /// Holder for internal parameters
    const PAR::ConstantPrestressStrategy* params_;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_PRESTRESS_STRATEGY_CONSTANT_H
