/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H
#define MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H

#include "mat_par_parameter.H"
#include "lib_pack_buffer.H"
#include "linalg_fixedsizematrix.H"
#include "mixture_prestress_strategy.H"

// forward declarations
namespace MAT
{
  class CoordinateSystemProvider;
}

namespace MIXTURE
{
  // forward declaration
  class IterativePrestressStrategy;
  class MixtureConstituent;

  namespace PAR
  {
    class IterativePrestressStrategy : public MIXTURE::PAR::PrestressStrategy
    {
      friend class MIXTURE::IterativePrestressStrategy;

     public:
      /// constructor
      explicit IterativePrestressStrategy(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// create prestress strategy instance of matching type with my parameters
      std::unique_ptr<MIXTURE::PrestressStrategy> CreatePrestressStrategy() override;


      /// @name parameters of the prestress strategy
      /// @{

      /// Flag whether the prestretch tensor is isochoric
      const bool isochoric_;
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Mixture prestress strategy to be used with MATERIAL_ITERATIVE prestressing
   *
   * This prestressing technique has to be used with the PRESTRESSTYPE MATERIAL_ITERATIVE. In each
   * prestress update step, the internal prestretch tensor is updated with the current stretch
   * tensor of the deformation.
   */
  class IterativePrestressStrategy : public PrestressStrategy
  {
   public:
    /// Constructor for the material given the material parameters
    explicit IterativePrestressStrategy(MIXTURE::PAR::IterativePrestressStrategy* params);

    void Setup(MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void EvaluatePrestress(const MixtureRule& mixtureRule,
        const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    void UpdatePrestress(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
        CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) override;

   private:
    /// Holder for internal parameters
    const PAR::IterativePrestressStrategy* params_;
  };
}  // namespace MIXTURE

#endif  // MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H
