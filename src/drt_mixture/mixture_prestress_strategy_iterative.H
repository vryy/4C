/*----------------------------------------------------------------------*/
/*! \file

\brief General prestress strategy for mixture constituents

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H
#define BACI_MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H

#include "../drt_mat/matpar_parameter.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../linalg/linalg_fixedsizematrix.H"
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

      /// default destructor
      ~IterativePrestressStrategy() override = default;

      /// create prestress strategy instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::PrestressStrategy> CreatePrestressStrategy() override;


      /// @name parameters of the prestress strategy
      /// @{
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief General mixture constituent prestress interface to be implemented by prestressing
   * techniques.
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

    void EvaluatePrestress(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, LINALG::Matrix<3, 3>& G,
        Teuchos::ParameterList& params, int gp, int eleGID) override;

    virtual void UpdatePrestress(const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
        MIXTURE::MixtureConstituent& constituent, const LINALG::Matrix<3, 3>& F,
        LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID) override;

   private:
    /// Holder for internal parameters
    const PAR::IterativePrestressStrategy* params_;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_PRESTRESS_STRATEGY_ITERATIVE_H
