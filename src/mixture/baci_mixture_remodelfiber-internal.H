/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a 1D remodel fiber implementation
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_REMODELFIBER_INTERNAL_H
#define BACI_MIXTURE_REMODELFIBER_INTERNAL_H

#include "baci_linalg_fixedsizematrix.H"

#include <Sacado_tradvec.hpp>

#include <memory>
#include <vector>

namespace DRT
{
  class PackBuffer;
}

namespace MIXTURE
{
  template <typename T>
  class RemodelFiberMaterial;
  template <typename T>
  class LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution;


  namespace IMPLEMENTATION
  {
    template <int numstates, typename T>
    class RemodelFiberImplementation
    {
      struct GRState
      {
        T growth_scalar = 1.0;
        T lambda_r = 1.0;
        T lambda_f = 1.0;
        T lambda_ext = 1.0;
      };

     public:
      RemodelFiberImplementation(std::shared_ptr<const RemodelFiberMaterial<T>> material,
          LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<T> growth_evolution, T lambda_pre);

      /*!
       * @brief Pack all internal data into tha #data
       *
       * @param data (out) : buffer to serialize data to.
       */
      void Pack(DRT::PackBuffer& data) const;

      /*!
       * @brief Unpack all internal data that was previously packed by #Pack(DRT::PackBuffer&)
       *
       * @param position (in/out) : Position, where to start reading
       * @param data (in) : Vector of chars to extract data from
       */
      void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

      /// @brief Updates previous history data
      void Update();

      /*!
       * @brief Sets the deposition (homeostatic) stretch.
       *
       * @param lambda_pre
       */
      void UpdateDepositionStretch(T lambda_pre);

      /*!
       * @brief Set deformation state of the fiber
       *
       * @note This method has to be called before any Evaluation or local integration
       *
       * @param lambda_f (in) : total stretch in fiber direction
       * @param lambda_ext (in) : external inelastic stretch in fiber direction
       */
      void SetState(T lambda_f, T lambda_ext);

      [[nodiscard]] T EvaluateGrowthEvolutionEquationDt(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T EvaluateDGrowthEvolutionEquationDtDSig(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T EvaluateDGrowthEvolutionEquationDtPartialDgrowth(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T EvaluateDGrowthEvolutionEquationDtPartialDRemodel(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T EvaluateDGrowthEvolutionEquationDtDGrowth(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T EvaluateDGrowthEvolutionEquationDtDRemodel(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;

      [[nodiscard]] T EvaluateRemodelEvolutionEquationDt(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtDSig(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtDI4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtPartialDGrowth(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtPartialDRemodel(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtDGrowth(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDRemodelEvolutionEquationDtDRemodel(
          T lambda_f, T lambda_r, T lambda_ext) const;

      [[nodiscard]] T EvaluateFiberCauchyStress(T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDFiberCauchyStressPartialDI4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDFiberCauchyStressPartialDI4DI4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T EvaluateDFiberCauchyStressDRemodel(
          T lambda_f, T lambda_r, T lambda_ext) const;

      /// @name Methods for doing explicit or implicit time integration
      /// @{
      /*!
       * @brief Integrate the local evolution equation with an implicit time integration scheme.
       *
       * @param dt (in) : timestep
       *
       * @return Derivative of the residuum of the time integration scheme w.r.t. growth scalar and
       * lambda_r
       */
      CORE::LINALG::Matrix<2, 2, T> IntegrateLocalEvolutionEquationsImplicit(T dt);

      /*!
       * @brief Integrate the local evolution equation with an explicit time integration scheme.
       *
       * @param dt (in) : timestep
       */
      void IntegrateLocalEvolutionEquationsExplicit(T dt);
      /// @}
      /// @brief Evaluation methods
      ///
      /// @note It is important to call #SetState(#T) first.
      ///
      /// @{
      [[nodiscard]] T EvaluateCurrentHomeostaticFiberCauchyStress() const;
      [[nodiscard]] T EvaluateCurrentFiberCauchyStress() const;
      [[nodiscard]] T EvaluateCurrentFiberPK2Stress() const;
      [[nodiscard]] T EvaluateDCurrentFiberPK2StressDLambdafsq() const;
      [[nodiscard]] T EvaluateDCurrentFiberPK2StressDLambdar() const;
      [[nodiscard]] T EvaluateDCurrentGrowthEvolutionImplicitTimeIntegrationResiduumDLambdafsq(
          T dt) const;
      [[nodiscard]] T EvaluateDCurrentRemodelEvolutionImplicitTimeIntegrationResiduumDLambdafsq(
          T dt) const;
      [[nodiscard]] T EvaluateCurrentGrowthScalar() const;
      [[nodiscard]] T EvaluateCurrentLambdar() const;
      /// @}

      /// array of G&R states (the last state in the array is the current state)
      std::array<GRState, numstates> states_;

      /// homeostatic quantities
      /// @{
      T sig_h_ = 0.0;
      T lambda_pre_ = 1.0;
      /// @}

      /// Strain energy function of the fiber
      const std::shared_ptr<const RemodelFiberMaterial<T>> fiber_material_;

      /// Growth evolution equation
      const LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<T> growth_evolution_;

#ifdef DEBUG
      bool state_is_set_ = false;
#endif
    };
  }  // namespace IMPLEMENTATION
}  // namespace MIXTURE

#endif  // MIXTURE_REMODELFIBER_INTERNAL_H