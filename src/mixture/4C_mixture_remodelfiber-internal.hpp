/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a 1D remodel fiber implementation
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MIXTURE_REMODELFIBER_INTERNAL_HPP
#define FOUR_C_MIXTURE_REMODELFIBER_INTERNAL_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Sacado_tradvec.hpp>

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::COMM
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
      void Pack(CORE::COMM::PackBuffer& data) const;

      /*!
       * @brief Unpack all internal data that was previously packed by
       * #Pack(CORE::COMM::PackBuffer&)
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
      void update_deposition_stretch(T lambda_pre);

      /*!
       * @brief Set deformation state of the fiber
       *
       * @note This method has to be called before any Evaluation or local integration
       *
       * @param lambda_f (in) : total stretch in fiber direction
       * @param lambda_ext (in) : external inelastic stretch in fiber direction
       */
      void set_state(T lambda_f, T lambda_ext);

      [[nodiscard]] T evaluate_growth_evolution_equation_dt(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T evaluate_d_growth_evolution_equation_dt_d_sig(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T evaluate_d_growth_evolution_equation_dt_partial_dgrowth(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T evaluate_d_growth_evolution_equation_dt_partial_d_remodel(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T evaluate_d_growth_evolution_equation_dt_d_growth(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;
      [[nodiscard]] T evaluate_d_growth_evolution_equation_dt_d_remodel(
          T lambda_f, T lambda_r, T lambda_ext, T growth_scalar) const;

      [[nodiscard]] T evaluate_remodel_evolution_equation_dt(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_d_sig(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_d_i4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_partial_d_growth(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_partial_d_remodel(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_d_growth(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_remodel_evolution_equation_dt_d_remodel(
          T lambda_f, T lambda_r, T lambda_ext) const;

      [[nodiscard]] T evaluate_fiber_cauchy_stress(T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_fiber_cauchy_stress_partial_d_i4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_fiber_cauchy_stress_partial_d_i4_d_i4(
          T lambda_f, T lambda_r, T lambda_ext) const;
      [[nodiscard]] T evaluate_d_fiber_cauchy_stress_d_remodel(
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
      CORE::LINALG::Matrix<2, 2, T> integrate_local_evolution_equations_implicit(T dt);

      /*!
       * @brief Integrate the local evolution equation with an explicit time integration scheme.
       *
       * @param dt (in) : timestep
       */
      void integrate_local_evolution_equations_explicit(T dt);
      /// @}
      /// @brief Evaluation methods
      ///
      /// @note It is important to call #set_state(#T) first.
      ///
      /// @{
      [[nodiscard]] T evaluate_current_homeostatic_fiber_cauchy_stress() const;
      [[nodiscard]] T evaluate_current_fiber_cauchy_stress() const;
      [[nodiscard]] T evaluate_current_fiber_p_k2_stress() const;
      [[nodiscard]] T evaluate_d_current_fiber_p_k2_stress_d_lambdafsq() const;
      [[nodiscard]] T evaluate_d_current_fiber_p_k2_stress_d_lambdar() const;
      [[nodiscard]] T
      evaluate_d_current_growth_evolution_implicit_time_integration_residuum_d_lambdafsq(
          T dt) const;
      [[nodiscard]] T
      evaluate_d_current_remodel_evolution_implicit_time_integration_residuum_d_lambdafsq(
          T dt) const;
      [[nodiscard]] T evaluate_current_growth_scalar() const;
      [[nodiscard]] T evaluate_current_lambdar() const;
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

#ifdef FOUR_C_ENABLE_ASSERTIONS
      bool state_is_set_ = false;
#endif
    };
  }  // namespace IMPLEMENTATION
}  // namespace MIXTURE

FOUR_C_NAMESPACE_CLOSE

#endif