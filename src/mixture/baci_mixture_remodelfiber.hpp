/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of a 1D remodel fiber
\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_REMODELFIBER_HPP
#define BACI_MIXTURE_REMODELFIBER_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.hpp"

#include <memory>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::COMM
{
  class PackBuffer;
}

namespace MIXTURE
{
  template <typename T>
  class RemodelFiberMaterial;


  namespace IMPLEMENTATION
  {
    template <int numstates, typename T>
    class RemodelFiberImplementation;
  }

  template <int numstates>
  class RemodelFiber
  {
    struct GRState
    {
      double growth_scalar = 1.0;
      double lambda_r = 1.0;
      double lambda_f = 1.0;
    };

   public:
    RemodelFiber(std::shared_ptr<const RemodelFiberMaterial<double>> material,
        LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution,
        double lambda_pre);

    /*!
     * @brief Pack all internal data into tha #data
     *
     * @param data (out) : buffer to serialize data to.
     */
    void Pack(CORE::COMM::PackBuffer& data) const;

    /*!
     * @brief Unpack all internal data that was previously packed by #Pack(CORE::COMM::PackBuffer&)
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
    void UpdateDepositionStretch(double lambda_pre);

    /*!
     * @brief Set deformation state of the fiber
     *
     * @note This method has to be called before any Evaluation or local integration
     *
     * @param lambda_f (in) : total stretch in fiber direction
     * @param lambda_ext (in) : inelastic external stretch
     */
    void SetState(double lambda_f, double lambda_ext);

   public:
    /// @brief Evaluation methods
    ///
    /// @note It is important to call #SetState(double) first.
    ///
    /// @{

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
    CORE::LINALG::Matrix<2, 2> IntegrateLocalEvolutionEquationsImplicit(double dt);

    /*!
     * @brief Integrate the local evolution equation with an explicit time integration scheme.
     *
     * @param dt (in) : timestep
     */
    void IntegrateLocalEvolutionEquationsExplicit(double dt);
    /// @}
    [[nodiscard]] double EvaluateCurrentHomeostaticFiberCauchyStress() const;
    [[nodiscard]] double EvaluateCurrentFiberCauchyStress() const;
    [[nodiscard]] double EvaluateCurrentFiberPK2Stress() const;
    [[nodiscard]] double EvaluateDCurrentFiberPK2StressDLambdafsq() const;
    [[nodiscard]] double EvaluateDCurrentFiberPK2StressDLambdar() const;
    [[nodiscard]] double EvaluateDCurrentGrowthEvolutionImplicitTimeIntegrationResiduumDLambdafsq(
        double dt) const;
    [[nodiscard]] double EvaluateDCurrentRemodelEvolutionImplicitTimeIntegrationResiduumDLambdafsq(
        double dt) const;
    [[nodiscard]] double EvaluateCurrentGrowthScalar() const;
    [[nodiscard]] double EvaluateCurrentLambdar() const;
    /// @}

   private:
    const std::shared_ptr<IMPLEMENTATION::RemodelFiberImplementation<2, double>> impl_;
  };
}  // namespace MIXTURE

BACI_NAMESPACE_CLOSE

#endif  // MIXTURE_REMODELFIBER_H
