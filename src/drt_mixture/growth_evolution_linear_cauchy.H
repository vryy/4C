/*----------------------------------------------------------------------*/
/*! \file

\brief 1D constitutive equation for growth (Linear in Cauchy fiber stress)

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_GROWTH_EVOLUTION_LINEAR_CAUCHY_H
#define BACI_GROWTH_EVOLUTION_LINEAR_CAUCHY_H
#include "growth_evolution.H"

namespace MIXTURE
{
  /*!
   * @brief Defines the evolution equation of net 1D mass production during growth and remodeling
   *
   * The discretization in time has been done with a simple FD method rho_dot =
   * (rho_{n+1}-rho_{n})/dt
   */
  class CauchyLinearGrowthEvolution : public GrowthEvolution
  {
   public:
    /*!
     * \brief Default constructor
     *
     * \param k_sig Gain parameter of linear growth equation
     */
    explicit CauchyLinearGrowthEvolution(double k_sig) : k_sig_(k_sig) {}

    /*!
     * \brief Evaluates the right hand side of the mass production differential equation of the form
     * rho_dot = f(delta_sig)
     *
     * \param delta_sig Difference in the fiber Cauchy stress
     * \return double Right hand side of the mass production differential equation
     */
    double EvaluateRHS(double delta_sig) const override { return k_sig_ * delta_sig; }

    /// Gain-Parameter controlling the speed of the growth process
    double k_sig_;
  };
}  // namespace MIXTURE

#endif  // BACI_GROWTH_EVOLUTION_LINEAR_CAUCHY_H
