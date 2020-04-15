/*----------------------------------------------------------------------*/
/*! \file

\brief General mass production pre interfave

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_GROWTH_EVOLUTION_H
#define BACI_GROWTH_EVOLUTION_H
namespace MIXTURE
{
  /*!
   * @brief Defines the evolution equation of net mass production during growth and remodeling
   *
   * The discretization in time has been done with a simple FD method (rho_dot =
   * (rho_{n+1}-rho_{n})/dt
   */
  class GrowthEvolution
  {
   public:
    /*!
     * \brief Default virtual destructor to ensure correct delete operations
     */
    virtual ~GrowthEvolution() = default;

    /*!
     * @brief Evaluates the right-hand side of the non-dimensionalized Growth ODE
     *
     * The ode is of the form
     * \f[
     *  \frac{\dot{\rho_0}}{\rho_0} = k ( \Delta \sigma )
     * \f]
     *
     * @param delta_sig Difference of the Cauchy fiber normal stress to the homeostatic stress
     *
     * @return The right hand side of the ode k ( \Delta \sigma )
     */
    virtual double EvaluateRHS(double delta_sig) const = 0;
  };
}  // namespace MIXTURE
#endif  // BACI_GROWTH_EVOLUTION_H
