/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_TIMEADA_ZIENXIE_HPP
#define FOUR_C_ADAPTER_STR_TIMEADA_ZIENXIE_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_timeada.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Solid::TimeInt
{
  class Base;
}  // namespace Solid::TimeInt

namespace Adapter
{
  /*! \brief Time step size adaptivity with Zienkiewicz-Xie error indicator
   *        only for lower or equal than second order accurate marching
   *        scheme
   *
   * References
   * - [1] OC Zienkiewicz and YM Xie, A simple error estimator and adaptive
   *   time stepping procedure for dynamic analysis, Earthquake Engrg.
   *   and Structural Dynamics, 20:871-887, 1991.
   */
  class StructureTimeAdaZienXie : public StructureTimeAda
  {
   public:
    /// constructor
    explicit StructureTimeAdaZienXie(Teuchos::RCP<Structure> structure)
        : StructureTimeAda(structure)
    {
    }

    //! Provide the name
    enum Inpar::Solid::TimAdaKind MethodName() const override
    {
      return Inpar::Solid::timada_kind_zienxie;
    }

    std::string MethodTitle() const override { return "ZienkiewiczXie"; }

    //! Provide local order of accuracy
    int method_order_of_accuracy_dis() const override { return 3; }

    //! Provide local order of accuracy
    int method_order_of_accuracy_vel() const override { return 2; }

    //! Return linear error coefficient of displacements
    double method_lin_err_coeff_dis() const override { return -1.0 / 24.0; }

    //! Return linear error coefficient of velocities
    double method_lin_err_coeff_vel() const override { return -1.0 / 12.0; }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const override { return ada_upward; }

   protected:
    /// setup of the auxiliary time integrator
    void setup_auxiliar() override {}

   private:
    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    void integrate_step_auxiliar() override;

    /*! \brief Update the auxiliar integrator
     */
    void update_auxiliar() override;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
