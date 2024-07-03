/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_TIMEADA_JOINT_HPP
#define FOUR_C_ADAPTER_STR_TIMEADA_JOINT_HPP

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
  /*! \brief Adaptive time loop for structural simulations
   *
   *  This adaptive scheme uses another time integrator to estimate the local error
   */
  class StructureTimeAdaJoint : public StructureTimeAda
  {
   public:
    /// constructor
    explicit StructureTimeAdaJoint(Teuchos::RCP<Structure> structure);

    //! Provide the name
    enum Inpar::Solid::TimAdaKind MethodName() const override
    {
      return Inpar::Solid::timada_kind_joint_explicit;
    }

    std::string MethodTitle() const override;

    //! Provide local order of accuracy
    int method_order_of_accuracy_dis() const override;

    //! Provide local order of accuracy
    int method_order_of_accuracy_vel() const override;

    //! Return linear error coefficient of displacements
    double method_lin_err_coeff_dis() const override;

    //! Return linear error coefficient of velocities
    double method_lin_err_coeff_vel() const override;

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const override;

   protected:
    /// setup of the auxiliary time integrator
    void setup_auxiliar() override;

   private:
    //! type of adaptivity algorithm
    enum AdaEnum ada_;

    //! the auxiliary integrator
    Teuchos::RCP<Solid::TimeInt::Base> sta_;

    //! wrapper of the auxiliary integrator
    Teuchos::RCP<Structure> sta_wrapper_;

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

    /*! \brief Prepare repetition of current time step
     *
     *  Print to screen and reset certain quantities in case that the current time
     *  step has to be repeated.
     *
     *  \author mayr.mt \date 12/2013
     */
    void reset_step() override;
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
