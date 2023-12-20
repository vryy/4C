/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_TIMEADA_JOINT_H
#define BACI_ADAPTER_STR_TIMEADA_JOINT_H

#include "baci_config.H"

#include "baci_adapter_str_timeada.H"

BACI_NAMESPACE_OPEN

// forward declaration
namespace STR::TIMINT
{
  class Base;
}  // namespace STR::TIMINT

namespace ADAPTER
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
    enum INPAR::STR::TimAdaKind MethodName() const override
    {
      return INPAR::STR::timada_kind_joint_explicit;
    }

    std::string MethodTitle() const override;

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyDis() const override;

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyVel() const override;

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const override;

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const override;

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const override;

   protected:
    /// setup of the auxiliary time integrator
    void SetupAuxiliar() override;

   private:
    //! type of adaptivity algorithm
    enum AdaEnum ada_;

    //! the auxiliary integrator
    Teuchos::RCP<STR::TIMINT::Base> sta_;

    //! wrapper of the auxiliary integrator
    Teuchos::RCP<Structure> sta_wrapper_;

    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    void IntegrateStepAuxiliar() override;

    /*! \brief Update the auxiliar integrator
     */
    void UpdateAuxiliar() override;

    /*! \brief Prepare repetition of current time step
     *
     *  Print to screen and reset certain quantities in case that the current time
     *  step has to be repeated.
     *
     *  \author mayr.mt \date 12/2013
     */
    void ResetStep() override;
  };

}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
