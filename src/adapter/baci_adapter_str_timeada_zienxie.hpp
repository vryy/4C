/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_TIMEADA_ZIENXIE_HPP
#define BACI_ADAPTER_STR_TIMEADA_ZIENXIE_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_timeada.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace STR::TIMINT
{
  class Base;
}  // namespace STR::TIMINT

namespace ADAPTER
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
    enum INPAR::STR::TimAdaKind MethodName() const override
    {
      return INPAR::STR::timada_kind_zienxie;
    }

    std::string MethodTitle() const override { return "ZienkiewiczXie"; }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyDis() const override { return 3; }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyVel() const override { return 2; }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const override { return -1.0 / 24.0; }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const override { return -1.0 / 12.0; }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const override { return ada_upward; }

   protected:
    /// setup of the auxiliary time integrator
    void SetupAuxiliar() override {}

   private:
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
  };

}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
