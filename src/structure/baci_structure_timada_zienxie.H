/*======================================================================*/
/*! \file
\brief Zienkiewicz-Xie time step indicator for time adaptivity
\level 1

*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef BACI_STRUCTURE_TIMADA_ZIENXIE_H
#define BACI_STRUCTURE_TIMADA_ZIENXIE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_structure_timada.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* belongs to structure namespace */
namespace STR
{
  /*====================================================================*/
  /*!
   * \brief Time step size adaptivity with Zienkiewicz-Xie error indicator
   *        only for lower or equal than second order accurate marching
   *        scheme
   *
   * References
   * - [1] OC Zienkiewicz and YM Xie, A simple error estimator and adaptive
   *   time stepping procedure for dynamic analysis, Earthquake Engrg.
   *   and Structural Dynamics, 20:871-887, 1991.
   *
   * \author bborn
   * \date 10/07
   */
  class TimAdaZienXie : public TimAda
  {
   public:
    //! Constructor
    TimAdaZienXie(const Teuchos::ParameterList& timeparams,  //!< TIS input parameters
        const Teuchos::ParameterList& adaparams,             //!< adaptive input flags
        Teuchos::RCP<TimInt> tis                             //!< marching time integrator
    );

    //! @name Actions
    //@{

    //! Finalize the class initialization (nothing to do here)
    void Init(Teuchos::RCP<TimInt>& sti) override {}

    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    void IntegrateStepAuxiliar() override;

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    enum INPAR::STR::TimAdaKind MethodName() const override
    {
      return INPAR::STR::timada_kind_zienxie;
    }

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

    //@}

   protected:
    //! not wanted: = operator
    TimAdaZienXie operator=(const TimAdaZienXie& old);

    //! not wanted: copy constructor
    TimAdaZienXie(const TimAdaZienXie& old);
  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // STRUCTURE_TIMADA_ZIENXIE_H
