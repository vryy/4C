/*======================================================================*/
/*!
\file strtimada_zienxie.H
\brief Zienkiewicz-Xie time step indicator for time adaptivity
\level 1

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef STRTIMADA_ZIENXIE_H
#define STRTIMADA_ZIENXIE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "strtimada.H"

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
    TimAdaZienXie
    (
      const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
      const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
      Teuchos::RCP<TimInt> tis  //!< marching time integrator
    );

    //! Destructor
    ~TimAdaZienXie() { ; }

    //! @name Actions
    //@{

    //! Finalize the class initialization (nothing to do here)
    void Init(Teuchos::RCP<TimInt>& sti) {}

    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    void IntegrateStepAuxiliar();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    enum INPAR::STR::TimAdaKind MethodName() const
    {
      return INPAR::STR::timada_kind_zienxie;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyDis() const
    {
      return 3;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyVel() const
    {
      return 2;
    }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis() const
    {
      return -1.0/24.0;
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel() const
    {
      return -1.0/12.0;
    }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis() const
    {
      return ada_upward;
    }

    //@}

  protected:

    //! not wanted: = operator
    TimAdaZienXie operator = (const TimAdaZienXie& old);

    //! not wanted: copy constructor
    TimAdaZienXie(const TimAdaZienXie& old);

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_ZIENXIE_H
