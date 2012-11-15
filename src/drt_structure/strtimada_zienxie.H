/*======================================================================*/
/*!
\file strtimada_zienxie.H

\brief ZienkiewiczXie time step indicator for time adaptivity

<pre>
Maintainer: Alexander Popp
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

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

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

    //! Provide local discretisation error
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
