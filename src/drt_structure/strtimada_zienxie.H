/*======================================================================*/
/*!
\file strtimada_zienxie.H

\brief ZienkiewiczXie time step indicator for time adaptivity

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET
#ifndef STRTIMADA_ZIENXIE_H
#define STRTIMADA_ZIENXIE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_contact/drt_contact_manager.H"

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
#include "strtimint_genalpha.H"

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
    enum TimAda::NameEnum MethodName()
    {
      return TimAda::name_zienxie;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyDis()
    {
      return 3;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracyVel()
    {
      return 2;
    }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis()
    {
      return -1.0/24.0;
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel()
    {
      return -1.0/12.0;
    }

    //! Provide type of algorithm
    enum AdaEnum MethodAdaptDis()
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
#endif  // #ifdef CCADISCRET
