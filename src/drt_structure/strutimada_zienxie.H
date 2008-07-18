/*======================================================================*/
/*!
\file strutimada_zienxie.H

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
#ifndef STRUTIMADA_ZIENXIE_H
#define STRUTIMADA_ZIENXIE_H

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

#include "strutimada.H"
#include "strutimint_genalpha.H"

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
   * [1] OC Zienkiewicz and YM Xie, A simple error estimator and adaptive
   *     time stepping procedure for dynamic analysis, Earthquake Engrg.
   *     Structural Dynamics, 20:871-887, 1991.
   *
   * \author bborn
   * \date 10/07
   */
  class StruTimAdaZienXie : public StruTimAda
  {

  public:

    //! Constructor
    StruTimAdaZienXie
    (
      const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
      const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
      Teuchos::RCP<StruTimInt> tis  //!< marching time integrator
    );

    //! Destructor
    ~StruTimAdaZienXie() { ; }

    //! @name Actions
    //@{

    //! Provide local discretisation error
    void EvaluateLocalErrorDis();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    enum StruTimAda::NameEnum MethodName()
    {
      return StruTimAda::name_zienxie;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracy()
    {
      return 3;
    }

    //@}

  protected:

    //! not wanted: = operator
    StruTimAdaZienXie operator = (const StruTimAdaZienXie& old);
  
    //! not wanted: copy constructor
    StruTimAdaZienXie(const StruTimAdaZienXie& old);

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMADA_ZIENXIE_H
#endif  // #ifdef CCADISCRET
