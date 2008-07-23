/*======================================================================*/
/*!
\file strtimada_ab2.H

\brief Adams-Bashforth2 indicator for time adaptivity

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
#ifndef STRUTIMADA_AB2_H
#define STRUTIMADA_AB2_H

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
#include "strtimint_ab2.H"

/*----------------------------------------------------------------------*/
/* belongs to structure namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Time step size adaptivity with Adams-Bashforth2 error indicator
   *
   * \author bborn
   * \date 10/07
   */
  class StruTimAdaAB2 : public StruTimAda
  {

  public:

    //! @name Life
    //@{

    //! Full-fledged constructor
    StruTimAdaAB2
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
      Teuchos::RCP<StruTimInt>& tis  //!< marching time integrator
    );

    //! Slender constructor
    StruTimAdaAB2
    (
      const Teuchos::ParameterList& sdynparams,  //!< TIS input parameters
      const Teuchos::ParameterList& adaparams,  //!< adaptive input flags
      Teuchos::RCP<StruTimInt>& tis  //!< marching time integrator
    );

    //! Destructor
    ~StruTimAdaAB2() { ; }

    //@}

    //! @name Actions
    //@{

    //! Provide local discretisation error
    void IntegrateStepAuxiliar();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    enum StruTimAda::NameEnum MethodName()
    {
      return StruTimAda::name_ab2;
    }

    //! Provide local order of accuracy
    int MethodOrderOfAccuracy()
    {
      return 2;
    }

    //@}

  protected:

    //! not wanted: = operator
    StruTimAdaAB2 operator = (const StruTimAdaAB2& old);
  
    //! not wanted: copy constructor
    StruTimAdaAB2(const StruTimAdaAB2& old);

    //! The AB2 integrator 
    Teuchos::RCP<StruTimIntAB2> ab2_;

  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMADA_AB2_H
#endif  // #ifdef CCADISCRET
