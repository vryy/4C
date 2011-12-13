/*----------------------------------------------------------------------*/
/*!
\file strtimint_create.H
\brief Creation of structural time integrators in accordance with user's wishes

<pre>
Maintainer: Thomas Klöppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef STRTIMINT_CREATE_H
#define STRTIMINT_CREATE_H

/*----------------------------------------------------------------------*/
/* headers */
#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"


/*----------------------------------------------------------------------*/
/* Structural Dynamics */
namespace STR
{

  /*====================================================================*/
  //! Create marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimInt> TimIntCreate
  (
    const Teuchos::ParameterList& ioflags,  //!< input-output-flags
    const Teuchos::ParameterList& sdyn,  //!< structural dynamic flags
    const Teuchos::ParameterList& xparams,  //!< extra flags
    Teuchos::RCP<DRT::Discretization>& actdis,  //!< discretisation
    Teuchos::RCP<LINALG::Solver>& solver,  //!< the solver
    Teuchos::RCP<LINALG::Solver>& contactsolver,  //!< the solver for contact/meshtying
    Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b implicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntImpl> TimIntImplCreate
  (
    const Teuchos::ParameterList& ioflags,  //!< input-output-flags
    const Teuchos::ParameterList& sdyn,  //!< structural dynamic flags
    const Teuchos::ParameterList& xparams,  //!< extra flags
    Teuchos::RCP<DRT::Discretization>& actdis,  //!< discretisation
    Teuchos::RCP<LINALG::Solver>& solver,  //!< the solver
    Teuchos::RCP<LINALG::Solver>& contactsolver,  //!< the contact solver
    Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b explicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntExpl> TimIntExplCreate
  (
    const Teuchos::ParameterList& ioflags,  //!< input-output-flags
    const Teuchos::ParameterList& sdyn,  //!< structural dynamic flags
    const Teuchos::ParameterList& xparams,  //!< extra flags
    Teuchos::RCP<DRT::Discretization>& actdis,  //!< discretisation
    Teuchos::RCP<LINALG::Solver>& solver,  //!< the solver
    Teuchos::RCP<LINALG::Solver>& contactsolver,  //!< the solver for contact/meshtying
    Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  /* #ifndef STRTIMINT_CREATE_H */
#endif  /* #ifdef CCADISCRET */
