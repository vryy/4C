/*----------------------------------------------------------------------*/
/*!
\file strtimint_create.H
\brief Creation of structural time integrators in accordance with user's wishes

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

#ifndef STRTIMINT_CREATE_H
#define STRTIMINT_CREATE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Epetra_MpiComm.h"
#include "Epetra_Time.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class Solver;
}

namespace IO
{
  class DiscretizationWriter;
}

/*----------------------------------------------------------------------*/
/* Structural Dynamics */
namespace STR
{
  // forward declarations
  class TimInt;
  class TimIntImpl;
  class TimIntExpl;

  /*====================================================================*/
  //! Create marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimInt> TimIntCreate(
      const Teuchos::ParameterList& timeparams,       //!< time parameters
      const Teuchos::ParameterList& ioflags,          //!< input-output-flags
      const Teuchos::ParameterList& sdyn,             //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,          //!< extra flags
      Teuchos::RCP<DRT::Discretization>& actdis,      //!< discretisation
      Teuchos::RCP<LINALG::Solver>& solver,           //!< the solver
      Teuchos::RCP<LINALG::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b implicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntImpl> TimIntImplCreate(
      const Teuchos::ParameterList& timeparams,       //!< time parameters
      const Teuchos::ParameterList& ioflags,          //!< input-output-flags
      const Teuchos::ParameterList& sdyn,             //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,          //!< extra flags
      Teuchos::RCP<DRT::Discretization>& actdis,      //!< discretisation
      Teuchos::RCP<LINALG::Solver>& solver,           //!< the solver
      Teuchos::RCP<LINALG::Solver>& contactsolver,    //!< the contact solver
      Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

  /*====================================================================*/
  //! Create \b explicit marching time integrator convenience routine
  //!
  //! \author bborn \date 07/08
  Teuchos::RCP<STR::TimIntExpl> TimIntExplCreate(
      const Teuchos::ParameterList& timeparams,       //!< time parameters
      const Teuchos::ParameterList& ioflags,          //!< input-output-flags
      const Teuchos::ParameterList& sdyn,             //!< structural dynamic flags
      const Teuchos::ParameterList& xparams,          //!< extra flags
      Teuchos::RCP<DRT::Discretization>& actdis,      //!< discretisation
      Teuchos::RCP<LINALG::Solver>& solver,           //!< the solver
      Teuchos::RCP<LINALG::Solver>& contactsolver,    //!< the solver for contact/meshtying
      Teuchos::RCP<IO::DiscretizationWriter>& output  //!< output writer
  );

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif /* #ifndef STRTIMINT_CREATE_H */
