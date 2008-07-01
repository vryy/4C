/*----------------------------------------------------------------------*/
/*!
\file strutimint.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_H
#define STRUTIMINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_constraint/uzawasolver.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_surfstress/drt_potential_manager.H"
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


/*======================================================================*/
/*!
\brief Front-end for structural dynamics by structure time integration. 
\author bborn
\date 06/08
*/
class StruTimInt
{

public:

  //! Structural result test is a befriended class
  friend class StruResultTest;

  //! Constructor
  StruTimInt()
  { 
    /*
    std::cout << "Welcome to Structural Time Integration _ " << std::endl;
    std::cout << " _ __       _       _      _ __    _  / |" << std::endl;
    std::cout << "| '_ \\     | |_____\\ \\    | '_ \\ _| |_| |" << std::endl;
    std::cout << "| | | |    |  _____   )   | | | |_   _| |" << std::endl;
    std::cout << "|_| |_|    |_|     /_/    |_| |_| |_| |_|" << std::endl;
    */
    std::cout << "Welcome to Structural Time Integration " << std::endl;
    std::cout << "     __o__                          __o__       " << std::endl;
    std::cout << "__  /-----\\__                  __  /-----\\__           " << std::endl;
    std::cout << "\\ \\/       \\ \\    |       \\    \\ \\/       \\ \\          " << std::endl;
    std::cout << " \\ |  tea  | |    |-------->    \\ |  tea  | |          " << std::endl;
    std::cout << "  \\|       |_/    |       /      \\|       |_/          " << std::endl;
    std::cout << "    \\_____/   ._                   \\_____/   ._ _|_ /| " << std::endl;
    std::cout << "              | |                            | | |   | " << std::endl;
    std::cout << std::endl;
  };

  //! Destructor
  virtual ~StruTimInt() { };

  //! Do time integration of multiple steps
  virtual void Integrate() = 0;

  //! Do time integration of single step
  virtual void IntegrateStep() = 0;

  //! Access Dicretisation
  virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

  //! Displacements
  virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;

  //! Velocities
  virtual Teuchos::RCP<Epetra_Vector> Vel() = 0;

  //!  Accelerations
  virtual Teuchos::RCP<Epetra_Vector> Acc() = 0;

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_H
#endif  // #ifdef CCADISCRET
