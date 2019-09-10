/*----------------------------------------------------------------------*/
/*! \file
\brief Method to deal with one dimensional artery inlet bcs

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/

#ifndef ART_TERMINAL_BC_H_
#define ART_TERMINAL_BC_H_


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include <Epetra_SerialDenseSolver.h>



namespace ART
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for inlet bc
    //--------------------------------------------------------------------

    /*!
    \brief 1d-artery inlet boundary condition, this class is meant to do
     solve the bc at the inlet of a one-dimensional arterial network
    \author ismail (ismail@lnm.mw.tum.de)
    */

    /*
    \brief Standard solver for 1d artery inlet
    */
    void SolvePrescribedTerminalBC(Teuchos::RCP<DRT::Discretization> actdis,
        const DRT::Condition* condition, Teuchos::ParameterList& params);

    /*
    \brief Standard solver for 1d artery reflective outlet
    */
    void SolveReflectiveTerminal(Teuchos::RCP<DRT::Discretization> actdis,
        const DRT::Condition* condition, Teuchos::ParameterList& params);

    /*
    \brief Standard solver for 1d artery explicit windkessel BC outlet
    */
    void SolveExplWindkesselBC(Teuchos::RCP<DRT::Discretization> actdis,
        const DRT::Condition* condition, Teuchos::ParameterList& params);

  }  // namespace UTILS
}  // namespace ART

#endif /* ART_TERMINAL_H_ */
