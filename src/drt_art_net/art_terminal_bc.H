/*!----------------------------------------------------------------------
\file art_terminal.H
\brief Method to deal with one dimensional artery inlet bcs 

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

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

using namespace std;
using namespace Teuchos;

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
    void SolvePrescribedTerminalBC(RCP<DRT::Discretization> actdis,
                                   const DRT::Condition *   condition,
                                   ParameterList &          params);

    /*
    \brief Standard solver for 1d artery reflective outlet
    */
    void SolveReflectiveTerminal(RCP<DRT::Discretization> actdis,
                                 const DRT::Condition *           condition,
                                 ParameterList &                  params);

    /*
    \brief Standard solver for 1d artery explicit windkessel BC outlet
    */
    void SolveExplWindkesselBC(RCP<DRT::Discretization> actdis,
                               const DRT::Condition *           condition,
                               ParameterList &                  params);

  } // namespace UTILS
} // namespace ART

#endif /* ART_TERMINAL_H_ */

