/*----------------------------------------------------------------------*/
/*! \file
\brief Method to deal with one dimensional artery inlet bcs


\level 3

*----------------------------------------------------------------------*/

#ifndef BACI_ART_NET_ART_TERMINAL_BC_H
#define BACI_ART_NET_ART_TERMINAL_BC_H

#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>



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

#endif  // ART_NET_ART_TERMINAL_BC_H
