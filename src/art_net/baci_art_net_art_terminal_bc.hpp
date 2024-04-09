/*----------------------------------------------------------------------*/
/*! \file
\brief Method to deal with one dimensional artery inlet bcs


\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_ART_TERMINAL_BC_HPP
#define FOUR_C_ART_NET_ART_TERMINAL_BC_HPP

#include "baci_config.hpp"

#include "baci_io.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN



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

BACI_NAMESPACE_CLOSE

#endif
