/*----------------------------------------------------------------------*/
/*! \file

\brief Interfacing FSI problems with NOX

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_MONOLITHICINTERFACE_HPP
#define FOUR_C_FSI_MONOLITHICINTERFACE_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// Interface of monolithic algorithms to NOX group
  class MonolithicInterface
  {
   public:
    virtual ~MonolithicInterface() = default;
    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void setup_rhs(Epetra_Vector& f, bool firstcall = false) = 0;

    /// setup composed system matrix from field solvers
    virtual void setup_system_matrix() = 0;

    //@}

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void scale_system(Epetra_Vector& b) = 0;

    /// undo infnorm scaling from scaled solution
    virtual void unscale_solution(Epetra_Vector& x, Epetra_Vector& b) = 0;

    //@}
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
