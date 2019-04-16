/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithicinterface.H

\brief Interfacing FSI problems with NOX

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHICINTERFACE_H
#define FSI_MONOLITHICINTERFACE_H

namespace FSI
{
  /// Interface of monolithic algorithms to NOX group
  class MonolithicInterface
  {
   public:
    virtual ~MonolithicInterface() {}

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall = false) = 0;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix() = 0;

    //@}

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(Epetra_Vector& b) = 0;

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(Epetra_Vector& x, Epetra_Vector& b) = 0;

    //@}
  };
}  // namespace FSI

#endif
