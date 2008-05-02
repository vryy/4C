#ifdef CCADISCRET

#ifndef FSI_DEBUGWRITER_H
#define FSI_DEBUGWRITER_H

#include <Teuchos_RCP.hpp>

#include "../drt_adapter/adapter_coupling.H"

namespace IO
{
  class OutputControl;
  class DiscretizationWriter;
}

namespace FSI
{

  /// Helper to write unconverged results at the FSI interface
  class DebugWriter
  {
  public:

    /// create FSI debug writer on given field discretization
    explicit DebugWriter(Teuchos::RCP<DRT::Discretization> dis);

    /// announce new time step
    /*!
      Results in a new control file
     */
    void NewTimeStep(int step);

    /// announce new interface iteration step
    void NewIteration();

    /// write interface vector
    /*!
      This includes an internal transfer of the interface vector to the
      interface discretization.
     */
    void WriteVector(const std::string& name, const Epetra_Vector& v);

  private:

    /// internal interface discretization
    Teuchos::RCP<DRT::Discretization> dis_;

    /// coupling of field discretization and interface discretization
    ADAPTER::Coupling coup_;

    /// current control file
    Teuchos::RCP<IO::OutputControl> control_;

    /// writer to control file
    Teuchos::RCP<IO::DiscretizationWriter> writer_;

    /// internal FSI iteration count
    int itnum_;
  };
}

#endif
#endif
