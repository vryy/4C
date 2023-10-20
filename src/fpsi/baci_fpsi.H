/*----------------------------------------------------------------------*/
/*! \file

\level 2

*/

#ifndef BACI_FPSI_H
#define BACI_FPSI_H

#include "baci_adapter_algorithmbase.H"

#include <Epetra_Comm.h>

namespace FPSI
{
  class FPSI_Base : public ADAPTER::AlgorithmBase
  {
   public:
    /// constructor of base class
    FPSI_Base(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams);

    /// setup
    virtual void SetupSystem() = 0;

    /// setup solver
    virtual void SetupSolver() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) = 0;

    /// read restart
    void ReadRestart(int restartstep) override = 0;

    /// redistribute FPSI interface if running on parallel
    void RedistributeInterface();
  };
}  // namespace FPSI

#endif
