/*----------------------------------------------------------------------*/
/*!
\file fpsi.H

\level 2

\maintainer  Christoph Ager
*/

#ifndef FPSI_H
#define FPSI_H

#include <Epetra_Comm.h>
#include "../drt_adapter/adapter_algorithmbase.H"

namespace FPSI
{
  class FPSI_Base : public ADAPTER::AlgorithmBase
  {
   public:
    /// constructor of base class
    FPSI_Base(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams);

    /// destructor of base class
    virtual ~FPSI_Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// setup solver
    virtual void SetupSolver() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) = 0;

    /// read restart
    virtual void ReadRestart(int restartstep) = 0;

    /// redistribute FPSI interface if running on parallel
    void RedistributeInterface();
  };
}  // namespace FPSI

#endif
