/*----------------------------------------------------------------------*/
/*! \file

\level 2

*/

#ifndef FOUR_C_FPSI_HPP
#define FOUR_C_FPSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

namespace FPSI
{
  class FpsiBase : public Adapter::AlgorithmBase
  {
   public:
    /// constructor of base class
    FpsiBase(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams);

    /// setup
    virtual void setup_system() = 0;

    /// setup solver
    virtual void setup_solver() = 0;

    /// timeloop of coupled problem
    virtual void timeloop() = 0;

    /// test results (if necessary)
    virtual void test_results(const Epetra_Comm& comm) = 0;

    /// read restart
    void read_restart(int restartstep) override = 0;

    /// redistribute FPSI interface if running on parallel
    void redistribute_interface();
  };
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
