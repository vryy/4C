// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FPSI_HPP
#define FOUR_C_FPSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"


FOUR_C_NAMESPACE_OPEN

namespace FPSI
{
  class FpsiBase : public Adapter::AlgorithmBase
  {
   public:
    /// constructor of base class
    FpsiBase(MPI_Comm comm, const Teuchos::ParameterList& fpsidynparams);

    /// setup
    virtual void setup_system() = 0;

    /// setup solver
    virtual void setup_solver() = 0;

    /// timeloop of coupled problem
    virtual void timeloop() = 0;

    /// test results (if necessary)
    virtual void test_results(MPI_Comm comm) = 0;

    /// read restart
    void read_restart(int restartstep) override = 0;

    /// redistribute FPSI interface if running on parallel
    void redistribute_interface();
  };
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
