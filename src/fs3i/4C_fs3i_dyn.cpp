// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fs3i_dyn.hpp"

#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fs3i.hpp"
#include "4C_fs3i_biofilm_fsi.hpp"
#include "4C_fs3i_fps3i_partitioned_1wc.hpp"
#include "4C_fs3i_partitioned_1wc.hpp"
#include "4C_fs3i_partitioned_2wc.hpp"
#include "4C_global_data.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
// entry point for all kinds of FS3I
/*----------------------------------------------------------------------*/
void fs3i_dyn()
{
  MPI_Comm comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  std::shared_ptr<FS3I::FS3IBase> fs3i;

  // what's the current problem type?
  Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();

  switch (probtype)
  {
    case Core::ProblemType::gas_fsi:
    {
      fs3i = std::make_shared<FS3I::PartFS3I1Wc>(comm);
    }
    break;
    case Core::ProblemType::thermo_fsi:
    {
      fs3i = std::make_shared<FS3I::PartFS3I2Wc>(comm);
    }
    break;
    case Core::ProblemType::biofilm_fsi:
    {
      fs3i = std::make_shared<FS3I::BiofilmFSI>(comm);
    }
    break;
    case Core::ProblemType::fps3i:
    {
      fs3i = std::make_shared<FS3I::PartFpS3I1Wc>(comm);
    }
    break;
    default:
      FOUR_C_THROW("solution of unknown problemtype %d requested", probtype);
      break;
  }

  fs3i->init();
  fs3i->setup();

  // read the restart information, set vectors and variables ---
  // be careful, dofmaps might be changed here in a redistribute() call
  fs3i->read_restart();

  // if running FPS3I in parallel one needs to redistribute the interface after restarting
  fs3i->redistribute_interface();

  // now do the coupling and create combined dofmaps
  fs3i->setup_system();

  fs3i->timeloop();

  fs3i->test_results(comm);

  std::shared_ptr<const Teuchos::Comm<int>> TeuchosComm =
      Core::Communication::to_teuchos_comm<int>(comm);
  Teuchos::TimeMonitor::summarize(Teuchos::Ptr(TeuchosComm.get()), std::cout, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
