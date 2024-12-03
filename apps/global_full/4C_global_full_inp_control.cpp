// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_global_data_read.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.hpp"

#include <utility>

void setup_parallel_output(std::string& outputfile_kenner, MPI_Comm lcomm, int group);

/*----------------------------------------------------------------------*
  | general input of the problem to be solved              m.gee 10/06  |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret(
    std::string& inputfile_name, std::string& outputfile_kenner, std::string& restartfile_kenner)
{
  using namespace FourC;

  Global::Problem* problem = Global::Problem::instance();
  MPI_Comm lcomm = problem->get_communicators()->local_comm();
  int group = problem->get_communicators()->group_id();
  Core::Communication::NestedParallelismType npType = problem->get_communicators()->np_type();

  // and now the actual reading
  Core::IO::InputFile reader(inputfile_name, lcomm);

  Global::read_parameter(*problem, reader);

  setup_parallel_output(outputfile_kenner, lcomm, group);

  // create control file for output and read restart data if required
  problem->open_control_file(lcomm, inputfile_name, outputfile_kenner, restartfile_kenner);

  // input of materials
  Global::read_materials(*problem, reader);

  // input of materials
  Global::read_contact_constitutive_laws(*problem, reader);

  // input of materials of cloned fields (if needed)
  Global::read_cloning_material_map(*problem, reader);

  {
    Core::Utils::FunctionManager function_manager;
    global_legacy_module_callbacks().AttachFunctionDefinitions(function_manager);
    function_manager.read_input(reader);
    problem->set_function_manager(std::move(function_manager));
  }

  // input of particles
  Global::read_particles(*problem, reader);

  switch (npType)
  {
    case Core::Communication::NestedParallelismType::no_nested_parallelism:
    case Core::Communication::NestedParallelismType::every_group_read_dat_file:
    case Core::Communication::NestedParallelismType::separate_dat_files:
      // input of fields
      Global::read_fields(*problem, reader);

      // read result tests
      Global::read_result(*problem, reader);

      // read all types of geometry related conditions (e.g. boundary conditions)
      // Also read time and space functions and local coord systems
      Global::read_conditions(*problem, reader);

      // read all knot information for isogeometric analysis
      // and add it to the (derived) nurbs discretization
      Global::read_knots(*problem, reader);
      break;
    default:
      FOUR_C_THROW("nptype (nested parallelity type) not recognized");
      break;
  }

  // all reading is done at this point!
  if (Core::Communication::my_mpi_rank(lcomm) == 0) problem->write_input_parameters();

  // before we destroy the reader we want to know about unused sections
  const bool all_ok = !reader.print_unknown_sections(std::cout);

  // we wait till all procs are here. Otherwise a hang up might occur where
  // one proc ended with FOUR_C_THROW but other procs were not finished and waited...
  // we also want to have the printing above being finished.
  Core::Communication::barrier(lcomm);
  FOUR_C_ASSERT_ALWAYS(all_ok,
      "Unknown sections detected. Correct this! Find hints on these unknown sections above.");
}


/*----------------------------------------------------------------------*
  | setup parallel output                                  ghamm 11/12  |
 *----------------------------------------------------------------------*/
void setup_parallel_output(std::string& outputfile_kenner, MPI_Comm lcomm, int group)
{
  using namespace FourC;

  // configure the parallel output environment
  const Teuchos::ParameterList& io = Global::Problem::instance()->io_params();
  bool screen = io.get<bool>("WRITE_TO_SCREEN");
  bool file = io.get<bool>("WRITE_TO_FILE");
  bool preGrpID = io.get<bool>("PREFIX_GROUP_ID");
  int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
  auto level = Teuchos::getIntegralValue<Core::IO::Verbositylevel>(io, "VERBOSITY");

  Core::IO::cout.setup(
      screen, file, preGrpID, level, std::move(lcomm), oproc, group, outputfile_kenner);
}
