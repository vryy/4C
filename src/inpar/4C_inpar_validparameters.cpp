// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_validparameters.hpp"

#include "4C_inpar_ale.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_inpar_beampotential.hpp"
#include "4C_inpar_binningstrategy.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_inpar_browniandyn.hpp"
#include "4C_inpar_cardiac_monodomain.hpp"
#include "4C_inpar_cardiovascular0d.hpp"
#include "4C_inpar_constraint_framework.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_cut.hpp"
#include "4C_inpar_ehl.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_inpar_elemag.hpp"
#include "4C_inpar_fbi.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_inpar_fpsi.hpp"
#include "4C_inpar_fs3i.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_geometric_search.hpp"
#include "4C_inpar_immersed.hpp"
#include "4C_inpar_io.hpp"
#include "4C_inpar_IO_monitor_structure_dbc.hpp"
#include "4C_inpar_IO_runtime_output.hpp"
#include "4C_inpar_IO_runtime_output_fluid.hpp"
#include "4C_inpar_IO_runtime_output_structure_beams.hpp"
#include "4C_inpar_IO_runtime_vtk_output_structure.hpp"
#include "4C_inpar_IO_runtime_vtp_output_structure.hpp"
#include "4C_inpar_levelset.hpp"
#include "4C_inpar_lubrication.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_mpc_rve.hpp"
#include "4C_inpar_particle.hpp"
#include "4C_inpar_pasi.hpp"
#include "4C_inpar_plasticity.hpp"
#include "4C_inpar_poroelast.hpp"
#include "4C_inpar_porofluidmultiphase.hpp"
#include "4C_inpar_poromultiphase.hpp"
#include "4C_inpar_poromultiphase_scatra.hpp"
#include "4C_inpar_poroscatra.hpp"
#include "4C_inpar_problemtype.hpp"
#include "4C_inpar_rebalance.hpp"
#include "4C_inpar_s2i.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_searchtree.hpp"
#include "4C_inpar_solver.hpp"
#include "4C_inpar_solver_nonlin.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_inpar_ssti.hpp"
#include "4C_inpar_sti.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_inpar_thermo.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_inpar_volmortar.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_io_input_file_utils.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_any.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <iostream>
#include <string>


FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void print_valid_parameters()
{
  std::shared_ptr<const Teuchos::ParameterList> list = Input::valid_parameters();
  list->print(std::cout,
      Teuchos::ParameterList::PrintOptions().showDoc(true).showFlags(false).indent(4).showTypes(
          false));
}


/*----------------------------------------------------------------------*/
//! Print help message
/*----------------------------------------------------------------------*/
void print_help_message()
{
  std::cout << "NAME\n"
            << "\t"
            << "4C - simulate just about anything\n"
            << "\n"
            << "SYNOPSIS\n"
            << "\t"
            << "4C [-h | --help] [-p | --parameters] [-d | --datfile] [-ngroup=<x>] \\ \n"
               "\t\t[-glayout=a,b,c,...] [-nptype=<parallelism_type>] \\ \n"
            << "\t\t<dat_name> <output_name> [restart=<y>] [restartfrom=restart_file_name] \\ \n"
               "\t\t[ <dat_name0> <output_name0> [restart=<y>] [restartfrom=restart_file_name] ... "
               "] \\ \n"
               "\t\t[--interactive]\n"
            << "\n"
            << "DESCRIPTION\n"
            << "\tThe am besten simulation tool in the world.\n"
            << "\n"
            << "OPTIONS\n"
            << "\t--help or -h\n"
            << "\t\tPrint this message.\n"
            << "\n"
            << "\t--parameters or -p\n"
            << "\t\tDumps information about the parameters for consumption by additional tools.\n"
            << "\n"
            << "\t--datfile or -d\n"
            << "\t\tPrint example dat_file with all available parameters.\n"
            << "\n"
            << "\t-ngroup=<x>\n"
            << "\t\tSpecify the number of groups for nested parallelism. (default: 1)\n"
            << "\n"
            << "\t-glayout=<a>,<b>,<c>,...\n"
            << "\t\tSpecify the number of processors per group. \n"
               "\t\tArgument \"-ngroup\" is mandatory and must be preceding. \n"
               "\t\t(default: equal distribution)\n"
            << "\n"
            << "\t-nptype=<parallelism_type>\n"
            << "\t\tAvailable options: \"separateDatFiles\" and \"everyGroupReadDatFile\"; \n"
               "\t\tMust be set if \"-ngroup\" > 1.\n"
            << "\t\t\"diffgroupx\" can be used to compare results from separate but parallel 4C "
               "runs; \n"
               "\t\tx must be 0 and 1 for the respective run\n"
            << "\n"
            << "\t<dat_name>\n"
            << "\t\tName of the input file, including the suffix (Usually *.dat)\n"
            << "\n"
            << "\t<output_name>\n"
            << "\t\tPrefix of your output files.\n"
            << "\n"
            << "\trestart=<y>\n"
            << "\t\tRestart the simulation from step <y>. \n"
               "\t\tIt always refers to the previously defined <dat_name> and <output_name>. \n"
               "\t\t(default: 0 or from <dat_name>)\n"
               "\t\tIf y=last_possible, it will restart from the last restart step defined in the "
               "control file.\n"
            << "\n"
            << "\trestartfrom=<restart_file_name>\n"
            << "\t\tRestart the simulation from the files prefixed with <restart_file_name>. \n"
               "\t\t(default: <output_name>)\n"
            << "\n"
            << "\t--interactive\n"
            << "\t\t4C waits at the beginning for keyboard input. \n"
               "\t\tHelpful for parallel debugging when attaching to a single job. \n"
               "\t\tMust be specified at the end in the command line.\n"
            << "\n"
            << "BUGS\n"
            << "\t100% bug free since 1964.\n"
            << "\n"
            << "TIPS\n"
            << "\tCan be obtain from a friendly colleague.\n"
            << "\n"
            << "\tAlso, espresso may be donated to room MW1236.\n";

  return;
}



void print_default_dat_header()
{
  std::shared_ptr<const Teuchos::ParameterList> list = Input::valid_parameters();
  Core::IO::InputFileUtils::print_dat(std::cout, *list);
}

std::shared_ptr<const Teuchos::ParameterList> Input::valid_parameters()
{
  std::shared_ptr<Teuchos::ParameterList> list = std::make_shared<Teuchos::ParameterList>();

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& discret = list->sublist("DISCRETISATION", false, "");

  Core::Utils::int_parameter("NUMFLUIDDIS", 1, "Number of meshes in fluid field", &discret);
  Core::Utils::int_parameter("NUMSTRUCDIS", 1, "Number of meshes in structural field", &discret);
  Core::Utils::int_parameter("NUMALEDIS", 1, "Number of meshes in ale field", &discret);
  Core::Utils::int_parameter(
      "NUMARTNETDIS", 1, "Number of meshes in arterial network field", &discret);
  Core::Utils::int_parameter("NUMTHERMDIS", 1, "Number of meshes in thermal field", &discret);
  Core::Utils::int_parameter("NUMAIRWAYSDIS", 1,
      "Number of meshes in reduced dimensional airways network field", &discret);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& size = list->sublist("PROBLEM SIZE", false, "");

  Core::Utils::int_parameter("DIM", 3, "2d or 3d problem", &size);

  // deactivate all the follwing (unused) parameters one day
  // they are nice as general info in the input file but should not
  // read into a parameter list. Misuse is possible
  Core::Utils::int_parameter("ELEMENTS", 0, "Total number of elements", &size);
  Core::Utils::int_parameter("NODES", 0, "Total number of nodes", &size);
  Core::Utils::int_parameter("NPATCHES", 0, "number of nurbs patches", &size);
  Core::Utils::int_parameter("MATERIALS", 0, "number of materials", &size);
  Core::Utils::int_parameter("NUMDF", 3, "maximum number of degrees of freedom", &size);

  Inpar::PROBLEMTYPE::set_valid_parameters(*list);

  /*----------------------------------------------------------------------*/

  Teuchos::ParameterList& nurbs_param = list->sublist(
      "NURBS", false, "Section to define information related to NURBS discretizations.");

  Core::Utils::bool_parameter("DO_LS_DBC_PROJECTION", "No",
      "Determines if a projection is needed for least square Dirichlet boundary conditions.",
      &nurbs_param);

  Core::Utils::int_parameter("SOLVER_LS_DBC_PROJECTION", -1,
      "Number of linear solver for the projection of least squares Dirichlet boundary conditions "
      "for NURBS "
      "discretizations",
      &nurbs_param);

  /*----------------------------------------------------------------------*/
  /* Finally call the problem-specific SetValidParameter functions        */
  /*----------------------------------------------------------------------*/

  Inpar::Solid::set_valid_parameters(*list);
  Inpar::IO::set_valid_parameters(*list);
  Inpar::IOMonitorStructureDBC::set_valid_parameters(*list);
  Inpar::IORuntimeOutput::set_valid_parameters(*list);
  Inpar::IORuntimeVTPStructure::set_valid_parameters(*list);
  Inpar::Mortar::set_valid_parameters(*list);
  Inpar::CONTACT::set_valid_parameters(*list);
  Inpar::VolMortar::set_valid_parameters(*list);
  Inpar::Wear::set_valid_parameters(*list);
  Inpar::IORuntimeOutput::FLUID::set_valid_parameters(*list);
  Inpar::IORuntimeOutput::Solid::set_valid_parameters(*list);
  Inpar::IORuntimeOutput::Beam::set_valid_parameters(*list);
  Inpar::BeamContact::set_valid_parameters(*list);
  Inpar::BEAMPOTENTIAL::set_valid_parameters(*list);
  Inpar::BEAMINTERACTION::set_valid_parameters(*list);
  Inpar::RveMpc::set_valid_parameters(*list);
  Inpar::BrownianDynamics::set_valid_parameters(*list);

  Inpar::Plasticity::set_valid_parameters(*list);

  Inpar::Thermo::set_valid_parameters(*list);
  Inpar::TSI::set_valid_parameters(*list);

  Inpar::FLUID::set_valid_parameters(*list);
  Inpar::LowMach::set_valid_parameters(*list);
  Inpar::Cut::set_valid_parameters(*list);
  Inpar::XFEM::set_valid_parameters(*list);
  Inpar::CONSTRAINTS::set_valid_parameters(*list);

  Inpar::Lubrication::set_valid_parameters(*list);
  Inpar::ScaTra::set_valid_parameters(*list);
  Inpar::LevelSet::set_valid_parameters(*list);
  Inpar::ElCh::set_valid_parameters(*list);
  Inpar::ElectroPhysiology::set_valid_parameters(*list);
  Inpar::STI::set_valid_parameters(*list);

  Inpar::S2I::set_valid_parameters(*list);
  Inpar::FS3I::set_valid_parameters(*list);
  Inpar::PoroElast::set_valid_parameters(*list);
  Inpar::PoroScaTra::set_valid_parameters(*list);
  Inpar::POROMULTIPHASE::set_valid_parameters(*list);
  Inpar::PoroMultiPhaseScaTra::set_valid_parameters(*list);
  Inpar::POROFLUIDMULTIPHASE::set_valid_parameters(*list);
  Inpar::EHL::set_valid_parameters(*list);
  Inpar::SSI::set_valid_parameters(*list);
  Inpar::SSTI::set_valid_parameters(*list);
  Inpar::ALE::set_valid_parameters(*list);
  Inpar::FSI::set_valid_parameters(*list);

  Inpar::ArtDyn::set_valid_parameters(*list);
  Inpar::ArteryNetwork::set_valid_parameters(*list);
  Inpar::BioFilm::set_valid_parameters(*list);
  Inpar::ReducedLung::set_valid_parameters(*list);
  Inpar::Cardiovascular0D::set_valid_parameters(*list);
  Inpar::Immersed::set_valid_parameters(*list);
  Inpar::FPSI::set_valid_parameters(*list);
  Inpar::FBI::set_valid_parameters(*list);

  Inpar::PARTICLE::set_valid_parameters(*list);

  Inpar::EleMag::set_valid_parameters(*list);

  Inpar::Geo::set_valid_parameters(*list);
  Inpar::BINSTRATEGY::set_valid_parameters(*list);
  Inpar::GeometricSearch::set_valid_parameters(*list);
  Inpar::PaSI::set_valid_parameters(*list);

  Inpar::Rebalance::set_valid_parameters(*list);
  Inpar::SOLVER::set_valid_parameters(*list);
  Inpar::NlnSol::set_valid_parameters(*list);

  return list;
}



FOUR_C_NAMESPACE_CLOSE
