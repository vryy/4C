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
            << "\t\tPrint a list of all available parameters for use in a dat_file.\n"
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


void Input::print_documentation(std::ostream& stream, const Teuchos::ParameterEntry& entry)
{
  // Helper function to print documentation
  std::string doc = entry.docString();
  if (!doc.empty())
  {
    Teuchos::StrUtils::printLines(stream, "// ", doc);
  }
}


void Input::print_sublist(std::ostream& stream, const std::string& parentname,
    const std::string& name, const Teuchos::ParameterList& list, bool comment)
{
  // Helper function to print a sublist
  std::string secname = parentname;
  if (!secname.empty()) secname += "/";
  secname += name;
  unsigned l = secname.length();
  stream << "--" << std::string(std::max<int>(65 - l, 0), '-');
  stream << secname << "\n";
  print_dat_header(stream, list.sublist(name), secname, comment);
}

void Input::print_parameter(std::ostream& stream, const Teuchos::ParameterEntry& entry,
    const std::string& name, const Teuchos::ParameterList& list, bool comment)
{
  // Retrieve the parameter entry's validator (if any)
  Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

  // Print comments if requested
  if (comment)
  {
    // Check if the validator has valid string values
    if (validator != Teuchos::null)
    {
      Teuchos::RCP<const Teuchos::Array<std::string>> validValues = validator->validStringValues();

      // If valid values exist, print them
      if (validValues != Teuchos::null)
      {
        unsigned totalLength = 0;
        // Calculate the total length of all valid values
        for (const auto& value : *validValues)
        {
          totalLength += value.length() + 1;  // Include space/comma
        }
        // Print valid values in a compact or expanded format based on total length
        if (totalLength < 74)
        {
          // Print all values in a single line, separated by commas
          stream << "//     ";
          for (auto it = validValues->begin(); it != validValues->end(); ++it)
          {
            stream << *it;
            if (std::next(it) != validValues->end())
            {
              stream << ",";  // Add a comma if it's not the last element
            }
          }
          stream << "\n";
        }
        else
        {
          // Print each value on a new line
          for (const auto& value : *validValues)
          {
            stream << "//     " << value << '\n';
          }
        }
      }
    }
  }

  // Print the parameter's name and value
  const Teuchos::any& value = entry.getAny(false);
  stream << name;
  unsigned nameLength = name.length();
  // Ensure proper spacing for alignment
  stream << std::string(std::max<int>(31 - nameLength, 0), ' ');

  // Optionally print an equal sign if needed
  if (need_to_print_equal_sign(list)) stream << " =";

  try
  {
    // print true/false for bool values to distinguish them from type int
    if (value.type() == typeid(bool))
    {
      stream << " " << (Teuchos::any_cast<bool>(value) ? "true" : "false") << "\n";
    }
    else
    {
      // For non-boolean types, print the value directly
      stream << " " << value << "\n";
    }
  }
  catch (const Teuchos::NonprintableTypeException&)
  {
    // Handle non-printable enum class types
    stream << value.typeName() << "\n";
  }
}

void Input::print_dat_header(
    std::ostream& stream, const Teuchos::ParameterList& list, std::string parentname, bool comment)
{
  // Main loop over the parameter list that calls the helper functions to print
  // documentation, sublists or parameters:
  //
  // Iterate through the parameter list in two distinct phases to ensure proper ordering and
  // handling:
  // - **Phase 0**:
  //    Print all parameters that are not sublists. This ensures that top-level parameters
  //    are written to stream first, without any nested content interfering.
  // - **Phase 1**:
  //    Recursively handle and print all sublists. This phase is executed after all non-sublists
  //    have been processed, allowing sublists to be printed in their hierarchical order.
  //
  // By separating the iteration into these phases, we avoid issues with alphabetical
  // ordering that could cause invalid output sequences for nested lists.
  for (int iterationPhase = 0; iterationPhase < 2; ++iterationPhase)
  {
    for (auto paramIter = list.begin(); paramIter != list.end(); ++paramIter)
    {
      const Teuchos::ParameterEntry& entry = list.entry(paramIter);
      const std::string& name = list.name(paramIter);

      if ((entry.isList() && iterationPhase == 0) || (!entry.isList() && iterationPhase == 1))
      {
        continue;
      }
      if (comment)
      {
        stream << "//\n";
        print_documentation(stream, entry);
      }
      if (entry.isList())
      {
        print_sublist(stream, parentname, name, list, comment);
      }
      else
      {
        print_parameter(stream, entry, name, list, comment);
      }
    }
  }
  stream << std::endl;
}

void print_default_dat_header()
{
  std::shared_ptr<const Teuchos::ParameterList> list = Input::valid_parameters();
  Input::print_dat_header(std::cout, *list);
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
  Inpar::IORuntimeOutput::BEAMS::set_valid_parameters(*list);
  Inpar::BEAMCONTACT::set_valid_parameters(*list);
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


bool Input::need_to_print_equal_sign(const Teuchos::ParameterList& list)
{
  // Helper function to check if string contains a space.
  const auto string_has_space = [](const std::string& s)
  { return std::any_of(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); }); };

  return std::any_of(list.begin(), list.end(),
      [&](const auto& it)
      {
        // skip entries that are lists: they are allowed to have spaces
        if (it.second.isList()) return false;

        const std::string& name = it.key;

        const Teuchos::RCP<const Teuchos::Array<std::string>>& values_ptr =
            it.second.validator()->validStringValues();

        const bool value_has_space =
            (values_ptr != Teuchos::null) &&
            std::any_of(values_ptr->begin(), values_ptr->end(), string_has_space);

        return value_has_space || string_has_space(name);
      });
}
FOUR_C_NAMESPACE_CLOSE
