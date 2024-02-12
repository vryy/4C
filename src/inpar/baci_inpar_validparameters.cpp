/*----------------------------------------------------------------------*/
/*! \file

\brief Setup of the list of valid input parameters

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_validparameters.hpp"

#include "baci_global_data_enums.hpp"
#include "baci_inpar.hpp"
#include "baci_inpar_ale.hpp"
#include "baci_inpar_beamcontact.hpp"
#include "baci_inpar_beaminteraction.hpp"
#include "baci_inpar_beampotential.hpp"
#include "baci_inpar_binningstrategy.hpp"
#include "baci_inpar_bio.hpp"
#include "baci_inpar_browniandyn.hpp"
#include "baci_inpar_cardiac_monodomain.hpp"
#include "baci_inpar_cardiovascular0d.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_inpar_cut.hpp"
#include "baci_inpar_ehl.hpp"
#include "baci_inpar_elch.hpp"
#include "baci_inpar_elemag.hpp"
#include "baci_inpar_fbi.hpp"
#include "baci_inpar_fluid.hpp"
#include "baci_inpar_fpsi.hpp"
#include "baci_inpar_fs3i.hpp"
#include "baci_inpar_fsi.hpp"
#include "baci_inpar_geometric_search.hpp"
#include "baci_inpar_immersed.hpp"
#include "baci_inpar_io.hpp"
#include "baci_inpar_IO_monitor_structure_dbc.hpp"
#include "baci_inpar_IO_runtime_output.hpp"
#include "baci_inpar_IO_runtime_output_fluid.hpp"
#include "baci_inpar_IO_runtime_output_structure_beams.hpp"
#include "baci_inpar_IO_runtime_vtk_output_structure.hpp"
#include "baci_inpar_IO_runtime_vtp_output_structure.hpp"
#include "baci_inpar_levelset.hpp"
#include "baci_inpar_lubrication.hpp"
#include "baci_inpar_mor.hpp"
#include "baci_inpar_mortar.hpp"
#include "baci_inpar_particle.hpp"
#include "baci_inpar_pasi.hpp"
#include "baci_inpar_plasticity.hpp"
#include "baci_inpar_poroelast.hpp"
#include "baci_inpar_porofluidmultiphase.hpp"
#include "baci_inpar_poromultiphase.hpp"
#include "baci_inpar_poromultiphase_scatra.hpp"
#include "baci_inpar_poroscatra.hpp"
#include "baci_inpar_problemtype.hpp"
#include "baci_inpar_rebalance.hpp"
#include "baci_inpar_s2i.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_inpar_searchtree.hpp"
#include "baci_inpar_solver.hpp"
#include "baci_inpar_solver_nonlin.hpp"
#include "baci_inpar_ssi.hpp"
#include "baci_inpar_ssti.hpp"
#include "baci_inpar_sti.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_inpar_thermo.hpp"
#include "baci_inpar_tsi.hpp"
#include "baci_inpar_turbulence.hpp"
#include "baci_inpar_volmortar.hpp"
#include "baci_inpar_wear.hpp"
#include "baci_inpar_xfem.hpp"
#include "baci_io_pstream.hpp"

#include <Teuchos_any.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintValidParameters()
{
  Teuchos::RCP<const Teuchos::ParameterList> list = INPUT::ValidParameters();
  list->print(std::cout,
      Teuchos::ParameterList::PrintOptions().showDoc(true).showFlags(false).indent(4).showTypes(
          false));
}


/*----------------------------------------------------------------------*/
//! Print help message
/*----------------------------------------------------------------------*/
void PrintHelpMessage()
{
#ifdef BACI_DEBUG
  char baci_build[] = "baci-debug";
#else
  char baci_build[] = "baci-release";
#endif

  std::cout << "NAME\n"
            << "\t" << baci_build << " - simulate just about anything\n"
            << "\n"
            << "SYNOPSIS\n"
            << "\t" << baci_build
            << " [-h | --help] [-p | --parameters] [-d | --datfile] [-ngroup=<x>] \\ \n"
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
            << "\t\t\"diffgroupx\" can be used to compare results from separate but parallel baci "
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
            << "\t\tBaci waits at the beginning for keyboard input. \n"
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::PrintDatHeader(
    std::ostream& stream, const Teuchos::ParameterList& list, std::string parentname, bool comment)
{
  // prevent invalid ordering of parameters caused by alphabetical output:
  // in the first run, print out all list elements that are not a sublist
  // in the second run, do the recursive call for all the sublists in the list
  for (int j = 0; j < 2; ++j)
  {
    for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i)
    {
      const Teuchos::ParameterEntry& entry = list.entry(i);
      if (entry.isList() && j == 0) continue;
      if ((!entry.isList()) && j == 1) continue;
      const std::string& name = list.name(i);
      Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = entry.validator();

      if (comment)
      {
        stream << "//" << '\n';

        std::string doc = entry.docString();
        if (doc != "")
        {
          Teuchos::StrUtils::printLines(stream, "// ", doc);
        }
      }

      if (entry.isList())
      {
        std::string secname = parentname;
        if (secname != "") secname += "/";
        secname += name;
        unsigned l = secname.length();
        stream << "--" << std::string(std::max<int>(65 - l, 0), '-');
        stream << secname << '\n';
        PrintDatHeader(stream, list.sublist(name), secname, comment);
      }
      else
      {
        if (comment)
          if (validator != Teuchos::null)
          {
            Teuchos::RCP<const Teuchos::Array<std::string>> values = validator->validStringValues();
            if (values != Teuchos::null)
            {
              unsigned len = 0;
              for (int i = 0; i < (int)values->size(); ++i)
              {
                len += (*values)[i].length() + 1;
              }
              if (len < 74)
              {
                stream << "//     ";
                for (int i = 0; i < static_cast<int>(values->size()) - 1; ++i)
                {
                  stream << (*values)[i] << ",";
                }
                stream << (*values)[values->size() - 1] << '\n';
              }
              else
              {
                for (int i = 0; i < (int)values->size(); ++i)
                {
                  stream << "//     " << (*values)[i] << '\n';
                }
              }
            }
          }
        const Teuchos::any& v = entry.getAny(false);
        stream << name;
        unsigned l = name.length();
        stream << std::string(std::max<int>(31 - l, 0), ' ');
        if (NeedToPrintEqualSign(list)) stream << " =";
        stream << ' ' << v << '\n';
      }
    }
  }
}

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintDefaultDatHeader()
{
  Teuchos::RCP<const Teuchos::ParameterList> list = INPUT::ValidParameters();
  INPUT::PrintDatHeader(std::cout, *list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::PrintDefaultParameters(IO::Pstream& stream, const Teuchos::ParameterList& list)
{
  bool hasDefault = false;
  for (Teuchos::ParameterList::ConstIterator i = list.begin(); i != list.end(); ++i)
  {
    const Teuchos::ParameterEntry& entry = list.entry(i);
    if (entry.isDefault())
    {
      if (not hasDefault)
      {
        hasDefault = true;
        stream << "default parameters in list '" << list.name() << "':\n";
      }
      const Teuchos::any& v = entry.getAny(false);
      int l = list.name(i).length();
      stream << "    " << list.name(i);
      stream << std::string(std::max<int>(31 - l, 0), ' ');
      stream << ' ' << v << "\n";
    }
  }
  if (hasDefault) stream << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::BoolParameter(std::string const& paramName, std::string const& value,
    std::string const& docString, Teuchos::ParameterList* paramList)
{
  Teuchos::Array<std::string> yesnotuple =
      Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = Teuchos::tuple<int>(true, false, true, false, true, false);
  Teuchos::setStringToIntegralParameter<int>(
      paramName, value, docString, yesnotuple, yesnovalue, paramList);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::IntParameter(std::string const& paramName, int const value,
    std::string const& docString, Teuchos::ParameterList* paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowInt(true);
  Teuchos::setIntParameter(paramName, value, docString, paramList, validator);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::DoubleParameter(std::string const& paramName, double const& value,
    std::string const& docString, Teuchos::ParameterList* paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowDouble(true);
  validator.allowInt(true);
  Teuchos::setDoubleParameter(paramName, value, docString, paramList, validator);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPUT::StringParameter(std::string const& paramName, std::string const& value,
    std::string const& docString, Teuchos::ParameterList* paramList)
{
  Teuchos::RCP<Teuchos::StringValidator> validator = Teuchos::rcp(new Teuchos::StringValidator());
  paramList->set(paramName, value, docString, validator);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> INPUT::ValidParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& discret = list->sublist("DISCRETISATION", false, "");

  IntParameter("NUMFLUIDDIS", 1, "Number of meshes in fluid field", &discret);
  IntParameter("NUMSTRUCDIS", 1, "Number of meshes in structural field", &discret);
  IntParameter("NUMALEDIS", 1, "Number of meshes in ale field", &discret);
  IntParameter("NUMARTNETDIS", 1, "Number of meshes in arterial network field", &discret);
  IntParameter("NUMTHERMDIS", 1, "Number of meshes in thermal field", &discret);
  IntParameter("NUMAIRWAYSDIS", 1, "Number of meshes in reduced dimensional airways network field",
      &discret);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& size = list->sublist("PROBLEM SIZE", false, "");

  IntParameter("DIM", 3, "2d or 3d problem", &size);

  // deactivate all the follwing (unused) parameters one day
  // they are nice as general info in the input file but should not
  // read into a parameter list. Misuse is possible
  IntParameter("ELEMENTS", 0, "Total number of elements", &size);
  IntParameter("NODES", 0, "Total number of nodes", &size);
  IntParameter("NPATCHES", 0, "number of nurbs patches", &size);
  IntParameter("MATERIALS", 0, "number of materials", &size);
  IntParameter("NUMDF", 3, "maximum number of degrees of freedom", &size);

  INPAR::PROBLEMTYPE::SetValidParameters(list);

  /*----------------------------------------------------------------------*/
  /* Finally call the problem-specific SetValidParameter functions        */
  /*----------------------------------------------------------------------*/

  INPAR::STR::SetValidParameters(list);
  INPAR::IO::SetValidParameters(list);
  INPAR::IO_MONITOR_STRUCTURE_DBC::SetValidParameters(list);
  INPAR::IO_RUNTIME_OUTPUT::SetValidParameters(list);
  INPAR::IO_RUNTIME_VTP_STRUCTURE::SetValidParameters(list);
  INPAR::MORTAR::SetValidParameters(list);
  INPAR::CONTACT::SetValidParameters(list);
  INPAR::VOLMORTAR::SetValidParameters(list);
  INPAR::WEAR::SetValidParameters(list);
  INPAR::IO_RUNTIME_OUTPUT::FLUID::SetValidParameters(list);
  INPAR::IO_RUNTIME_OUTPUT::STRUCTURE::SetValidParameters(list);
  INPAR::IO_RUNTIME_OUTPUT::BEAMS::SetValidParameters(list);
  INPAR::BEAMCONTACT::SetValidParameters(list);
  INPAR::BEAMPOTENTIAL::SetValidParameters(list);
  INPAR::BEAMINTERACTION::SetValidParameters(list);
  INPAR::BROWNIANDYN::SetValidParameters(list);

  INPAR::PLASTICITY::SetValidParameters(list);

  INPAR::THR::SetValidParameters(list);
  INPAR::TSI::SetValidParameters(list);

  INPAR::FLUID::SetValidParameters(list);
  INPAR::LOMA::SetValidParameters(list);
  INPAR::CUT::SetValidParameters(list);
  INPAR::XFEM::SetValidParameters(list);

  INPAR::LUBRICATION::SetValidParameters(list);
  INPAR::SCATRA::SetValidParameters(list);
  INPAR::LEVELSET::SetValidParameters(list);
  INPAR::ELCH::SetValidParameters(list);
  INPAR::EP::SetValidParameters(list);
  INPAR::STI::SetValidParameters(list);

  INPAR::S2I::SetValidParameters(list);
  INPAR::FS3I::SetValidParameters(list);
  INPAR::POROELAST::SetValidParameters(list);
  INPAR::PORO_SCATRA::SetValidParameters(list);
  INPAR::POROMULTIPHASE::SetValidParameters(list);
  INPAR::POROMULTIPHASESCATRA::SetValidParameters(list);
  INPAR::POROFLUIDMULTIPHASE::SetValidParameters(list);
  INPAR::EHL::SetValidParameters(list);
  INPAR::SSI::SetValidParameters(list);
  INPAR::SSTI::SetValidParameters(list);
  INPAR::ALE::SetValidParameters(list);
  INPAR::FSI::SetValidParameters(list);

  INPAR::ARTDYN::SetValidParameters(list);
  INPAR::ARTNET::SetValidParameters(list);
  INPAR::BIOFILM::SetValidParameters(list);
  INPAR::REDAIRWAYS::SetValidParameters(list);
  INPAR::CARDIOVASCULAR0D::SetValidParameters(list);
  INPAR::IMMERSED::SetValidParameters(list);
  INPAR::FPSI::SetValidParameters(list);
  INPAR::FBI::SetValidParameters(list);

  INPAR::PARTICLE::SetValidParameters(list);

  INPAR::MOR::SetValidParameters(list);

  INPAR::ELEMAG::SetValidParameters(list);

  INPAR::GEO::SetValidParameters(list);
  INPAR::BINSTRATEGY::SetValidParameters(list);
  INPAR::GEOMETRICSEARCH::SetValidParameters(list);
  INPAR::PASI::SetValidParameters(list);

  INPAR::REBALANCE::SetValidParameters(list);
  INPAR::SOLVER::SetValidParameters(list);
  INPAR::NLNSOL::SetValidParameters(list);

  return list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool INPUT::NeedToPrintEqualSign(const Teuchos::ParameterList& list)
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
BACI_NAMESPACE_CLOSE
