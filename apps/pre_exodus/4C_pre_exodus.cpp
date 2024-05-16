/*----------------------------------------------------------------------*/
/*! \file

\brief preprocessor for exodusII format

\level 1


Pre_exodus contains classes to open and preprocess exodusII files into the
drt of 4C. It uses the "valid-parameters"-list defined in 4C for preparing
a up-to-date 4C header and another file specifying element and boundary
specifications based on "valid-conditions". As result either a preliminary
input file set is suggestioned, or the well-known .dat file is created.
Addionally, specify an already existing 4C input file in order to validate
its parameters and conditions.

*/
/*----------------------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_inpar_validconditions.hpp"
#include "4C_inpar_validmaterials.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_lib_conditiondefinition.hpp"
#include "4C_lib_elementdefinition.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_pre_exodus_readbc.hpp"
#include "4C_pre_exodus_reader.hpp"
#include "4C_pre_exodus_validate.hpp"
#include "4C_pre_exodus_writedat.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>


using namespace FourC;

namespace
{
  /*----------------------------------------------------------------------*/
  /* create default bc file                                               */
  /*----------------------------------------------------------------------*/
  int CreateDefaultBCFile(EXODUS::Mesh& mymesh)
  {
    using namespace FourC;

    std::string defaultbcfilename = "default.bc";
    std::cout << "found no BC specification file --> creating " << defaultbcfilename << std::endl;

    // open default bc specification file
    std::ofstream defaultbc(defaultbcfilename.c_str());
    if (!defaultbc) FOUR_C_THROW("failed to open file: %s", defaultbcfilename.c_str());

    // write mesh verbosely
    defaultbc << "----------- Mesh contents -----------" << std::endl << std::endl;
    mymesh.Print(defaultbc, false);

    // give examples for element and boundary condition syntax
    defaultbc << "---------- Syntax examples ----------" << std::endl
              << std::endl
              << "Element Block, named: " << std::endl
              << "of Shape: TET4" << std::endl
              << "has 9417816 Elements" << std::endl
              << "'*eb0=\"ELEMENT\"'" << std::endl
              << "sectionname=\"FLUID\"" << std::endl
              << "description=\"MAT 1 NA Euler\"" << std::endl
              << "elementname=\"FLUID\" \n"
              << std::endl
              << "Element Block, named: " << std::endl
              << "of Shape: HEX8" << std::endl
              << "has 9417816 Elements" << std::endl
              << "'*eb0=\"ELEMENT\"'" << std::endl
              << "sectionname=\"STRUCTURE\"" << std::endl
              << "description=\"MAT 1 KINEM nonlinear EAS none\"" << std::endl
              << "elementname=\"SOLIDH8\" \n"
              << std::endl
              << "Node Set, named:" << std::endl
              << "Property Name: INFLOW" << std::endl
              << "has 45107 Nodes" << std::endl
              << "'*ns0=\"CONDITION\"'" << std::endl
              << "sectionname=\"DESIGN SURF DIRICH CONDITIONS\"" << std::endl
              << "description=\"NUMDOF 6 ONOFF 1 1 1 0 0 0 VAL 2.0 0.0 0.0 0.0 0.0 0.0 FUNCT 1 0 0 "
                 "0 0 0\""
              << std::endl
              << std::endl;

    defaultbc << "MIND that you can specify a condition also on an ElementBlock, just replace "
                 "'ELEMENT' with 'CONDITION'"
              << std::endl;
    defaultbc << "The 'E num' in the dat-file depends on the order of the specification below"
              << std::endl;
    defaultbc << "------------------------------------------------BCSPECS" << std::endl
              << std::endl;

    // write ElementBlocks with specification proposal
    const std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> myblocks = mymesh.GetElementBlocks();
    std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator it;
    for (it = myblocks.begin(); it != myblocks.end(); ++it)
    {
      it->second->Print(defaultbc);
      defaultbc
          << "*eb" << it->first << "=\"ELEMENT\"" << std::endl
          << "sectionname=\"\"" << std::endl
          << "description=\"\"" << std::endl
          << "elementname=\"\""
          << std::endl
          //<<"elementshape=\""
          //<< CORE::FE::CellTypeToString(PreShapeToDrt(it->second.GetShape()))<<"\""<<std::endl
          << std::endl;
    }

    // write NodeSets with specification proposal
    const std::map<int, EXODUS::NodeSet> mynodesets = mymesh.GetNodeSets();
    std::map<int, EXODUS::NodeSet>::const_iterator ins;
    for (ins = mynodesets.begin(); ins != mynodesets.end(); ++ins)
    {
      ins->second.Print(defaultbc);
      defaultbc << "*ns" << ins->first << "=\"CONDITION\"" << std::endl
                << "sectionname=\"\"" << std::endl
                << "description=\"\"" << std::endl
                << std::endl;
    }

    // write SideSets with specification proposal
    const std::map<int, EXODUS::SideSet> mysidesets = mymesh.GetSideSets();
    std::map<int, EXODUS::SideSet>::const_iterator iss;
    for (iss = mysidesets.begin(); iss != mysidesets.end(); ++iss)
    {
      iss->second.Print(defaultbc);
      defaultbc << "*ss" << iss->first << "=\"CONDITION\"" << std::endl
                << "sectionname=\"\"" << std::endl
                << "description=\"\"" << std::endl
                << std::endl;
    }

    // print validconditions as proposal
    defaultbc << "-----------------------------------------VALIDCONDITIONS" << std::endl;
    Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>> condlist =
        INPUT::ValidConditions();
    INPUT::PrintEmptyConditionDefinitions(defaultbc, *condlist);

    // print valid element lines as proposal (parobjects have to be registered for doing this!)
    defaultbc << std::endl << std::endl;
    INPUT::ElementDefinition ed;
    ed.PrintElementDatHeaderToStream(defaultbc);

    // close default bc specification file
    if (defaultbc.is_open()) defaultbc.close();

    return 0;
  }
}  // namespace

/**
 *
 * Pre_exodus contains classes to open and preprocess exodusII files into the
 * discretization of 4C. It uses the "valid-parameters"-list defined in 4C for preparing
 * a up-to-date 4C header and another file specifying element and boundary
 * specifications. As result either a preliminary input file set is suggestioned,
 * or the well-known .dat file is created.
 *
 */
int main(int argc, char** argv)
{
  using namespace FourC;

  // communication
  MPI_Init(&argc, &argv);

  GlobalLegacyModuleCallbacks().RegisterParObjectTypes();

  // create a problem instance
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  // create default communicators
  Teuchos::RCP<CORE::COMM::Communicators> communicators = CORE::COMM::CreateComm({});
  GLOBAL::Problem::Instance()->SetCommunicators(communicators);
  Teuchos::RCP<Epetra_Comm> comm = communicators->GlobalComm();

  try
  {
    if ((comm->NumProc() > 1)) FOUR_C_THROW("Using more than one processor is not supported.");

    std::string exofile;
    std::string bcfile;
    std::string headfile;
    std::string datfile;
    std::string cline;

    bool twodim = false;

    // related to quad->tri conversion
    bool quadtri = false;

    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString(
        "This preprocessor converts Exodus2 files to the proprietary FourC format\n");
    My_CLP.throwExceptions(false);
    My_CLP.setOption("exo", &exofile, "exodus file to open");
    My_CLP.setOption("bc", &bcfile, "bc's and ele's file to open");
    My_CLP.setOption("head", &headfile, "4C header file to open");
    My_CLP.setOption("dat", &datfile, "output .dat file name [defaults to exodus file name]");

    // switch for genarating a 2d .dat - file
    My_CLP.setOption("d2", "d3", &twodim, "space dimensions in .dat-file: d2: 2D, d3: 3D");

    // check for quad->tri conversion
    My_CLP.setOption(
        "quadtri", "noquadtri", &quadtri, "transform quads to tris by cutting in two halves");

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc, argv);

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
      // free the global problem instance
      problem->Done();
      return 0;
    }
    if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
      FOUR_C_THROW("CommandLineProcessor reported an error");
    }

    if (datfile != "")
    {
      const std::string basename = datfile.substr(0, datfile.find_last_of(".")) + "_pre";
      IO::cout.setup(
          true, false, false, IO::standard, comm, 0, 0, basename);  // necessary setup of IO::cout
    }
    else
    {
      IO::cout.setup(
          true, false, false, IO::standard, comm, 0, 0, "xxx_pre");  // necessary setup of IO::cout
    }


    /**************************************************************************
     * Start with the preprocessing
     **************************************************************************/
    if (exofile == "")
    {
      if (datfile != "")
      {
        // just validate a given 4C input file
        EXODUS::ValidateInputFile(comm, datfile);
        return 0;
      }
      else
      {
        My_CLP.printHelpMessage(argv[0], std::cout);
        FOUR_C_THROW("No Exodus II file was found");
      }
    }

    // create mesh object based on given exodus II file
    EXODUS::Mesh mymesh(exofile);
    // print infos to std::cout
    mymesh.Print(std::cout);

    /**************************************************************************
     * Edit a existing Mesh, e.g. extrusion of surface
     **************************************************************************/

    // transform quads->tris
    if (quadtri)
    {
      EXODUS::Mesh trimesh = EXODUS::QuadtoTri(mymesh);
      trimesh.WriteMesh("tri_" + exofile);
      exit(0);
    }

    /**************************************************************************
     * Read ControlFile for Boundary and Element descriptions
     **************************************************************************/

    // declare empty vectors for holding "boundary" conditions
    std::vector<EXODUS::ElemDef> eledefs;
    std::vector<EXODUS::CondDef> condefs;

    if (bcfile == "")
    {
      int error = CreateDefaultBCFile(mymesh);
      if (error != 0) FOUR_C_THROW("Creation of default bc-file not successful.");
    }
    else
    {
      // read provided bc-file
      EXODUS::ReadBCFile(bcfile, eledefs, condefs);

      int sum = mymesh.GetNumElementBlocks() + mymesh.GetNumNodeSets() + mymesh.GetNumSideSets();
      int test = eledefs.size() + condefs.size();
      if (test != sum)
        std::cout
            << "Your " << test << " definitions do not match the " << sum
            << " entities in your mesh!" << std::endl
            << "(This is OK, if more than one BC is applied to an entity, e.g in FSI simulations)"
            << std::endl;
    }

    /**************************************************************************
     * Read HeaderFile for 'header' parameters, e.g. solver, dynamic, material
     * or create a default HeaderFile
     **************************************************************************/
    if (headfile == "")
    {
      const std::string defaultheadfilename = "default.head";
      std::cout << "found no header file           --> creating " << defaultheadfilename
                << std::endl;

      // open default header file
      std::ofstream defaulthead(defaultheadfilename.c_str());
      if (!defaulthead) FOUR_C_THROW("failed to open file: %s", defaultheadfilename.c_str());

      // get valid input parameters
      Teuchos::RCP<const Teuchos::ParameterList> list = INPUT::ValidParameters();

      // write default .dat header into file
      std::stringstream prelimhead;
      INPUT::PrintDatHeader(prelimhead, *list);
      std::string headstring = prelimhead.str();
      size_t size_section =
          headstring.find("-------------------------------------------------------PROBLEM SIZE");
      if (size_section != std::string::npos)
      {
        size_t typ_section =
            headstring.find("--------------------------------------------------------PROBLEM TYP");
        headstring.erase(size_section, typ_section - size_section);
      }
      defaulthead << headstring;

      // get valid input materials
      {
        Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>> mlist =
            INPUT::ValidMaterials();
        INPUT::PrintEmptyMaterialDefinitions(defaulthead, *mlist);
      }

      // print cloning material map default lines (right after the materials)
      INPUT::Lines lines = DRT::UTILS::ValidCloningMaterialMapLines();
      lines.Print(defaulthead);

      // print spatial functions
      defaulthead << "-------------------------------------------------------------FUNCT1"
                  << std::endl
                  << "-------------------------------------------------------------FUNCT2"
                  << std::endl
                  << "-------------------------------------------------------------FUNCT3"
                  << std::endl
                  << "-------------------------------------------------------------FUNCT4"
                  << std::endl;
      {
        std::stringstream tmp;
        CORE::UTILS::FunctionManager functionmanager;
        INPUT::Lines flines = functionmanager.ValidFunctionLines();
        flines.Print(tmp);
        std::string tmpstring = tmp.str();
        std::string removeit =
            "--------------------------------------------------------------FUNCT\n";
        size_t pos = tmpstring.find(removeit);
        if (pos != std::string::npos)
        {
          tmpstring.erase(pos, removeit.length());
        }
        defaulthead << tmpstring;
      }

      // default result-test lines
      {
        CORE::UTILS::ResultTestManager resulttestmanager;
        INPUT::Lines lines("RESULT DESCRIPTION",
            "The result of the simulation with respect to specific quantities at concrete points "
            "can be tested against particular values with a given tolerance.");
        GlobalLegacyModuleCallbacks().AttachResultLines(lines);
        lines.Print(defaulthead);
      }

      // close default header file
      if (defaulthead.is_open()) defaulthead.close();
    }

    /**************************************************************************
     * Finally, create and validate the 4C input file
     **************************************************************************/
    if ((headfile != "") && (bcfile != "") && (exofile != ""))
    {
      // set default dat-file name if needed
      if (datfile == "")
      {
        const std::string exofilebasename = exofile.substr(0, exofile.find_last_of("."));
        datfile = exofilebasename + ".dat";
      }

      // screen info
      std::cout << "creating and checking 4C input file       --> " << datfile << std::endl;
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("pre-exodus timer");

      // check for positive Element-Center-Jacobians and otherwise rewind them
      {
        std::cout << "...Ensure positive element jacobians";
        timer->start();
        ValidateMeshElementJacobians(mymesh);
        timer->stop();
        std::cout << "        in...." << timer->totalElapsedTime(true) << " secs" << std::endl;
        timer->reset();
      }

      // in case of periodic boundary conditions :
      // ensure that the two coordinates of two matching nodes,
      // which should be the same are exactly the same
      // in order to keep the Krylov norm below 1e-6 :-)
      // only supported for angle 0.0
      {
        if (PeriodicBoundaryConditionsFound(condefs))
        {
          std::cout << "...Ensure high quality p.b.c.";
          timer->start();
          CorrectNodalCoordinatesForPeriodicBoundaryConditions(mymesh, condefs);
          timer->stop();
          std::cout << "               in...." << timer->totalElapsedTime(true) << " secs"
                    << std::endl;
          timer->reset();
        }
      }

      // write the 4C input file
      {
        if (twodim) mymesh.SetNsd(2);
        std::cout << "...Writing dat-file";
        timer->start();
        EXODUS::WriteDatFile(datfile, mymesh, headfile, eledefs, condefs);
        timer->stop();
        std::cout << "                         in...." << timer->totalElapsedTime(true) << " secs"
                  << std::endl;
        timer->reset();
      }

      // validate the generated 4C input file
      EXODUS::ValidateInputFile(comm, datfile);
    }
  }
  catch (CORE::Exception& err)
  {
    char line[] = "=========================================================================\n";
    std::cout << "\n\n" << line << err.what_with_stacktrace() << "\n" << line << "\n" << std::endl;

    // free the global problem instance and the communicator
    problem->Done();
    comm = Teuchos::null;

#ifdef FOUR_C_ENABLE_CORE_DUMP
    abort();
#endif

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // free the global problem instance
  problem->Done();

  MPI_Finalize();

  return 0;

}  // main.cpp
