/*----------------------------------------------------------------------*/
/*! \file

\brief preprocessor for exodusII format

\level 1


Pre_exodus contains classes to open and preprocess exodusII files into the
drt of Baci. It uses the "valid-parameters"-list defined in Baci for preparing
a up-to-date Baci header and another file specifying element and boundary
specifications based on "valid-conditions". As result either a preliminary
input file set is suggestioned, or the well-known .dat file is created.
Addionally, specify an already existing BACI input file in order to validate
its parameters and conditions.

*/
/*----------------------------------------------------------------------*/

#include "pre_exodus.H"
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/drt_validmaterials.H"
#include "../drt_inpar/drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "../drt_lib/drt_elementdefinition.H"
#include "../drt_lib/drt_parobjectregister.H"
#include "../drt_comm/comm_utils.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_soshextrusion.H"
#include "pre_exodus_writedat.H"
#include "pre_exodus_readbc.H"
#include "pre_exodus_validate.H"
#include "pre_exodus_centerline.H"

#include <Epetra_MpiComm.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  // communication
  MPI_Init(&argc, &argv);

  // create a problem instance
  DRT::Problem* problem = DRT::Problem::Instance();
  // create "dummy" NP group which only sets the correct communicators
  COMM_UTILS::CreateComm(0, NULL);
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(problem->GetNPGroup()->GlobalComm().get(), false);

  try
  {
    if ((comm->NumProc() > 1)) dserror("Using more than one processor is not supported.");

    std::string exofile;
    std::string bcfile;
    std::string headfile;
    std::string datfile;
    std::string cline;

    int printparobjecttypes = 0;

    // related to solid shell extrusion
    double soshthickness = 0.0;
    int soshnumlayer = 1;
    int soshseedid = 0;
    int soshgmsh = -1;
    int concat2loose = 0;
    int diveblocks = 0;

    // related to centerline
    std::vector<double> cline_coordcorr(3);
    double clinedx = 0.0;
    double clinedy = 0.0;
    double clinedz = 0.0;

    bool twodim = false;

    // related to quad->tri conversion
    bool quadtri = false;

    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Preprocessor Exodus2Baci \n");
    My_CLP.throwExceptions(false);
    My_CLP.setOption("exo", &exofile, "exodus file to open");
    My_CLP.setOption("bc", &bcfile, "bc's and ele's file to open");
    My_CLP.setOption("head", &headfile, "baci header file to open");
    My_CLP.setOption("dat", &datfile, "output .dat file name [defaults to exodus file name]");
    // here options related to solid shell extrusion are defined
    My_CLP.setOption("gensosh", &soshthickness, "generate solid-shell body with given thickness");
    My_CLP.setOption("numlayer", &soshnumlayer, "number of layers of generated solid-shell body");
    My_CLP.setOption("diveblocks", &diveblocks,
        "if larger 0 the xxx inner elements of generated layers are put into first eblock, the "
        "rest into second");
    My_CLP.setOption("seedid", &soshseedid, "id where to start extrusion, default is first");
    My_CLP.setOption("gmsh", &soshgmsh, "gmsh output of xxx elements, default off, 0 all eles");
    My_CLP.setOption("concf", &concat2loose,
        "concatenate extruded volume with base, however loose every xxx'th node, default "
        "0=off=fsi");

    // switch for genarating a 2d .dat - file
    My_CLP.setOption("d2", "d3", &twodim, "space dimensions in .dat-file: d2: 2D, d3: 3D");

    // centerline related
    My_CLP.setOption("cline", &cline,
        "generate local element coordinate systems based on centerline file, or mesh line (set to "
        "'mesh'");
    My_CLP.setOption("clinedx", &clinedx, "move centerline coords to align with mesh: delta x");
    My_CLP.setOption("clinedy", &clinedy, "move centerline coords to align with mesh: delta y");
    My_CLP.setOption("clinedz", &clinedz, "move centerline coords to align with mesh: delta z");
    std::map<int, std::map<int, std::vector<std::vector<double>>>> elecenterlineinfo;

    // check for quad->tri conversion
    My_CLP.setOption(
        "quadtri", "noquadtri", &quadtri, "transform quads to tris by cutting in two halves");

    // print parobject types (needed for making automatic object registration working)
    My_CLP.setOption("printparobjecttypes", &printparobjecttypes,
        "print names of parobject types (registration hack)");

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc, argv);

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
      // free the global problem instance
      problem->Done();
      return 0;
    }
    if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    {
      dserror("CommandLineProcessor reported an error");
    }

    if (printparobjecttypes)
    {
      // hack so that the parobject types are registered!!!
      PrintParObjectList();

      // free the global problem instance
      problem->Done();
      return 0;
    }

    // create error file (enforce the file opening!)
    if (datfile != "")
    {
      const std::string basename = datfile.substr(0, datfile.find_last_of(".")) + "_pre";
      IO::cout.setup(
          true, false, false, IO::standard, comm, 0, 0, basename);  // necessary setup of IO::cout
      problem->OpenErrorFile(*comm, basename, true);
    }
    else
    {
      IO::cout.setup(
          true, false, false, IO::standard, comm, 0, 0, "xxx_pre");  // necessary setup of IO::cout
      problem->OpenErrorFile(*comm, "xxx_pre", true);
    }

    // centerline related: transfer separate doubles into vector
    cline_coordcorr[0] = clinedx;
    cline_coordcorr[1] = clinedy;
    cline_coordcorr[2] = clinedz;

    /**************************************************************************
     * Start with the preprocessing
     **************************************************************************/
    if (exofile == "")
    {
      if (datfile != "")
      {
        // just validate a given BACI input file
        EXODUS::ValidateInputFile(comm, datfile);
        return 0;
      }
      else
      {
        My_CLP.printHelpMessage(argv[0], std::cout);
        dserror("No Exodus II file was found");
      }
    }

    // create mesh object based on given exodus II file
    EXODUS::Mesh mymesh(exofile);
    // print infos to std::cout
    mymesh.Print(std::cout);

    /**************************************************************************
     * Edit a existing Mesh, e.g. extrusion of surface
     **************************************************************************/

    // generate solid shell extrusion based on exodus file
    if (soshthickness != 0.0)
    {
      if (exofile == "") dserror("no exofile specified for extrusion");
      if (soshnumlayer <= diveblocks)
        dserror(
            "number of layers and inner-layer elements mismatch, check if numlayer>diveblocks!");
      EXODUS::Mesh mysosh = EXODUS::SolidShellExtrusion(mymesh, soshthickness, soshnumlayer,
          soshseedid, soshgmsh, concat2loose, diveblocks, cline, cline_coordcorr);
      mysosh.WriteMesh("extr_" + exofile);

      exit(0);
    }

    // generate local element coordinate systems based on centerline
    if (cline != "")
    {
      elecenterlineinfo = EleCenterlineInfo(cline, mymesh, cline_coordcorr);
    }

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
    std::vector<EXODUS::elem_def> eledefs;
    std::vector<EXODUS::cond_def> condefs;

    if (bcfile == "")
    {
      int error = EXODUS::CreateDefaultBCFile(mymesh);
      if (error != 0) dserror("Creation of default bc-file not successful.");
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
      if (!defaulthead) dserror("failed to open file: %s", defaultheadfilename.c_str());

      // get valid input parameters
      Teuchos::RCP<const Teuchos::ParameterList> list = DRT::INPUT::ValidParameters();

      // write default .dat header into file
      std::stringstream prelimhead;
      DRT::INPUT::PrintDatHeader(prelimhead, *list);
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
        Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::MaterialDefinition>>> mlist =
            DRT::INPUT::ValidMaterials();
        DRT::INPUT::PrintEmptyMaterialDefinitions(defaulthead, *mlist, false);
      }

      // print cloning material map default lines (right after the materials)
      Teuchos::RCP<DRT::INPUT::Lines> lines = DRT::UTILS::ValidCloningMaterialMapLines();
      lines->Print(defaulthead);

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
        DRT::UTILS::FunctionManager functionmanager;
        Teuchos::RCP<DRT::INPUT::Lines> flines = functionmanager.ValidFunctionLines();
        flines->Print(tmp);
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
        DRT::ResultTestManager resulttestmanager;
        Teuchos::RCP<DRT::INPUT::Lines> lines = resulttestmanager.ValidResultLines();
        lines->Print(defaulthead);
      }

      // close default header file
      if (defaulthead.is_open()) defaulthead.close();
    }

    /**************************************************************************
     * Finally, create and validate the BACI input file
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
      std::cout << "creating and checking BACI input file       --> " << datfile << std::endl;
      Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("pre-exodus timer");

      // check for positive Element-Center-Jacobians and otherwise rewind them
      {
        std::cout << "...Ensure positive element jacobians";
        timer->start();
        ValidateMeshElementJacobians(mymesh);
        timer->stop();
        std::cout << "        in...." << timer->totalElapsedTime() << " secs" << std::endl;
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
          std::cout << "               in...." << timer->totalElapsedTime() << " secs" << std::endl;
          timer->reset();
        }
      }

      // write the BACI input file
      {
        if (twodim) mymesh.SetNsd(2);
        std::cout << "...Writing dat-file";
        timer->start();
        EXODUS::WriteDatFile(datfile, mymesh, headfile, eledefs, condefs, elecenterlineinfo);
        timer->stop();
        std::cout << "                         in...." << timer->totalElapsedTime() << " secs"
                  << std::endl;
        timer->reset();
      }

      // validate the generated BACI input file
      EXODUS::ValidateInputFile(comm, datfile);
    }
  }
  catch (std::runtime_error& err)
  {
    char line[] = "=========================================================================\n";
    std::cout << "\n\n" << line << err.what() << "\n" << line << "\n" << std::endl;

    // free the global problem instance and the communicator
    problem->Done();
    comm = Teuchos::null;

#ifdef DSERROR_DUMP
    abort();
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
    exit(1);
#endif
  }

  // free the global problem instance
  problem->Done();

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;

}  // main.cpp


/*----------------------------------------------------------------------*/
/* create default bc file                                               */
/*----------------------------------------------------------------------*/
int EXODUS::CreateDefaultBCFile(EXODUS::Mesh& mymesh)
{
  std::string defaultbcfilename = "default.bc";
  std::cout << "found no BC specification file --> creating " << defaultbcfilename << std::endl;

  // open default bc specification file
  std::ofstream defaultbc(defaultbcfilename.c_str());
  if (!defaultbc) dserror("failed to open file: %s", defaultbcfilename.c_str());

  // write mesh verbosely
  defaultbc << "----------- Mesh contents -----------" << std::endl << std::endl;
  mymesh.Print(defaultbc, false);

  // give examples for element and boundary condition syntax
  defaultbc
      << "---------- Syntax examples ----------" << std::endl
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
      << "description=\"NUMDOF 6 ONOFF 1 1 1 0 0 0 VAL 2.0 0.0 0.0 0.0 0.0 0.0 FUNCT 1 0 0 0 0 0\""
      << std::endl
      << std::endl;

  defaultbc << "MIND that you can specify a condition also on an ElementBlock, just replace "
               "'ELEMENT' with 'CONDITION'"
            << std::endl;
  defaultbc << "The 'E num' in the dat-file depends on the order of the specification below"
            << std::endl;
  defaultbc << "------------------------------------------------BCSPECS" << std::endl << std::endl;

  // write ElementBlocks with specification proposal
  const std::map<int, Teuchos::RCP<EXODUS::ElementBlock>> myblocks = mymesh.GetElementBlocks();
  std::map<int, Teuchos::RCP<EXODUS::ElementBlock>>::const_iterator it;
  for (it = myblocks.begin(); it != myblocks.end(); ++it)
  {
    it->second->Print(defaultbc);
    defaultbc << "*eb" << it->first << "=\"ELEMENT\"" << std::endl
              << "sectionname=\"\"" << std::endl
              << "description=\"\"" << std::endl
              << "elementname=\"\""
              << std::endl
              //<<"elementshape=\""
              //<< DRT::DistypeToString(PreShapeToDrt(it->second.GetShape()))<<"\""<<std::endl
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
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>> condlist =
      DRT::INPUT::ValidConditions();
  DRT::INPUT::PrintEmptyConditionDefinitions(defaultbc, *condlist, false);

  // print valid element lines as proposal (parobjects have to be registered for doing this!)
  defaultbc << std::endl << std::endl;
  DRT::INPUT::ElementDefinition ed;
  ed.PrintElementDatHeaderToStream(defaultbc);

  // close default bc specification file
  if (defaultbc.is_open()) defaultbc.close();

  return 0;
}
