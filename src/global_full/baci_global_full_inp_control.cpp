/*----------------------------------------------------------------------*/
/*! \file

\brief Global control routine of baci

\level 0

*----------------------------------------------------------------------*/

#include "baci_comm_utils.H"
#include "baci_global_legacy_module.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_io_inputreader.H"
#include "baci_io_pstream.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"

#include <utility>

void SetupParallelOutput(
    std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group);

/*----------------------------------------------------------------------*
  | general input of the problem to be solved              m.gee 10/06  |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret(
    std::string& inputfile_name, std::string& outputfile_kenner, std::string& restartfile_kenner)
{
  using namespace BACI;

  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetCommunicators()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetCommunicators()->GlobalComm();
  int group = problem->GetCommunicators()->GroupId();
  CORE::COMM::NestedParallelismType npType = problem->GetCommunicators()->NpType();



  // and now the actual reading
  INPUT::DatFileReader reader(inputfile_name, lcomm);

  problem->ReadParameter(reader);

  SetupParallelOutput(outputfile_kenner, lcomm, group);

  // create control file for output and read restart data if required
  problem->OpenControlFile(*lcomm, inputfile_name, outputfile_kenner, restartfile_kenner);

  // input of materials
  problem->ReadMaterials(reader);

  // input of materials
  problem->ReadContactConstitutiveLaws(reader);

  // input of materials of cloned fields (if needed)
  problem->ReadCloningMaterialMap(reader);

  {
    CORE::UTILS::FunctionManager function_manager;
    BACI::GlobalLegacyModuleCallbacks().AttachFunctionDefinitions(function_manager);
    function_manager.ReadInput(reader);
    problem->SetFunctionManager(std::move(function_manager));
  }

  problem->ReadResult(reader);

  // input of particles
  problem->ReadParticles(reader);

  switch (npType)
  {
    case CORE::COMM::NestedParallelismType::no_nested_parallelism:
    case CORE::COMM::NestedParallelismType::every_group_read_dat_file:
    case CORE::COMM::NestedParallelismType::separate_dat_files:
      // input of fields
      problem->ReadFields(reader);

      // read all types of geometry related conditions (e.g. boundary conditions)
      // Also read time and space functions and local coord systems
      problem->ReadConditions(reader);

      // read all knot information for isogeometric analysis
      // and add it to the (derived) nurbs discretization
      problem->ReadKnots(reader);
      break;
    case CORE::COMM::NestedParallelismType::copy_dat_file:
      // group 0 only reads discretization etc
      if (group == 0)
      {
        // input of fields
        problem->ReadFields(reader);

        // read all types of geometry related conditions (e.g. boundary conditions)
        // Also read time and space functions and local coord systems
        problem->ReadConditions(reader);

        // read all knot information for isogeometric analysis
        // and add it to the (derived) nurbs discretization
        problem->ReadKnots(reader);
      }
      gcomm->Barrier();
      // group 0 broadcasts the discretizations to the other groups
      DRT::BroadcastDiscretizations(*problem);
      gcomm->Barrier();
      break;
    default:
      dserror("nptype (nested parallelity type) not recognized");
      break;
  }

  // all reading is done at this point!

  if (lcomm->MyPID() == 0 && npType != CORE::COMM::NestedParallelismType::copy_dat_file)
    problem->WriteInputParameters();


  // before we destroy the reader we want to know about unused sections
  if (npType == CORE::COMM::NestedParallelismType::copy_dat_file)
  {
    if (group == 0) reader.PrintUnknownSections();
  }
  else
    reader.PrintUnknownSections();
}  // end of ntainp_ccadiscret()


/*----------------------------------------------------------------------*
  | setup parallel output                                  ghamm 11/12  |
 *----------------------------------------------------------------------*/
void SetupParallelOutput(std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group)
{
  using namespace BACI;

  // configure the parallel output environment
  const Teuchos::ParameterList& io = DRT::Problem::Instance()->IOParams();
  bool screen = INPUT::IntegralValue<int>(io, "WRITE_TO_SCREEN");
  bool file = INPUT::IntegralValue<int>(io, "WRITE_TO_FILE");
  bool preGrpID = INPUT::IntegralValue<int>(io, "PREFIX_GROUP_ID");
  int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
  auto level = INPUT::IntegralValue<IO::verbositylevel>(io, "VERBOSITY");

  IO::cout.setup(screen, file, preGrpID, level, std::move(lcomm), oproc, group, outputfile_kenner);
}
