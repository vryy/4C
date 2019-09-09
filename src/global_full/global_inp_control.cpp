/*----------------------------------------------------------------------*/
/*! \file

\brief Global control routine of baci

\level 0

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_io/io_pstream.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"

void SetupParallelOutput(
    std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group);

/*----------------------------------------------------------------------*
  | general input of the problem to be solved              m.gee 10/06  |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret(
    std::string& inputfile_name, std::string& outputfile_kenner, std::string& restartfile_kenner)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetNPGroup()->GlobalComm();
  int group = problem->GetNPGroup()->GroupId();
  NP_TYPE npType = problem->GetNPGroup()->NpType();



  // and now the actual reading
  DRT::INPUT::DatFileReader reader(inputfile_name, lcomm);

  problem->ReadParameter(reader);

  SetupParallelOutput(outputfile_kenner, lcomm, group);

  // create error files
  problem->OpenErrorFile(*lcomm, outputfile_kenner);

  // create control file for output and read restart data if required
  problem->OpenControlFile(*lcomm, inputfile_name, outputfile_kenner, restartfile_kenner);

  // input of materials
  problem->ReadMaterials(reader);

  // input of materials
  problem->ReadContactConstitutiveLaws(reader);

  // input of materials of cloned fields (if needed)
  problem->ReadCloningMaterialMap(reader);

  // input of time curves, functions and result tests
  problem->ReadTimeFunctionResult(reader);

  // input of particles
  problem->ReadParticles(reader);

  switch (npType)
  {
    case no_nested_parallelism:
    case every_group_read_dat_file:
    case separate_dat_files:
      // input of fields
      problem->ReadFields(reader);

      // read all types of geometry related conditions (e.g. boundary conditions)
      // Also read time and space functions and local coord systems
      problem->ReadConditions(reader);

      // read all knot information for isogeometric analysis
      // and add it to the (derived) nurbs discretization
      problem->ReadKnots(reader);
      break;
    case copy_dat_file:
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
      COMM_UTILS::BroadcastDiscretizations();
      gcomm->Barrier();
      break;
    default:
      dserror("nptype (nested parallelity type) not recognized");
      break;
  }

  // all reading is done at this point!

  if (lcomm->MyPID() == 0 && npType != copy_dat_file) problem->WriteInputParameters();

  /// dump input file contents to error file (DEBUG-mode only)
  if (npType == copy_dat_file)
  {
    if (group == 0) reader.DumpInput();
  }
  else
    reader.DumpInput();

  // before we destroy the reader we want to know about unused sections
  if (npType == copy_dat_file)
  {
    if (group == 0) reader.PrintUnknownSections();
  }
  else
    reader.PrintUnknownSections();

  return;
}  // end of ntainp_ccadiscret()


/*----------------------------------------------------------------------*
  | setup parallel output                                  ghamm 11/12  |
 *----------------------------------------------------------------------*/
void SetupParallelOutput(std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group)
{
  // configure the parallel output environment
  const Teuchos::ParameterList& io = DRT::Problem::Instance()->IOParams();
  bool screen = DRT::INPUT::IntegralValue<int>(io, "WRITE_TO_SCREEN");
  bool file = DRT::INPUT::IntegralValue<int>(io, "WRITE_TO_FILE");
  bool preGrpID = DRT::INPUT::IntegralValue<int>(io, "PREFIX_GROUP_ID");
  int oproc = io.get<int>("LIMIT_OUTP_TO_PROC");
  IO::verbositylevel level = DRT::INPUT::IntegralValue<IO::verbositylevel>(io, "VERBOSITY");

  IO::cout.setup(screen, file, preGrpID, level, lcomm, oproc, group, outputfile_kenner);

  return;
}
