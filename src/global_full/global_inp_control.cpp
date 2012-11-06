/*!----------------------------------------------------------------------
\file global_inp_control.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_io/io_pstream.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"


/*----------------------------------------------------------------------*
  | general input of the problem to be solved              m.gee 10/06  |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret(
  std::string& inputfile_name,
  std::string& outputfile_kenner,
  std::string& restartfile_kenner
  )
{

  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetNPGroup()->GlobalComm();
  int group = problem->GetNPGroup()->GroupId();
  NP_TYPE npType = problem->GetNPGroup()->NpType();

  // create error files
  problem->OpenErrorFile(*lcomm,outputfile_kenner);

  // and now the actual reading
  DRT::INPUT::DatFileReader reader(inputfile_name,
                                   lcomm);

  problem->ReadParameter(reader);

  // input of materials
  problem->ReadMaterials(reader);
  
  // input of time curves, functions and result tests
  problem->ReadTimeFunctionResult(reader);

  switch(npType)
  {
  case no_nested_parallelism:
  case every_group_read_dat_file:
  case separate_dat_files:
    // input of fields
    problem->ReadFields(reader);

    // input of materials of cloned fields (if needed)
    problem->ReadCloningMaterialMap(reader);

    // read all types of geometry related conditions (e.g. boundary conditions)
    // Also read time and space functions and local coord systems
    problem->ReadConditions(reader);

    // read all knot information for isogeometric analysis
    // and add it to the (derived) nurbs discretization
    problem->ReadKnots(reader);
  break;
  case copy_dat_file:
    // group 0 only reads discretization etc
    if (group==0) 
    {
      
      // input of fields
      problem->ReadFields(reader);

      // input of materials of cloned fields (if needed)
      problem->ReadCloningMaterialMap(reader);

      // read all types of geometry related conditions (e.g. boundary conditions)
      // Also read time and space functions and local coord systems
      problem->ReadConditions(reader);

      // read all knot information for isogeometric analysis
      // and add it to the (derived) nurbs discretization
      problem->ReadKnots(reader);
    }
    gcomm->Barrier();
    COMM_UTILS::BroadcastDiscretizations(0); // group 0 broadcasts the discretizations
  break;
  default:
    dserror("nptype (nested parallelity type) not recognized");
    break;
  }
  
  // all reading is done at this point!

  // create control file for output and read restart data if required
  problem->OpenControlFile(*lcomm,
                           inputfile_name,
                           outputfile_kenner,
                           restartfile_kenner);

  // configure the parallel output environment
 {
   const Teuchos::ParameterList& io = problem->IOParams();
   bool screen   = DRT::INPUT::IntegralValue<int>(io,"WRITE_TO_SCREEN");
   bool file     = DRT::INPUT::IntegralValue<int>(io,"WRITE_TO_FILE");
   bool preGrpID = DRT::INPUT::IntegralValue<int>(io,"PREFIX_GROUP_ID");
   int  oproc    = io.get<int>("LIMIT_OUTP_TO_PROC");

   IO::cout.setup(screen, file, preGrpID, lcomm, oproc, group, outputfile_kenner);
 }

  if (lcomm->MyPID()==0)
    problem->WriteInputParameters();

  // before we destroy the reader we want to know about unused sections
  if (npType==copy_dat_file)
  {
    if (group==0) reader.PrintUnknownSections();
  }
  else
    reader.PrintUnknownSections();

  return;
} // end of ntainp_ccadiscret()


