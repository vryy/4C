/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <string>

#include "global_inp_control2.H"
#include "drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "drt_inputreader.H"
#include "standardtypes_cpp.H"


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
  | input of control, element and load information         m.gee 10/06  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret(
  std::string& inputfile_name,
  std::string& outputfile_kenner,
  std::string& restartfile_kenner
  )
{

  Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();

  // create error files
  // call this one rather early, since ReadConditions etc
  // underlying methods may try to write to allfiles.out_err
  // this old-style global variable is set as well
  problem->OpenErrorFile(*lcomm,outputfile_kenner);

  // and now the actual reading
  DRT::INPUT::DatFileReader reader(inputfile_name,
                                   lcomm);

  problem->ReadParameter(reader);

  /* input of not mesh or time based problem data  */
  problem->InputControl();

  /* input of materials */
  problem->ReadMaterials(reader);

  /* input of fields */
  problem->ReadFields(reader);

  /* input of materials of cloned fields (if needed) */
  problem->ReadClonedMaterials(reader);

  // read all types of geometry related conditions (e.g. boundary conditions)
  // Also read time and space functions and local coord systems
  problem->ReadConditions(reader);

  // read all knot information for isogeometric analysis
  // and add it to the (derived) nurbs discretization
  problem->ReadKnots(reader);

  // all reading is done at this point!

  // create control file for output and read restart data if required
  problem->OpenControlFile(*lcomm,
                           inputfile_name,
                           outputfile_kenner,
                           restartfile_kenner);

  if (lcomm->MyPID()==0)
    problem->WriteInputParameters();

  // before we destroy the reader we want to know about unused sections
  reader.PrintUnknownSections();

  return;
} // end of ntainp_ccadiscret()


#endif  // #ifdef CCADISCRET
