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

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "global_inp_control2.H"

#include "drt_inputreader.H"



/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;


 extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*
  | input of control, element and load information         m.gee 10/06  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret()
{
#ifdef PARALLEL
  int myrank = 0;
  int nproc  = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#else
  Epetra_SerialComm* com = new Epetra_SerialComm();
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#endif

  Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // and now the actual reading
  DRT::DatFileReader reader(allfiles.inputfile_name, comm);
  reader.Activate();

  problem->ReadParameter(reader);

  /* input of not mesh or time based problem data  */
  //inpctr();
  problem->InputControl();

  /* input of materials */
  problem->ReadMaterial();

  /* input of fields */
  problem->ReadFields(reader);

  // read dynamic control data
  if (genprob.timetyp==time_dynamic)
    inpctrdyn();

  // read static control data
  else inpctrstat();

  // read input of eigensolution control data
  inpctreig();

  // read all types of geometry related conditions (e.g. boundary conditions)
  // Also read time and space functions and local coord systems
  Epetra_Time time(*comm);
  if (comm->MyPID()==0)
  {
    cout << "Read conditions                          in....";
    fflush(stdout);
  }
  DRT::Problem::Instance()->ReadConditions();
  if (comm->MyPID()==0)
  {
    cout << time.ElapsedTime() << " secs\n";
  }

#ifdef RESULTTEST
  /*---------------------------------------- input of result descriptions */
  inp_resultdescr();
#endif

  // all reading is done at this point!

  return;
} // end of ntainp_ccadiscret()


#endif  // #ifdef CCADISCRET
