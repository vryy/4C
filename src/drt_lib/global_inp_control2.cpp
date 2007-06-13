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
#ifdef TRILINOS_PACKAGE

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_SerialDenseMatrix.h"
#include "global_inp_control2.H"
#include "Epetra_SerialComm.h"

#include "drt_inputreader.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

#ifdef DEBUG
/*!----------------------------------------------------------------------
  \brief the tracing variable

  <pre>                                                         m.gee 8/00
  defined in pss_ds.c, declared in tracing.h
  </pre>
 *----------------------------------------------------------------------*/
extern struct _CCA_TRACE         trace;
#endif

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

  // and now the actual reading
  DRT::DatFileReader reader(allfiles.inputfile_name, comm);
  reader.Activate();

  /* input of not mesh or time based problem data  */
  inpctr();

  /* input of materials */
  input_material_ccadiscret(*DRT::Problem::Instance());

  /* input of fields */
  inpfield_ccadiscret(*DRT::Problem::Instance(), reader);

  // read dynamic control data
  if (genprob.timetyp==time_dynamic) inpctrdyn();

  // read static control data
  else inpctrstat();

  // read input of eigensolution control data
  inpctreig();

  // read all types of geometry related conditions (e.g. boundary conditions)
  // Also read time and space functions and local coord systems
  input_conditions(*DRT::Problem::Instance());

#ifdef RESULTTEST
  /*---------------------------------------- input of result descriptions */
  inp_resultdescr();
#endif

  // all reading is done at this point!

  return;
} // end of ntainp_ccadiscret()


/*----------------------------------------------------------------------*
  | input of fields                                        m.gee 03/07  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void inpfield_ccadiscret(DRT::Problem& problem, DRT::DatFileReader& reader)
{
  fflush(stdout);

  genprob.create_dis = 0;
  genprob.create_ale = 0;
  genprob.maxnode    = 0;
  genprob.nodeshift  = genprob.nnode;

  // read elements the first time to create graph object
  // row distribution of nodes
  // column distribution of nodes
  // graph of problem
  RefCountPtr<Epetra_Map> rownodes   = null;
  RefCountPtr<Epetra_Map> colnodes   = null;
  RefCountPtr<Epetra_Map> roweles    = null;
  RefCountPtr<Epetra_Map> coleles    = null;
  RefCountPtr<Epetra_CrsGraph> graph = null;

  RefCountPtr<DRT::Discretization> structdis = null;
  RefCountPtr<DRT::Discretization> fluiddis  = null;
  RefCountPtr<DRT::Discretization> aledis    = null;
  RefCountPtr<DRT::Discretization> structdis_macro = null;
  RefCountPtr<DRT::Discretization> structdis_micro = null;


  if (genprob.probtyp == prb_fsi)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=3) dserror("numfld != 3 for fsi problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));
    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));
    aledis = rcp(new DRT::Discretization("Ale",reader.Comm()));

    problem.AddDis(genprob.numsf, structdis);
    problem.AddDis(genprob.numff, fluiddis);
    problem.AddDis(genprob.numaf, aledis);

    DRT::NodeReader nodereader(reader, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::ElementReader(aledis, reader, "--ALE ELEMENTS")));

    nodereader.Read();
  }

  else if (genprob.probtyp==prb_ale)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));

    aledis = rcp(new DRT::Discretization("Ale",reader.Comm()));
    problem.AddDis(genprob.numaf, aledis);

    DRT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader(aledis, reader, "--ALE ELEMENTS")));
    nodereader.Read();
  }

  else if (genprob.probtyp==prb_fluid  || genprob.probtyp==prb_condif)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));

    fluiddis = rcp(new DRT::Discretization("Fluid",reader.Comm()));
    problem.AddDis(genprob.numff, fluiddis);

    DRT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS")));
    nodereader.Read();
  }

  else if (genprob.probtyp==prb_fluid_pm)
    dserror("prb_fluid_pm not yet impl.");

  else if (genprob.probtyp == prb_tsi)
    dserror("prb_tsi not yet impl.");

  else if (genprob.probtyp==prb_structure)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));

    structdis = rcp(new DRT::Discretization("Structure",reader.Comm()));
    problem.AddDis(genprob.numsf, structdis);

    DRT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS")));
    nodereader.Read();
  } // end of else if (genprob.probtyp==prb_structure)

  else if (genprob.probtyp==prb_struct_multi)
  {
    // allocate and input general old stuff....

    if (genprob.numfld!=1) dserror("numfld != 1 for structural multi-scale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));


    // read macroscale fields from main inputfile

    structdis_macro = rcp(new DRT::Discretization("Macro Structure",reader.Comm()));
    problem.AddDis(genprob.numsf, structdis_macro);

    DRT::NodeReader nodereader(reader, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader(structdis_macro, reader, "--STRUCTURE ELEMENTS")));
    nodereader.Read();


    // read microscale fields from second inputfile

    RefCountPtr<DRT::Problem> micro_problem = DRT::Problem::Instance(1);
    RefCountPtr<Epetra_SerialComm> serialcomm = rcp(new Epetra_SerialComm());

    string micro_inputfile_name = inp_micro_filename(structdis_macro->Comm().MyPID());

    DRT::DatFileReader micro_reader(micro_inputfile_name, serialcomm, 1);
    micro_reader.Activate();

    structdis_micro = rcp(new DRT::Discretization("Micro Structure", micro_reader.Comm()));
    micro_problem->AddDis(genprob.numsf, structdis_micro);

    DRT::NodeReader micronodereader(micro_reader, "--NODE COORDS");
    micronodereader.AddElementReader(rcp(new DRT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS")));
    micronodereader.Read();


    // read materials of microscale

    input_material_ccadiscret(*micro_problem);


    // read conditions of microscale -> note that no time curves and
    // spatial functions can be read!

    input_conditions(*micro_problem);


    // reactivate reader of macroscale as well as macroscale material

    reader.Activate();
    problem.ActivateMaterial();

    if (0)
    {
    for (int i=0; i<structdis_macro->Comm().NumProc(); ++i)
    {
      cout << "\n";
      if (structdis_macro->Comm().MyPID()==i)
      {
        cout << "Proc " << i << "meine serielle Diskretisierung:\n";
        cout << *structdis_micro;
      }
      structdis_macro->Comm().Barrier();
    }
    cout << "\n";
    }
    exit(0);
  } // end of else if (genprob.probtyp==prb_struct_multi)

  else dserror("Type of problem unknown");

  return;
} // void inpfield_ccadiscret()


string inp_micro_filename(int mypid)
{
  string microfile = "MICROFILE";
  string tmp;
  string micro_inputfile_name;
  ifstream file(allfiles.inputfile_name);
  int filecount=0;
  int found=0;

  for (; file; ++filecount)
  {
    file >> tmp;

    if (tmp == microfile)
    {
      file >> micro_inputfile_name;
      found = 1;
    }
  }

  if (found == 0)
    dserror("No inputfile for microstructure given!\n");

  if (!mypid)
    cout << "input for microscale is read from        " << micro_inputfile_name << "\n";

  return micro_inputfile_name;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
