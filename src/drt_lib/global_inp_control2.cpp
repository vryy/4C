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
#include <mpi.h>
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
  /* the input of the tracing option has not been done yet, so
     we have to make the dstrc_enter 'by hand'
     */
#ifdef DEBUG
  trace.actroutine = trace.actroutine->next;
  trace.actroutine->name = "ntainp";
  trace.actroutine->dsroutcontrol=TRACEROUT::dsin;
  trace.deepness++;
#endif

  /* input of not mesh or time based problem data  */
  inpctr();

  /* input of materials */
  inp_material();
  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();

  /* input of fields */
  inpfield_ccadiscret_jumbo();

  // read dynamic control data
  if (genprob.timetyp==time_dynamic) inpctrdyn();

  // read static control data
  else inpctrstat();

  // read input of eigensolution control data
  inpctreig();

  // read all types of geometry related conditions (e.g. boundary conditions)
  // Also read time and space functions and local coord systems
  input_conditions();

  /*-------------------------------------------- input of monitoring data */
  inp_monitor();

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
void inpfield_ccadiscret_jumbo()
{
  DSTraceHelper dst("inpfield_ccadiscret_jumbo");

  fflush(stdout);

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
  RefCountPtr<DRT::Discretization> structdis_micro_serial = null;

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

    DRT::NodeReader nodereader(comm, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new DRT::ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Fluid", genprob.numff, comm, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Ale", genprob.numaf, comm, "--ALE ELEMENTS")));

    nodereader.Read();
  }

  else if (genprob.probtyp==prb_ale)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));

    DRT::NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Ale", genprob.numaf, comm, "--ALE ELEMENTS")));
    nodereader.Read();
  }

  else if (genprob.probtyp==prb_fluid)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));

    DRT::NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Fluid", genprob.numff, comm, "--FLUID ELEMENTS")));
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

    DRT::NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.Read();
  } // end of else if (genprob.probtyp==prb_structure)

  else if (genprob.probtyp==prb_struct_multi)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for structural multi-scale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));

    DRT::NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new DRT::ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.Read();

    DRT::NodeReader micronodereader(comm, "--MICROSTRUCTURE NODE COORDS");
    micronodereader.AddElementReader(rcp(new DRT::ElementReader("Micro Structure", genprob.numsf, comm, "--MICROSTRUCTURE ELEMENTS")));
    micronodereader.Read();

    // microscale discretization is distributed over processors but it
    // is needed on every processor redundantly

    vector<RefCountPtr<DRT::Discretization> >* discretization;
    discretization = (vector<RefCountPtr<DRT::Discretization> >*)field[genprob.numsf].ccadis;
    RefCountPtr<DRT::Discretization> structdis_micro = (*discretization)[1];

    RefCountPtr<Epetra_SerialComm> serialcomm = rcp(new Epetra_SerialComm());
    structdis_micro_serial = rcp(new DRT::Discretization("Micro Structure Serial",serialcomm));
    (*discretization)[1] = structdis_micro_serial;

    RefCountPtr<Epetra_Map> parallel_rownodes = rcp(new Epetra_Map(*structdis_micro->NodeRowMap()));
    RefCountPtr<Epetra_Map> parallel_roweles  = rcp(new Epetra_Map(*structdis_micro->ElementRowMap()));

    // build redundant colnodes
    vector<int> mygid(parallel_rownodes->NumMyElements());
    for (int i=0; i<parallel_rownodes->NumMyElements(); ++i) mygid[i] = parallel_rownodes->MyGlobalElements()[i];
    vector<int> rmygid(0);
    vector<int> targetprocs(structdis_micro->Comm().NumProc());
    for (int i=0; i<structdis_micro->Comm().NumProc(); ++i) targetprocs[i] = i;
    LINALG::Gather(mygid,rmygid,structdis_micro->Comm().NumProc(),&targetprocs[0], structdis_micro->Comm());
    RefCountPtr<Epetra_Map> redundantmap = rcp(new Epetra_Map(-1,(int)rmygid.size(),&rmygid[0],0,structdis_micro->Comm()));

    // build redundant coleles
    vector<int> mygid_ele(parallel_roweles->NumMyElements());
    for (int i=0; i<parallel_roweles->NumMyElements(); ++i) mygid_ele[i] = parallel_roweles->MyGlobalElements()[i];
    vector<int> rmygid_ele(0);
    LINALG::Gather(mygid_ele,rmygid_ele,structdis_micro->Comm().NumProc(),&targetprocs[0], structdis_micro->Comm());
    RefCountPtr<Epetra_Map> redundantmap_ele = rcp(new Epetra_Map(-1,(int)rmygid_ele.size(),&rmygid_ele[0],0,structdis_micro->Comm()));

    structdis_micro->ExportColumnNodes(*redundantmap);
    structdis_micro->ExportColumnElements(*redundantmap_ele);
    int err = structdis_micro->FillComplete();
    if (err) dserror("structdis_micro->FillComplete() returned %d",err);

    for (int i=0; i<structdis_micro->NumMyColElements(); ++i)
    {
      DRT::Element* actele = structdis_micro->lColElement(i);
      RefCountPtr<DRT::Element> newele = rcp(actele->Clone());
      newele->SetOwner(0);
      structdis_micro_serial->AddElement(newele);
    }
    for (int i=0; i<structdis_micro->NumMyColNodes(); ++i)
    {
      DRT::Node* actnode = structdis_micro->lColNode(i);
      RefCountPtr<DRT::Node> newnode = rcp(actnode->Clone());
      newnode->SetOwner(0);
      structdis_micro_serial->AddNode(newnode);
    }
    err = structdis_micro_serial->FillComplete();

    if (1)
    {
    for (int i=0; i<structdis_micro->Comm().NumProc(); ++i)
    {
      cout << "\n";
      if (structdis_micro->Comm().MyPID()==i)
      {
        cout << "Proc " << i << "meine serielle Diskretisierung:\n";
        cout << *structdis_micro_serial;
      }
      structdis_micro->Comm().Barrier();
    }
    cout << "\n";
    }
  } // end of else if (genprob.probtyp==prb_struct_multi)

  else dserror("Type of problem unknown");

  return;
} // void inpfield_ccadiscret_jumbo()



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
