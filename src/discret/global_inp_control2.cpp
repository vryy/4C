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

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Epetra_SerialDenseMatrix.h"

#include "global_inp_control2.H"



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
#ifdef PERF
  perf_begin(3);
#endif
  inpctr();
#ifdef PERF
  perf_end(3);
#endif

  /* input of materials */
#ifdef PERF
  perf_begin(6);
#endif
  inp_material();
  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();
#ifdef PERF
  perf_end(6);
#endif

  /* input of fields */
#ifdef PERF
  perf_begin(7);
#endif
  inpfield_ccadiscret();
#ifdef PERF
  perf_end(7);
#endif


cout << "I'm here" << endl;
exit(0);


  return;
} // end of ntainp_ccadiscret()


/*----------------------------------------------------------------------*
  | input of fields                                        m.gee 10/06  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void inpfield_ccadiscret()
{
  DSTraceHelper dst("inpfield_ccadiscret");
  
  int myrank = 0;
  int nproc  = 1;

#ifdef PARALLEL
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

  // create the discretization on proc 0 only
  // later on we'll use metis to partition the whole thing

  // read nodal coords in temporary array (proc 0 only)
  // allocate temporary array for nodal coords
  RefCountPtr<Epetra_SerialDenseMatrix> tmpnodes;
  if (myrank==0)
  {
    tmpnodes = rcp(new Epetra_SerialDenseMatrix(genprob.nnode,3));
    // read nodal coords
    inpnodes_ccadiscret(*tmpnodes);
  }
  
  // read elements
  if (genprob.probtyp == prb_fsi)
    dserror("prb_fsi not yet impl.");
    
  if (genprob.probtyp==prb_fluid)
    dserror("prb_fluid not yet impl.");

  if (genprob.probtyp==prb_fluid_pm)
    dserror("prb_fluid_pm not yet impl.");
    
  if (genprob.probtyp == prb_tsi)
    dserror("prb_tsi not yet impl.");
    
  if (genprob.probtyp==prb_structure)
  {
    if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    input_structural_field(&(field[genprob.numsf]),comm);
  }




  comm->Barrier(); // everybody wait for proc 0
  return;
} // void inpfield_ccadiscret()


/*-----------------------------------------------------------------------*/
/*!
  \brief input of structure field

  Create the structure field: allocate the discretizations, the required
  number of elements and then read and create the elements

  \param structfield    FIELD  (i) pointer to the structure field

  \return void

  \author m.gee
  \date   11/06

 */
/*-----------------------------------------------------------------------*/
void input_structural_field(FIELD *structfield, RefCountPtr<Epetra_Comm> comm)
{
  DSTraceHelper dst("input_structural_field");
  
  /* allocate discretizations */
  //structfield->dis = (DISCRET*)CCACALLOC(structfield->ndis,sizeof(DISCRET));
  structfield->dis = NULL; // not using this here!
  structfield->dis[0].disclass = dc_normal;
  /* initialize array positions with -1 */
  memset(&structfield->dis[0].ipos, 0xff, sizeof(ARRAY_POSITION));

  // allocate the discretizations
  vector<RefCountPtr<CCADISCRETIZATION::Discretization> >* discretization = 
    new vector<RefCountPtr<CCADISCRETIZATION::Discretization> >(structfield->ndis);
  structfield->ccadis = (void*)discretization;
  for (int i=0; i<structfield->ndis; ++i)
    (*discretization)[i] = rcp(new CCADISCRETIZATION::Discretization(comm));

  // count number of elements
  {
    int counter=0;
    if (frfind("--STRUCTURE ELEMENTS")==1)
    {
      frread();
      while(strncmp(allfiles.actplace,"------",6)!=0)
      {
        counter++;
        frread();
      }
    }
    structfield->dis[0].numele = counter;
  }

  // read elements  
  RefCountPtr<CCADISCRETIZATION::Discretization> actdis = (*discretization)[0];
  if (actdis->Comm().MyPID()==0)
  {
    if (frfind("--STRUCTURE ELEMENTS")==0) return;
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      char *colpointer = allfiles.actplace;
      int elenumber    = strtol(colpointer,&colpointer,10);
      --elenumber;
      int ierr=0;
      // elementtyp shell8
      frchk("SHELL8",&ierr);
      if (ierr==1)
      {
#ifndef D_SHELL8
        dserror("SHELL8 needed but not defined in Makefile");
#else  
        RefCountPtr<CCADISCRETIZATION::Shell8> shell8 = 
                              rcp(new CCADISCRETIZATION::Shell8(elenumber));
        shell8->ReadElement();
#endif        
      }


      // elementtyp brick1
      
      
      

      frread();
    } // while(strncmp(allfiles.actplace,"------",6)!=0)  
  } // if (actdis->Comm().MyPID()==0)

  frrewind();
  return;
} // void input_structural_field

/*----------------------------------------------------------------------*
  | input of nodal coords (proc 0 only)                    m.gee 10/06  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void inpnodes_ccadiscret(Epetra_SerialDenseMatrix& tmpnodes)
{
  DSTraceHelper dst("inpnodes_ccadiscret");

  if (frfind("--NODE COORDS")==0) dserror("frfind: NODE COORDS is not in input file");
  frread();
  int counter=0;
  int nodeid=0;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    int ierr;
    frint("NODE",&(nodeid),&ierr);
    if (ierr!=1) dserror("reading of nodes failed");
    if (nodeid-1 != counter) dserror("Reading of nodes failed: Nodes must be numbered consecutive!!");
    double nodes[3];
    frdouble_n("COORD",nodes,3,&ierr);
    if (ierr!=1) dserror("reading of nodes failed");
    for (int i=0; i<3; ++i) tmpnodes(counter,i) = nodes[i];
    counter++;
    frread();
  }
  frrewind();
  return;
} // void inpnodes_ccadiscret


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
