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
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

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
  
  // input of design if desired
  design = NULL;
  //input_design();

  /* input of materials */
  inp_material();
  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();

  /* input of fields */
  inpfield_ccadiscret();

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
  // Note that all discretizations and designs have everything on proc 0 here
  // Note that all discretizations and all design use MPI_COMM_WORLD here
  // These things will be fixed in create_communicators_ccadiscret

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
  RefCountPtr<Epetra_SerialDenseMatrix> tmpnodes = null;
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
  
  // assign nodes to the fields
  int nnode_total = 0;
  for (int i=0; i<genprob.numfld; i++)
  {
    vector<RefCountPtr<DRT::Discretization> >* discretization = 
                 (vector<RefCountPtr<DRT::Discretization> >*)field[i].ccadis;
    for (int j=0;j<field[i].ndis;j++)
    {
      RefCountPtr<DRT::Discretization> actdis = (*discretization)[j];
      input_assign_nodes(*actdis,tmpnodes.get());
      int err = actdis->FillComplete();
      if (err) dserror("Fillcomplete() returned %d",err);
      nnode_total += actdis->NumGlobalNodes();
    }
  }
  // store total number of nodes
  genprob.nnode = nnode_total;

  comm->Barrier(); // everybody wait for proc 0
  return;
} // void inpfield_ccadiscret()


/*-----------------------------------------------------------------------*/
/*!
  \brief sort nodes to the fields

  \author m.gee
  \date   11/06

 */
/*-----------------------------------------------------------------------*/
void input_assign_nodes(DRT::Discretization& actdis,Epetra_SerialDenseMatrix* tmpnodes)
{
  DSTraceHelper dst("input_assign_nodes");
  
  // assign nodes on proc 0 only
  if (actdis.Comm().MyPID()==0)
  {
  
    // allocate a temporary flag array
    vector<int> nodeflag(genprob.nnode);
    for (int i=0; i< genprob.nnode; ++i) nodeflag[i] = -1;
    
    // set flag for each node in this discretization
    for (int i=0; i<actdis.NumMyColElements(); ++i)
    {
      const DRT::Element* actele = actdis.gElement(i);
      const int  nnode = actele->NumNode();
      const int* nodes = actele->NodeIds();
      for (int j=0; j<nnode; ++j)
        nodeflag[nodes[j]] = nodes[j];
    }
    
    // create the nodes and add them to actdis
    for (int i=0; i<genprob.nnode; ++i)
    {
      if (nodeflag[i]==-1) continue;
      double coords[3];
      for (int j=0; j<3; ++j) coords[j] = (*tmpnodes)(i,j);
      RefCountPtr<DRT::Node> node = 
                 rcp(new DRT::Node(nodeflag[i],coords,actdis.Comm().MyPID()));
      actdis.AddNode(node);
    }
  } // if (actdis.Comm().MyPID()==0)
  return;
}


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
  
  structfield->dis = NULL; // not using this here!

  // allocate the discretizations
  vector<RefCountPtr<DRT::Discretization> >* discretization = 
            new vector<RefCountPtr<DRT::Discretization> >(structfield->ndis);
  structfield->ccadis = (void*)discretization;
  for (int i=0; i<structfield->ndis; ++i)
    (*discretization)[i] = rcp(new DRT::Discretization("Structure",comm));

  // count number of elements
  int numele = 0;
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
    numele = counter;
  }

  // read elements (proc 0 only) 
  RefCountPtr<DRT::Discretization> actdis = (*discretization)[0];
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
        RefCountPtr<DRT::Elements::Shell8> ele = 
                         rcp(new DRT::Elements::Shell8(elenumber,actdis->Comm().MyPID()));
        
        // read input for this element
        ele->ReadElement();
        
        // add element to discretization (discretization takes ownership)
        actdis->AddElement(ele);
#endif        
      }

      // elementtyp brick1
      // not yet impl...

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
