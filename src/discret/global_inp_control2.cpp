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
  input_design();

  /* input of materials */
  inp_material();
  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();

  /* input of fields */
  inpfield_ccadiscret();

  // read the design topology
  input_design_topology_discretization();


  cout << "Reached regular exit\n"; fflush(stdout);
  exit(0);


  return;
} // end of ntainp_ccadiscret()


/*-----------------------------------------------------------------------*/
/*!
  \brief read design 

  \author m.gee
  \date   11/06

 */
/*-----------------------------------------------------------------------*/
void input_design()
{
  DSTraceHelper dst("input_design");
  

  // if there's no design description, do nothing
  if (frfind("--DESIGN DESCRIPTION")==0)
  {
    design = NULL;
    return;
  }
  // allocate the design structure
  design = (DESIGN*)CCACALLOC(1,sizeof(DESIGN));
#ifdef PARALLEL
  Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#else
  Epetra_SerialComm* com = new Epetra_SerialComm();
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#endif
  //------- allocate 3 discretizations for lines, surfaces and volumes
  vector<RefCountPtr<CCADISCRETIZATION::Design> >* dptr = 
                 new vector<RefCountPtr<CCADISCRETIZATION::Design> >(3);  
  for (int i=0; i<3; ++i)
    (*dptr)[i] = rcp(new CCADISCRETIZATION::Design(comm));
  design->ccadesign = (void*)dptr;
  RefCountPtr<CCADISCRETIZATION::Design> designlines = (*dptr)[0];
  RefCountPtr<CCADISCRETIZATION::Design> designsurfs = (*dptr)[1];
  RefCountPtr<CCADISCRETIZATION::Design> designvols  = (*dptr)[2];

  int ierr=0;
  
  //------------------------------------------------ read design sizes
  frread();
  int numdnode=0;
  frint("NDPOINT",&numdnode,&ierr);
  if (!ierr) dserror("Cannot read design");
  frread();
  
  int numdline=0;
  frint("NDLINE",&numdline,&ierr);
  if (!ierr) dserror("Cannot read design");
  frread();
  
  int numdsurf=0;
  frint("NDSURF",&numdsurf,&ierr);
  if (!ierr) dserror("Cannot read design");
  frread();
      
  int numdvol=0;
  frint("NDVOL",&numdvol,&ierr);
  if (!ierr) dserror("Cannot read design");
  
  design->ndnode = numdnode;
  design->ndline = numdline;
  design->ndsurf = numdsurf;
  design->ndvol  = numdvol;

  frrewind();

  // create whole design on proc 0 only
  if (comm->MyPID()==0)
  {
    //--------------------------------------------- input of design nodes
    if (frfind("--DESIGN POINTS")==0)
      dserror("frfind: DESIGN POINTS is not in input file");
    frread();
    for (int i=0; i<numdnode; ++i)
    {
      ierr=0;
      int readid = i+1;
      frchk("POINT",&ierr);
      while (!ierr) {frread(); frchk("POINT",&ierr);}
      
      int num=0;
      frchk("Num:",&ierr);
      while (!ierr) {frread(); frchk("Num:",&ierr);}
      frint("Num:",&num,&ierr);
      if (!ierr) dserror("Cannot read DNODE");
      if (num != readid) dserror("DNODEs got mixed up");
      
      int ncond=0;
      frint("conditions:",&ncond,&ierr);
      if (!ierr) dserror("Cannot read DNODE");
      
      frchk("Coord:",&ierr);
      while (!ierr) {frread(); frchk("Coord:",&ierr);}
      double x[3];
      frdouble_n("Coord:",x,3,&ierr);
      if (!ierr) dserror("Cannot read DNODE");
      
      frchk("END POINT",&ierr);
      while (!ierr) {frread(); frchk("END POINT",&ierr);}
      frread();
    
      // create the design node and store it in ccadesign
      RefCountPtr<CCADISCRETIZATION::DesignNode> node = 
                      rcp(new CCADISCRETIZATION::DesignNode(i,x));
      
      // we add the nodes to all of the design descriptions
      // this is not a problem as they are refcountpointed
      designlines->AddNode(node);                      
    } // for (int i=0; i<numdnode; ++i)

    //----------------------------------------------- input of design lines
    if (frfind("--DESIGN LINES")==0)
      dserror("frfind: DESIGN LINES is not in input file");
    frread();
    for (int i=0; i<numdline; ++i)
    {
      ierr=0;
      int readid = i+1;
      frchk("LINE",&ierr);
      while (!ierr) {frread(); frchk("LINE",&ierr);}
      
      frchk("Num:",&ierr);
      while (!ierr) {frread(); frchk("Num:",&ierr);}
      int num=0;
      frint("Num:",&num,&ierr);
      if (!ierr) dserror("Cannot read DLINE");
      if (num != readid) dserror("DLINEs got mixed up");
      
      // this is currently only reading 2 nodes
      int nodeids[2];
      frchk("Points:",&ierr);
      while (!ierr) {frread(); frchk("Points:",&ierr);}
      frint_n("Points:",nodeids,2,&ierr);
      if (!ierr) dserror("Cannot read DLINE");
      nodeids[0]--;
      nodeids[1]--;
      
      frchk("END",&ierr);
      while (!ierr) {frread(); frchk("END",&ierr);}
      frread();
      
      // create the design line element and store it in ccadesign
      RefCountPtr<CCADISCRETIZATION::DesignElement> line = 
        rcp(new CCADISCRETIZATION::DesignElement(i,CCADISCRETIZATION::Element::element_designline)); 
      line->SetNodeIds(2,nodeids);     
      
      // Add the line to the lines description
      designlines->AddElement(line);
    } // for (int i=0; i<numdline; ++i)
    
    //------------------------------------------ input of design surfaces
    if (frfind("--DESIGN SURFACES")==0)
      dserror("frfind: DESIGN SURFACES is not in input file");
    frread();
    for (int i=0; i<numdsurf; ++i)
    {
      ierr=0;
      int readid = i+1;
      frchk("NURBSURFACE",&ierr);
      while (!ierr) { frread(); frchk("NURBSURFACE",&ierr); }
      
      // read surface number
      int num=0;
      frchk("Num:",&ierr);
      while (!ierr) {frread(); frchk("Num:",&ierr);}
      frint("Num:",&num,&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      if (num != readid) dserror("DSURF got mixed up");
      
      // read number of my lines
      frchk("NumLines:",&ierr);
      while (!ierr) {frread(); frchk("NumLines:",&ierr);}
      int numlines=0;
      frint("NumLines:",&numlines,&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      
      vector<int> lineids(numlines);
      vector<int> orientation(numlines);
      
      // read lines
      frchk("Line:",&ierr);
      while (!ierr) {frread(); frchk("Line:",&ierr);}
      for (int j=0; j<numlines; ++j)
      {
        frint("Line:",&(lineids[j]),&ierr);
        if (!ierr) dserror("Cannot read DSURF");
        lineids[j]--;
        char buffer[100];
        frchar("Orientation:",buffer,&ierr);
        if (!ierr) dserror("Cannot read DSURF");
        if (strncmp("SAME1ST",buffer,7)==0) orientation[j]=0;
        else                                orientation[j]=1;
        frread();
      }
      
      // create the design surface
      RefCountPtr<CCADISCRETIZATION::DesignElement> surf = 
        rcp(new CCADISCRETIZATION::DesignElement(i,CCADISCRETIZATION::Element::element_designsurface));
      surf->SetLowerEntities(numlines,&lineids[0],&orientation[0]);
      
      // Add the surface to the surfaces description
      designsurfs->AddElement(surf);
    } // for (int i=0; i<numdsurf; ++i)
    
    //-------------------------------------------------- input of design volumes
    if (frfind("--DESIGN VOLUMES")==0)
      dserror("frfind: DESIGN VOLUMES is not in input file");
    frread();
    for (int i=0; i<numdvol; ++i)
    {
      ierr=0;
      int readid = i+1;
      frchk("VOLUME",&ierr);
      while (!ierr) {frread(); frchk("VOLUME",&ierr);}
      int num=0;
      frint("Num:",&num,&ierr);
      if (!ierr) dserror("Cannot read DVOL");
      if (num != readid) dserror("DVOLs got mixed up");
      
      frchk("NumSurfaces:",&ierr);
      while (!ierr) {frread(); frchk("NumSurfaces:",&ierr);}
      int numsurfs;
      frint("NumSurfaces:",&numsurfs,&ierr);
      if (!ierr) dserror("Cannot read DVOL");
      
      vector<int> surfids(numsurfs);
      vector<int> orientation(numsurfs);
      
      // read surfs adjacent to this volume
      frchk("Surface:",&ierr);
      while (!ierr) {frread(); frchk("Surface:",&ierr);}
      for (int j=0; j<numsurfs; ++j)
      {
        frint("Surface:",&(surfids[j]),&ierr);
        if (!ierr) dserror("Cannot read DVOL");
        surfids[j]--;
        char buffer[100];
        frchar("Orientation:",buffer,&ierr);
        if (!ierr) dserror("Cannot read DVOL");
        if (strncmp("SAME1ST",buffer,7)==0) orientation[j]=0;
        else                                orientation[j]=1;
        frread();
      }
      frchk("END VOLUME",&ierr);
      while (!ierr) {frread(); frchk("END VOLUME",&ierr);}
      frread();
      
      // create the design volume
      RefCountPtr<CCADISCRETIZATION::DesignElement> vol = 
        rcp(new CCADISCRETIZATION::DesignElement(i,CCADISCRETIZATION::Element::element_designsurface));
      vol->SetLowerEntities(numsurfs,&surfids[0],&orientation[0]);
      
      // Add the volume to the volumes description
      designvols->AddElement(vol);
    } // for (int i=0; i<numdvol; ++i)

  } // if (comm->MyPID()==0)

  // all design line, surface and volume discretizations have been read 
  // at this point by proc 0
  // Now, call FillComplete() on them by all procs
  // We call a special Fillcomplete on volumes, surfaces and surfaces
  // that create pointers to the other discretizations:
  // volumes <-> surfaces
  // surfaces <-> lines
  // lines <-> nodes (standard Discretization::Fillcomplete()) 
  ierr=0;
  ierr += designvols->FillComplete(NULL,designsurfs.get());
  ierr += designsurfs->FillComplete(designvols.get(),designlines.get());
  ierr += designlines->FillComplete(designsurfs.get(),NULL);
  if (ierr)
    dserror("FillComplete of Design returned %d",ierr);
  
  return;
}

/*-----------------------------------------------------------------------*/
/*!
  \brief read design <-> discretization topology

  \author m.gee
  \date   11/06

 */
/*-----------------------------------------------------------------------*/
void input_design_topology_discretization()
{
  DSTraceHelper dst("input_design_topology_discretization");
  






  return;
}

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
    vector<RefCountPtr<CCADISCRETIZATION::Discretization> >* discretization = 
      (vector<RefCountPtr<CCADISCRETIZATION::Discretization> >*)field[i].ccadis;
    for (int j=0;j<field[i].ndis;j++)
    {
      RefCountPtr<CCADISCRETIZATION::Discretization> actdis = (*discretization)[i];
      input_assign_nodes(*actdis,*tmpnodes);
      nnode_total += actdis->NumGlobalNodes();
      int err = actdis->FillComplete();
      if (err)
        dserror("Fillcomplete() returned %d",err);
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
void input_assign_nodes(CCADISCRETIZATION::Discretization& actdis,Epetra_SerialDenseMatrix& tmpnodes)
{
  DSTraceHelper dst("input_assign_nodes");
  
  // assign nodes on proc 0 only
  if (actdis.Comm().MyPID()==0)
  {
  
    // allocate a temporary flag array
    vector<int> nodeflag(genprob.nnode);
    for (int i=0; i< genprob.nnode; ++i) nodeflag[i] = -1;
    
    // set flag for each node in this discretization
    for (int i=0; i<actdis.NumMyElements(); ++i)
    {
      const CCADISCRETIZATION::Element* actele = actdis.gElement(i);
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
      for (int j=0; j<3; ++j) coords[j] = tmpnodes(i,j);
      RefCountPtr<CCADISCRETIZATION::Node> node = 
        rcp(new CCADISCRETIZATION::Node(nodeflag[i],coords));
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
  vector<RefCountPtr<CCADISCRETIZATION::Discretization> >* discretization = 
    new vector<RefCountPtr<CCADISCRETIZATION::Discretization> >(structfield->ndis);
  structfield->ccadis = (void*)discretization;
  for (int i=0; i<structfield->ndis; ++i)
    (*discretization)[i] = rcp(new CCADISCRETIZATION::Discretization(comm));

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
        RefCountPtr<CCADISCRETIZATION::Shell8> ele = 
                              rcp(new CCADISCRETIZATION::Shell8(elenumber));
        
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
