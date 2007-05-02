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



static void input_point_dirich(multimap<int,RefCountPtr<DRT::Condition> >& pdmap);
static void input_line_dirich(multimap<int,RefCountPtr<DRT::Condition> >& ldmap);
static void input_surf_dirich(multimap<int,RefCountPtr<DRT::Condition> >& sdmap);
static void input_vol_dirich(multimap<int,RefCountPtr<DRT::Condition> >& vdmap);

static void input_point_neum(multimap<int,RefCountPtr<DRT::Condition> >& pnmap);
static void input_line_neum(multimap<int,RefCountPtr<DRT::Condition> >& lnmap);
static void input_surf_neum(multimap<int,RefCountPtr<DRT::Condition> >& snmap);
static void input_vol_neum(multimap<int,RefCountPtr<DRT::Condition> >& vnmap);

static void add_nodeids_to_condition(const int id, RefCountPtr<DRT::Condition> cond,
                                     const vector<int> nd_fenode,
                                     const vector<vector<int> > d_fenode);

/*----------------------------------------------------------------------*
 | input of conditions                                    m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_conditions()
{
  DSTraceHelper dst("input_conditions");
  /*---------------------------------------------- input of time curves */
  inp_cond_curve();
  /*---------------------------------------- input of spatial functions */
  inp_cond_funct();
  //------------------------------- read number of design objects we have
  // this currently serves to determine how many node sets we might have
  int ndnode=0;
  int ndline=0;
  int ndsurf=0;
  int ndvol=0;
  int ierr=0;
  frrewind();
  if (frfind("--DESIGN DESCRIPTION")==0)
    dserror("Cannot find design description");
  frread();
  frint("NDPOINT",&ndnode,&ierr); if (!ierr) dserror("Cannot read design");
  frread();
  frint("NDLINE",&ndline,&ierr);  if (!ierr) dserror("Cannot read design");
  frread();
  frint("NDSURF",&ndsurf,&ierr);  if (!ierr) dserror("Cannot read design");
  frread();
  frint("NDVOL",&ndvol,&ierr);    if (!ierr) dserror("Cannot read design");
  frrewind();
  //------------------------------- myrank and nproc from discretization
  vector<RefCountPtr<DRT::Discretization> >* discretization =
    (vector<RefCountPtr<DRT::Discretization> >*)field[0].ccadis;
  RefCountPtr<DRT::Discretization> actdis = (*discretization)[0];
  //--------------------------------------------- read generic node sets
  // read design nodes <-> nodes
  vector<int> ndnode_fenode(ndnode);
  vector<vector<int> > dnode_fenode(ndnode);
  for (int i=0; i<ndnode; ++i)
    ndnode_fenode[i] = 0;
  input_design_dpoint_fenode_read(dnode_fenode,ndnode_fenode);

  // read design lines <-> nodes
  vector<int> ndline_fenode(ndline);
  vector<vector<int> > dline_fenode(ndline);
  for (int i=0; i<ndline; ++i)
    ndline_fenode[i] = 0;
  input_design_dline_fenode_read(dline_fenode,ndline_fenode);

  // read design surfaces <-> nodes
  vector<int> ndsurf_fenode(ndsurf);
  vector<vector<int> > dsurf_fenode(ndsurf);
  for (int i=0; i<ndsurf; ++i)
    ndsurf_fenode[i] = 0;
  input_design_dsurf_fenode_read(dsurf_fenode,ndsurf_fenode);

  // read design volumes <-> nodes
  vector<int> ndvol_fenode(ndvol);
  vector<vector<int> > dvol_fenode(ndvol);
  for (int i=0; i<ndvol; ++i)
    ndvol_fenode[i] = 0;
  input_design_dvol_fenode_read(dvol_fenode,ndvol_fenode);

  //-------------------------------------read point dirichlet conditions
  multimap<int,RefCountPtr<DRT::Condition> > pointdirich;
  input_point_dirich(pointdirich);
  //-------------------------------------read line dirichlet conditions
  multimap<int,RefCountPtr<DRT::Condition> > linedirich;
  input_line_dirich(linedirich);
  //-------------------------------------read surface dirichlet conditions
  multimap<int,RefCountPtr<DRT::Condition> > surfdirich;
  input_surf_dirich(surfdirich);
  //-------------------------------------read volume dirichlet conditions
  multimap<int,RefCountPtr<DRT::Condition> > voldirich;
  input_vol_dirich(voldirich);
  //--------------------------------------- read point neumann conditions
  multimap<int,RefCountPtr<DRT::Condition> > pointneum;
  input_point_neum(pointneum);
  //--------------------------------------- read line neumann conditions
  multimap<int,RefCountPtr<DRT::Condition> > lineneum;
  input_line_neum(lineneum);
  //--------------------------------------- read surface neumann conditions
  multimap<int,RefCountPtr<DRT::Condition> > surfneum;
  input_surf_neum(surfneum);
  //--------------------------------------- read vol neumann conditions
  multimap<int,RefCountPtr<DRT::Condition> > volneum;
  input_vol_neum(volneum);

  // Iterate through all conditions and add the finite element node ids
  // to the condition itself
  multimap<int,RefCountPtr<DRT::Condition> >::iterator curr;
  // iterate through point dirichlet conditions and add fe nodes
  for (curr=pointdirich.begin(); curr!=pointdirich.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndnode_fenode,dnode_fenode);
  // iterate through line dirichlet conditions and add fe nodes
  for (curr=linedirich.begin(); curr!=linedirich.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndline_fenode,dline_fenode);
  // iterate through surface dirichlet conditions and add fe nodes
  for (curr=surfdirich.begin(); curr!=surfdirich.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndsurf_fenode,dsurf_fenode);
  // iterate through surface dirichlet conditions and add fe nodes
  for (curr=voldirich.begin(); curr!=voldirich.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndvol_fenode,dvol_fenode);

  // iterate through point neumann conditions and add fe nodes
  for (curr=pointneum.begin(); curr!=pointneum.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndnode_fenode,dnode_fenode);
  // iterate through line neumann conditions and add fe nodes
  for (curr=lineneum.begin(); curr!=lineneum.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndline_fenode,dline_fenode);
  // iterate through surface neumann conditions and add fe nodes
  for (curr=surfneum.begin(); curr!=surfneum.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndsurf_fenode,dsurf_fenode);
  // iterate through surface neumann conditions and add fe nodes
  for (curr=volneum.begin(); curr!=volneum.end(); ++curr)
    add_nodeids_to_condition(curr->first,curr->second,ndvol_fenode,dvol_fenode);

  // Iterate through all discretizations and sort the appropiate condition into
  // the correct discretization it applies to
  for (int i=0; i<genprob.numfld; i++)
  {
    vector<RefCountPtr<DRT::Discretization> >* discretization =
              (vector<RefCountPtr<DRT::Discretization> >*)field[i].ccadis;
    for (int j=0;j<field[i].ndis;j++)
    {
      RefCountPtr<DRT::Discretization> actdis = (*discretization)[j];
      const Epetra_Map* noderowmap = actdis->NodeRowMap();

      // point dirichlet
      for (curr=pointdirich.begin(); curr!=pointdirich.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("Dirichlet",curr->second);
      }
      // line dirichlet
      for (curr=linedirich.begin(); curr!=linedirich.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("Dirichlet",curr->second);
      }
      // surface dirichlet
      for (curr=surfdirich.begin(); curr!=surfdirich.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("Dirichlet",curr->second);
      }
      // volume dirichlet
      for (curr=voldirich.begin(); curr!=voldirich.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("Dirichlet",curr->second);
      }

      // point neumann
      for (curr=pointneum.begin(); curr!=pointneum.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("PointNeumann",curr->second);
      }
      // line neumann
      for (curr=lineneum.begin(); curr!=lineneum.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("LineNeumann",curr->second);
      }
      // surface neumann
      for (curr=surfneum.begin(); curr!=surfneum.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("SurfaceNeumann",curr->second);
      }
      // volume neumann
      for (curr=volneum.begin(); curr!=volneum.end(); ++curr)
      {
        const vector<int>* nodes = curr->second->Get<vector<int> >("Node Ids");
        if (!(int)nodes->size()) dserror("Condition has no nodal cloud");
        const int firstnode = (*nodes)[0];
        int foundit = 0;
        if (noderowmap->MyGID(firstnode)) foundit = 1;
        int found=0;
        noderowmap->Comm().SumAll(&foundit,&found,1);
        if (found)
          actdis->SetCondition("VolumeNeumann",curr->second);
      }
    }  // for (int j=0;j<field[i].ndis;j++)
  } // for (int i=0; i<genprob.numfld; i++)
  return;
} /* end of input_conditions */

/*----------------------------------------------------------------------*
 | add node ids to a condition                            m.gee 01/07   |
 *----------------------------------------------------------------------*/
void add_nodeids_to_condition(const int id, RefCountPtr<DRT::Condition> cond,
                                     const vector<int> nd_fenode,
                                     const vector<vector<int> > d_fenode)
{
  DSTraceHelper dst("add_nodeids_to_condition");

  // vector of finite element node ids in this node set
  const vector<int>& nodes = d_fenode[id];

  // add the list of nodal ids to the condition
  cond->Add("Node Ids",nodes);
  return;
} // add_nodeids_to_condition



/*----------------------------------------------------------------------*
 | input of design node neumann conditions                m.gee 01/07   |
 *----------------------------------------------------------------------*/
void input_point_neum(multimap<int,RefCountPtr<DRT::Condition> >& pnmap)
{
  DSTraceHelper dst("input_point_neum");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of nodal dirichlet conditions */
  if (frfind("--DESIGN POINT NEUMANN CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design points with conditions */
  int ierr=0;
  int ndnode=0;
  frint("DPOINT",&ndnode,&ierr);
  dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
  frread();

  /*-------------------------------------- start reading the design nodes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
    int dnodeid = -1;
    frint("E",&dnodeid,&ierr);
    dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
    dnodeid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-nodal neumann conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    neum_onoff(MAXDOFPERNODE);
    vector<double> neum_val(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      neum_onoff[i] = 0;
      neum_val[i]   = 0.0;
    }

    //---------------------------------- read the curve number of "none"
    int curve=0;
    char buffer[200];
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
    if (strncmp(buffer,"none",4)==0)
    {
       curve = -1;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-nodal neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
       curve--;
       dsassert(ierr==1,"Cannot read design-nodal neumann conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
    }


    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
             rcp(new DRT::Condition(dnodeid,DRT::Condition::PointNeumann,false,
                                    DRT::Condition::Point));

    // read whether load is on surface (shells)
    frchk("Mid",&ierr);
    if (ierr) condition->Add("surface","mid");
    frchk("Top",&ierr);
    if (ierr) condition->Add("surface","top");
    frchk("Bot",&ierr);
    if (ierr) condition->Add("surface","bot");

    // add stuff to boundary condition
    condition->Add("onoff",neum_onoff);
    condition->Add("val",neum_val);
    condition->Add("curve",&curve,1);

    //------------------------------- put condition in map of conditions
    pnmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dnodeid,condition));
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_point_neum

/*----------------------------------------------------------------------*
 | input of design line neumann conditions                m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_line_neum(multimap<int,RefCountPtr<DRT::Condition> >& lnmap)
{
  DSTraceHelper dst("input_line_neum");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of line dirichlet conditions */
  if (frfind("--DESIGN LINE NEUMANN CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design lines with conditions */
  int ierr=0;
  int ndline=0;
  frint("DLINE",&ndline,&ierr);
  dsassert(ierr==1,"Cannot read design-line neumann conditions");
  frread();

  /*------------------------------------- start reading the design lines */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design line Id */
    int dlineid = -1;
    frint("E",&dlineid,&ierr);
    dsassert(ierr==1,"Cannot read design-line neumann conditions");
    dlineid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-line neumann conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    neum_onoff(MAXDOFPERNODE);
    vector<double> neum_val(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      neum_onoff[i] = 0;
      neum_val[i]   = 0.0;
    }

    //---------------------------------- read the curve number of "none"
    int curve=0;
    char buffer[200];
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design-line neumann conditions");
    if (strncmp(buffer,"none",4)==0)
    {
       curve = -1;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-line neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
       curve--;
       dsassert(ierr==1,"Cannot read design-line neumann conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
    }


    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
              rcp(new DRT::Condition(dlineid,DRT::Condition::LineNeumann,true,
                                     DRT::Condition::Line));

    // read whether load is on surface (shells)
    condition->Add("type","neum_live");
    frchk("Live",&ierr);
    if (ierr) condition->Add("type","neum_live");
    frchk("Dead",&ierr);
    if (ierr) condition->Add("type","neum_dead");
    frchk("PrescribedDomainLoad",&ierr);
    if (ierr) condition->Add("type","pres_domain_load");
    frchk("constHydro_z",&ierr);
    if (ierr) condition->Add("type","neum_consthydro_z");
    frchk("increaseHydro_z",&ierr);
    if (ierr) condition->Add("type","neum_increhydro_z");
    frchk("orthopressure",&ierr);
    if (ierr) condition->Add("type","neum_orthopressure");
    frchk("LAS",&ierr);
    if (ierr) condition->Add("type","neum_LAS");
    /*----------- read if load is applied on surface -> shell elements */
    frchk("Mid",&ierr);
    if (ierr) condition->Add("surface","mid");
    frchk("Top",&ierr);
    if (ierr) condition->Add("surface","top");
    frchk("Bot",&ierr);
    if (ierr) condition->Add("surface","bot");

    // add stuff to boundary condition
    condition->Add("onoff",neum_onoff);
    condition->Add("val",neum_val);
    condition->Add("curve",&curve,1);

    //------------------------------- put condition in map of conditions
    lnmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dlineid,condition));

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_line_neum




/*----------------------------------------------------------------------*
 | input of design surface neumann conditions             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_surf_neum(multimap<int,RefCountPtr<DRT::Condition> >& snmap)
{
  DSTraceHelper dst("input_surf_neum");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*----------------- find the beginning of surface dirichlet conditions */
  if (frfind("--DESIGN SURF NEUMANN CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design surfs with conditions */
  int ierr=0;
  int ndsurf=0;
  frint("DSURF",&ndsurf,&ierr);
  dsassert(ierr==1,"Cannot read design-surface neumann conditions");
  frread();

  /*------------------------------------- start reading the design surfs */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design surf Id */
    int dsurfid = -1;
    frint("E",&dsurfid,&ierr);
    dsassert(ierr==1,"Cannot read design-surface neumann conditions");
    dsurfid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-surface neumann conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    neum_onoff(MAXDOFPERNODE);
    vector<double> neum_val(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      neum_onoff[i] = 0;
      neum_val[i]   = 0.0;
    }

    //---------------------------------- read the curve number of "none"
    int curve=0;
    char buffer[200];
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design-surface neumann conditions");
    if (strncmp(buffer,"none",4)==0)
    {
       curve = -1;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-surface neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
       curve--;
       dsassert(ierr==1,"Cannot read design-surface neumann conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
    }


    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
           rcp(new DRT::Condition(dsurfid,DRT::Condition::SurfaceNeumann,true,
                                  DRT::Condition::Surface));

    // read whether load is on surface (shells)
    condition->Add("type","neum_live");
    frchk("Live",&ierr);
    if (ierr) condition->Add("type","neum_live");
    frchk("Dead",&ierr);
    if (ierr) condition->Add("type","neum_dead");
    frchk("PrescribedDomainLoad",&ierr);
    if (ierr) condition->Add("type","pres_domain_load");
    frchk("constHydro_z",&ierr);
    if (ierr) condition->Add("type","neum_consthydro_z");
    frchk("increaseHydro_z",&ierr);
    if (ierr) condition->Add("type","neum_increhydro_z");
    frchk("orthopressure",&ierr);
    if (ierr) condition->Add("type","neum_orthopressure");
    frchk("LAS",&ierr);
    if (ierr) condition->Add("type","neum_LAS");
    /*----------- read if load is applied on surface -> shell elements */
    frchk("Mid",&ierr);
    if (ierr) condition->Add("surface","mid");
    frchk("Top",&ierr);
    if (ierr) condition->Add("surface","top");
    frchk("Bot",&ierr);
    if (ierr) condition->Add("surface","bot");

    // add stuff to boundary condition
    condition->Add("onoff",neum_onoff);
    condition->Add("val",neum_val);
    condition->Add("curve",&curve,1);

    //------------------------------- put condition in map of conditions
    snmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dsurfid,condition));

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_surf_neum




/*----------------------------------------------------------------------*
 | input of design volume neumann conditions             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_vol_neum(multimap<int,RefCountPtr<DRT::Condition> >& vnmap)
{
  DSTraceHelper dst("input_vol_neum");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*------------------ find the beginning of volume dirichlet conditions */
  if (frfind("--DESIGN VOL NEUMANN CONDITIONS")==0) return;
  frread();

  /*---------------------- read number of design volumes with conditions */
  int ierr=0;
  int ndvol=0;
  frint("DVOL",&ndvol,&ierr);
  dsassert(ierr==1,"Cannot read design-volume neumann conditions");
  frread();

  /*----------------------------------- start reading the design volumes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*---------------------------------------- read the design volume Id */
    int dvolid = -1;
    frint("E",&dvolid,&ierr);
    dsassert(ierr==1,"Cannot read design-volume neumann conditions");
    dvolid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-volume neumann conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    neum_onoff(MAXDOFPERNODE);
    vector<double> neum_val(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      neum_onoff[i] = 0;
      neum_val[i]   = 0.0;
    }

    //---------------------------------- read the curve number of "none"
    int curve=0;
    char buffer[200];
    ierr=sscanf(colptr," %s ",buffer);
    dsassert(ierr==1,"Cannot read design-volume neumann conditions");
    if (strncmp(buffer,"none",4)==0)
    {
       curve = -1;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-volume neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
       curve--;
       dsassert(ierr==1,"Cannot read design-volume neumann conditions");
       colptr = strpbrk(colptr,"1234567890");
       colptr++;
    }


    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        neum_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
            rcp(new DRT::Condition(dvolid,DRT::Condition::VolumeNeumann,true,
                                   DRT::Condition::Volume));

    // read whether load is on surface (shells)
    condition->Add("type","neum_dead");
    frchk("Dead",&ierr);
    if (ierr) condition->Add("type","neum_dead");
    frchk("LAS",&ierr);
    if (ierr) condition->Add("type","neum_LAS");

    // add stuff to boundary condition
    condition->Add("onoff",neum_onoff);
    condition->Add("val",neum_val);
    condition->Add("curve",&curve,1);

    //------------------------------- put condition in map of conditions
    vnmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dvolid,condition));

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_vol_neum


/*----------------------------------------------------------------------*
 | input of design node dirichlet conditions              m.gee 01/07   |
 *----------------------------------------------------------------------*/
void input_point_dirich(multimap<int,RefCountPtr<DRT::Condition> >& pdmap)
{
  DSTraceHelper dst("input_point_dirich");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of nodal dirichlet conditions */
  if (frfind("--DESIGN POINT DIRICH CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design points with conditions */
  int ierr=0;
  int ndnode=0;
  frint("DPOINT",&ndnode,&ierr);
  dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
  frread();

  /*-------------------------------------- start reading the design nodes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
    int dnodeid = -1;
    frint("E",&dnodeid,&ierr);
    dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
    dnodeid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    dirich_onoff(MAXDOFPERNODE);
    vector<double> dirich_val(MAXDOFPERNODE);
    vector<int>    dirich_curve(MAXDOFPERNODE);
    vector<int>    dirich_funct(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      dirich_onoff[i] = 0;
      dirich_val[i]   = 0.0;
      dirich_curve[i] = -1;
      dirich_funct[i] = 0;
    }

    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // read curve number or 'none'
    for (int i=0; i<numread; ++i)
    {
      char buffer[200];
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
      if (strncmp(buffer,"none",4)==0)
      {
        colptr = strstr(colptr,"none");
        dsassert(colptr!=NULL,"Cannot read design-nodal dirichlet conditions");
        colptr += 4;
      }
      else
      {
        ierr=1;
        if (i < MAXDOFPERNODE)
        {
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
          dirich_curve[i]--;
        }
        dsassert(ierr==1,"Cannot read design-nodal dirichlet conditions");
        colptr = strpbrk(colptr,"1234567890");
        colptr++;
      }
    }

    // read function number
    for (int i=0; i<numread; ++i)
    {
      if (i < MAXDOFPERNODE)
        dirich_funct[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);
    }

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
          rcp(new DRT::Condition(dnodeid,DRT::Condition::PointDirichlet,false,
                                 DRT::Condition::Point));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);

    //---------------------- add the condition to the map of all conditions
    pdmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dnodeid,condition));

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_point_dirich

/*----------------------------------------------------------------------*
 | input of design line dirichlet conditions              m.gee 01/07   |
 *----------------------------------------------------------------------*/
void input_line_dirich(multimap<int,RefCountPtr<DRT::Condition> >& ldmap)
{
  DSTraceHelper dst("input_line_dirich");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of line dirichlet conditions */
  if (frfind("--DESIGN LINE DIRICH CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design lines with conditions */
  int ierr=0;
  int ndline=0;
  frint("DLINE",&ndline,&ierr);
  dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
  frread();

  /*-------------------------------------- start reading the design lines */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design line Id */
    int dlineid = -1;
    frint("E",&dlineid,&ierr);
    dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
    dlineid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    dirich_onoff(MAXDOFPERNODE);
    vector<double> dirich_val(MAXDOFPERNODE);
    vector<int>    dirich_curve(MAXDOFPERNODE);
    vector<int>    dirich_funct(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      dirich_onoff[i] = 0;
      dirich_val[i]   = 0.0;
      dirich_curve[i] = -1;
      dirich_funct[i] = 0;
    }

    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // read curve number or 'none'
    for (int i=0; i<numread; ++i)
    {
      char buffer[200];
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
      if (strncmp(buffer,"none",4)==0)
      {
        colptr = strstr(colptr,"none");
        dsassert(colptr!=NULL,"Cannot read design-line dirichlet conditions");
        colptr += 4;
      }
      else
      {
        ierr=1;
        if (i < MAXDOFPERNODE)
        {
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
          dirich_curve[i]--;
        }
        dsassert(ierr==1,"Cannot read design-line dirichlet conditions");
        colptr = strpbrk(colptr,"1234567890");
        colptr++;
      }
    }

    // read function number
    for (int i=0; i<numread; ++i)
      if (i < MAXDOFPERNODE)
        dirich_funct[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
           rcp(new DRT::Condition(dlineid,DRT::Condition::LineDirichlet,false,
                                  DRT::Condition::Line));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);

    //---------------------- add the condition to the map of all conditions
    ldmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dlineid,condition));
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_line_dirich

/*----------------------------------------------------------------------*
 | input of design surface dirichlet conditions           m.gee 01/07   |
 *----------------------------------------------------------------------*/
void input_surf_dirich(multimap<int,RefCountPtr<DRT::Condition> >& sdmap)
{
  DSTraceHelper dst("input_surf_dirich");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*----------------- find the beginning of surface dirichlet conditions */
  if (frfind("--DESIGN SURF DIRICH CONDITIONS")==0) return;
  frread();

  /*------------------------ read number of design surfs with conditions */
  int ierr=0;
  int ndsurf=0;
  frint("DSURF",&ndsurf,&ierr);
  dsassert(ierr==1,"Cannot read design-surface dirichlet conditions");
  frread();

  /*------------------------------------- start reading the design surfs */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design surf Id */
    int dsurfid = -1;
    frint("E",&dsurfid,&ierr);
    dsassert(ierr==1,"Cannot read design-surface dirichlet conditions");
    dsurfid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-surface dirichlet conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    dirich_onoff(MAXDOFPERNODE);
    vector<double> dirich_val(MAXDOFPERNODE);
    vector<int>    dirich_curve(MAXDOFPERNODE);
    vector<int>    dirich_funct(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      dirich_onoff[i] = 0;
      dirich_val[i]   = 0.0;
      dirich_curve[i] = -1;
      dirich_funct[i] = 0;
    }

    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // read curve number or 'none'
    for (int i=0; i<numread; ++i)
    {
      char buffer[200];
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-surface dirichlet conditions");
      if (strncmp(buffer,"none",4)==0)
      {
        colptr = strstr(colptr,"none");
        dsassert(colptr!=NULL,"Cannot read design-surface dirichlet conditions");
        colptr += 4;
      }
      else
      {
        ierr=1;
        if (i < MAXDOFPERNODE)
        {
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
          dirich_curve[i]--;
        }
        dsassert(ierr==1,"Cannot read design-line surface conditions");
        colptr = strpbrk(colptr,"1234567890");
        colptr++;
      }
    }

    // read function number
    for (int i=0; i<numread; ++i)
      if (i < MAXDOFPERNODE)
        dirich_funct[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
        rcp(new DRT::Condition(dsurfid,DRT::Condition::SurfaceDirichlet,false,
                               DRT::Condition::Surface));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);

    //--------------------------------- add condition to map of conditions
    sdmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dsurfid,condition));
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_surf_dirich

/*----------------------------------------------------------------------*
 | input of design volume dirichlet conditions            m.gee 01/07   |
 *----------------------------------------------------------------------*/
void input_vol_dirich(multimap<int,RefCountPtr<DRT::Condition> >& vdmap)
{
  DSTraceHelper dst("input_vol_dirich");

  // currently, we always have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*------------------ find the beginning of volume dirichlet conditions */
  if (frfind("--DESIGN VOL DIRICH CONDITIONS")==0) return;
  frread();

  /*---------------------- read number of design volumes with conditions */
  int ierr=0;
  int ndvol=0;
  frint("DVOL",&ndvol,&ierr);
  dsassert(ierr==1,"Cannot read volume dirichlet conditions");
  frread();

  /*----------------------------------- start reading the design volumes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*---------------------------------------- read the design volume Id */
    int dvolid = -1;
    frint("E",&dvolid,&ierr);
    dsassert(ierr==1,"Cannot read design-volume dirichlet conditions");
    dvolid--;
    /*--------------------------------- move pointer behind the "-" sign */
    char* colptr = strstr(allfiles.actplace,"-");
    dsassert(colptr!=NULL,"Cannot read design-volume dirichlet conditions");
    colptr++;

    //------------------------------- define some temporary reading vectors
    vector<int>    dirich_onoff(MAXDOFPERNODE);
    vector<double> dirich_val(MAXDOFPERNODE);
    vector<int>    dirich_curve(MAXDOFPERNODE);
    vector<int>    dirich_funct(MAXDOFPERNODE);
    for (int i=0; i<MAXDOFPERNODE; ++i)
    {
      dirich_onoff[i] = 0;
      dirich_val[i]   = 0.0;
      dirich_curve[i] = -1;
      dirich_funct[i] = 0;
    }

    /* NOTE: number of read values = 6  does not need to be */
    /*       equivalent to the MAXDOFPERNODE -> e.g. for shell9! sh 12/02 */
    // read on/off toggles
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_onoff[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // read values
    for (int i=0; i<numread; ++i)
      if (i<MAXDOFPERNODE)
        dirich_val[i] = strtod(colptr,&colptr);
      else
        strtod(colptr,&colptr);

    // read curve number or 'none'
    for (int i=0; i<numread; ++i)
    {
      char buffer[200];
      ierr=sscanf(colptr," %s ",buffer);
      dsassert(ierr==1,"Cannot read design-volume dirichlet conditions");
      if (strncmp(buffer,"none",4)==0)
      {
        colptr = strstr(colptr,"none");
        dsassert(colptr!=NULL,"Cannot read design-volume dirichlet conditions");
        colptr += 4;
      }
      else
      {
        ierr=1;
        if (i < MAXDOFPERNODE)
        {
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
          dirich_curve[i]--;
        }
        dsassert(ierr==1,"Cannot read design-volume conditions");
        colptr = strpbrk(colptr,"1234567890");
        colptr++;
      }
    }

    // read function number
    for (int i=0; i<numread; ++i)
      if (i < MAXDOFPERNODE)
        dirich_funct[i] = strtol(colptr,&colptr,10);
      else
        strtol(colptr,&colptr,10);

    // create boundary condition
    RefCountPtr<DRT::Condition> condition =
          rcp(new DRT::Condition(dvolid,DRT::Condition::VolumeDirichlet,false,
                                 DRT::Condition::Volume));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);

    //--------------------------------- put condition in map of conditions
    vdmap.insert(pair<int,RefCountPtr<DRT::Condition> >(dvolid,condition));
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_vol_dirich


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
