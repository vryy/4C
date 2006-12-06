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


// local methods
static void input_designnode_dirich(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designline_dirich(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designsurf_dirich(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designvol_dirich(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designnode_neum(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designline_neum(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designsurf_neum(RefCountPtr<DRT::DesignDiscretization> designdis);
static void input_designvol_neum(RefCountPtr<DRT::DesignDiscretization> designdis);

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

  if (design)
  {
    RefCountPtr<DRT::Design>* tmp = (RefCountPtr<DRT::Design>*)design->ccadesign;
    DRT::Design& ccadesign = **tmp;
    RefCountPtr<DRT::DesignDiscretization> designlines = ccadesign[0];
    RefCountPtr<DRT::DesignDiscretization> designsurfs = ccadesign[1];
    RefCountPtr<DRT::DesignDiscretization> designvols  = ccadesign[2];
  
    if (ccadesign.Comm().MyPID()==0)
    {
      //------------------------- input of design node dirichlet conditions
      input_designnode_dirich(designlines);
      //------------------------- input of design line dirichlet conditions
      input_designline_dirich(designlines);
      //---------------------- input of design surface dirichlet conditions
      input_designsurf_dirich(designsurfs);
      //---------------------- input of design volume dirichlet conditions
      input_designvol_dirich(designvols);
      //---------------------------input of design node neumann conditions
      input_designnode_neum(designlines);
      //---------------------------input of design line neumann conditions
      input_designline_neum(designlines);
      //------------------------input of design surface neumann conditions
      input_designsurf_neum(designsurfs);
      //------------------------input of design volume neumann conditions
      input_designvol_neum(designvols);
      // add reading of more conditions here (coupling, fsi etc)
    
    } // if (ccadesign.Comm().MyPID()==0)
    
  }
  else
  {
    dserror("DESIGN-free boundary conditions not yet impl.");
  }


  return;
} /* end of input_conditions */


/*----------------------------------------------------------------------*
 | input of design node dirichlet conditions              m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designnode_dirich(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designnode_dirich");

  // currently, we alsways have 6 values read from file
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
      dirich_curve[i] = 0;
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
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
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
                rcp(new DRT::Condition(DRT::Condition::condition_Dirichlet));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);
    
    //------------------------------- get the dnode from the discretization
    DRT::Node* node = designdis->gNode(dnodeid);
    if (node==NULL) dserror("Cannot find design node");

    //-------------------------------------------- attach condition to node
    node->SetCondition("Dirichlet",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designnode_dirich


/*----------------------------------------------------------------------*
 | input of design line dirichlet conditions              m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designline_dirich(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designline_dirich");

  // currently, we alsways have 6 values read from file
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
    /*------------------------------------------ read the design node Id */
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
      dirich_curve[i] = 0;
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
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
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
                rcp(new DRT::Condition(DRT::Condition::condition_Dirichlet));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);
    
    //------------------------------- get the dline from the discretization
    
    DRT::Element* line = designdis->gElement(dlineid);
    if (line==NULL) dserror("Cannot find design line");

    //-------------------------------------------- attach condition to node
    line->SetCondition("Dirichlet",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designline_dirich



/*----------------------------------------------------------------------*
 | input of design surface dirichlet conditions           m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designsurf_dirich(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designsurf_dirich");

  // currently, we alsways have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of line dirichlet conditions */
  if (frfind("--DESIGN SURF DIRICH CONDITIONS")==0) return;
  frread();
  
  /*------------------------ read number of design surfs with conditions */
  int ierr=0;
  int ndsurf=0;
  frint("DSURF",&ndsurf,&ierr);
  dsassert(ierr==1,"Cannot read design-surface dirichlet conditions");
  frread();
  
  /*-------------------------------------- start reading the design lines */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
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
      dirich_curve[i] = 0;
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
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
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
                 rcp(new DRT::Condition(DRT::Condition::condition_Dirichlet));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);
    
    //------------------------------- get the dsurf from the discretization
    
    DRT::Element* surf = designdis->gElement(dsurfid);
    if (surf==NULL) dserror("Cannot find design surface");

    //-------------------------------------------- attach condition to node
    surf->SetCondition("Dirichlet",condition);
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designsurf_dirich




/*----------------------------------------------------------------------*
 | input of design volume dirichlet conditions            m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designvol_dirich(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designvol_dirich");

  // currently, we alsways have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of line dirichlet conditions */
  if (frfind("--DESIGN VOL DIRICH CONDITIONS")==0) return;
  frread();
  
  /*------------------------ read number of design surfs with conditions */
  int ierr=0;
  int ndvol=0;
  frint("DVOL",&ndvol,&ierr);
  dsassert(ierr==1,"Cannot read design-surface dirichlet conditions");
  frread();
  
  /*-------------------------------------- start reading the design lines */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
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
      dirich_curve[i] = 0;
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
          ierr=sscanf(colptr," %d ",&dirich_curve[i]);
        dsassert(ierr==1,"Cannot read design-line volume conditions");
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
                rcp(new DRT::Condition(DRT::Condition::condition_Dirichlet));
    condition->Add("onoff",dirich_onoff);
    condition->Add("val",dirich_val);
    condition->Add("curve",dirich_curve);
    condition->Add("funct",dirich_funct);
    
    //------------------------------- get the dvol from the discretization
    DRT::Element* vol = designdis->gElement(dvolid);
    if (vol==NULL) dserror("Cannot find design volume");

    //-------------------------------------------- attach condition to node
    vol->SetCondition("Dirichlet",condition);
    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designvol_dirich


/*----------------------------------------------------------------------*
 | input of design node neumann conditions                m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designnode_neum(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designnode_neum");

  // currently, we alsways have 6 values read from file
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
       curve = 0;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-nodal neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
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
                    rcp(new DRT::Condition(DRT::Condition::condition_Neumann));
    
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
    //------------------------------- get the dnode from the discretization
    DRT::Node* node = designdis->gNode(dnodeid);
    if (node==NULL) dserror("Cannot find design node");

    //-------------------------------------------- attach condition to node
    node->SetCondition("Neumann",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designnode_neum



/*----------------------------------------------------------------------*
 | input of design line neumann conditions                m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designline_neum(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designline_neum");

  // currently, we alsways have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of nodal dirichlet conditions */
  if (frfind("--DESIGN LINE NEUMANN CONDITIONS")==0) return;
  frread();
  
  /*------------------------ read number of design points with conditions */
  int ierr=0;
  int ndline=0;
  frint("DLINE",&ndline,&ierr);
  dsassert(ierr==1,"Cannot read design-line neumann conditions");
  frread();
  
  /*-------------------------------------- start reading the design nodes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
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
       curve = 0;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-line neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
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
                  rcp(new DRT::Condition(DRT::Condition::condition_Neumann));
    
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

    //------------------------------- get the dline from the discretization
    DRT::Element* line = designdis->gElement(dlineid);
    if (line==NULL) dserror("Cannot find design line");

    //-------------------------------------------- attach condition to node
    line->SetCondition("Neumann",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designline_neum




/*----------------------------------------------------------------------*
 | input of design surface neumann conditions             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designsurf_neum(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designsurf_neum");

  // currently, we alsways have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of nodal dirichlet conditions */
  if (frfind("--DESIGN SURF NEUMANN CONDITIONS")==0) return;
  frread();
  
  /*------------------------ read number of design points with conditions */
  int ierr=0;
  int ndsurf=0;
  frint("DSURF",&ndsurf,&ierr);
  dsassert(ierr==1,"Cannot read design-surface neumann conditions");
  frread();
  
  /*-------------------------------------- start reading the design nodes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
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
       curve = 0;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-surface neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
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
                  rcp(new DRT::Condition(DRT::Condition::condition_Neumann));
    
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

    //------------------------------- get the dsurf from the discretization
    DRT::Element* surf = designdis->gElement(dsurfid);
    if (surf==NULL) dserror("Cannot find design surface");

    //-------------------------------------------- attach condition to node
    surf->SetCondition("Neumann",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designsurf_neum




/*----------------------------------------------------------------------*
 | input of design volume neumann conditions             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_designvol_neum(RefCountPtr<DRT::DesignDiscretization> designdis)
{
  DSTraceHelper dst("input_designvol_neum");

  // currently, we alsways have 6 values read from file
  // this might change at some point
  const int numread = 6;

  /*-------------------- find the beginning of nodal dirichlet conditions */
  if (frfind("--DESIGN VOL NEUMANN CONDITIONS")==0) return;
  frread();
  
  /*------------------------ read number of design points with conditions */
  int ierr=0;
  int ndvol=0;
  frint("DVOL",&ndvol,&ierr);
  dsassert(ierr==1,"Cannot read design-volume neumann conditions");
  frread();
  
  /*-------------------------------------- start reading the design nodes */
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /*------------------------------------------ read the design node Id */
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
       curve = 0;
       colptr = strstr(allfiles.actplace,"none");
       dsassert(colptr!=NULL,"Cannot read design-volume neumann conditions");
       colptr += 4;
    }
    else
    {
       ierr=sscanf(colptr," %d ",&curve);
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
                  rcp(new DRT::Condition(DRT::Condition::condition_Neumann));
    
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

    //------------------------------- get the dvol from the discretization
    DRT::Element* vol = designdis->gElement(dvolid);
    if (vol==NULL) dserror("Cannot find design volume");

    //-------------------------------------------- attach condition to node
    vol->SetCondition("Neumann",condition);

    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designvol_neum

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
