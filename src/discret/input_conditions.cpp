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
static void input_designnode_dirich(RefCountPtr<CCADISCRETIZATION::DesignDiscretization> designdis);

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
  /*-------------------------------- input of local co-ordinate systems */
  // not yet implemented as these are stored inside design objects

  if (design)
  {
    RefCountPtr<CCADISCRETIZATION::Design>* tmp = 
               (RefCountPtr<CCADISCRETIZATION::Design>*)design->ccadesign;
    CCADISCRETIZATION::Design& ccadesign = **tmp;
    RefCountPtr<CCADISCRETIZATION::DesignDiscretization> designlines = ccadesign[0];
    RefCountPtr<CCADISCRETIZATION::DesignDiscretization> designsurfs = ccadesign[1];
    RefCountPtr<CCADISCRETIZATION::DesignDiscretization> designvols  = ccadesign[2];
  
    if (ccadesign.Comm().MyPID()==0)
    {
      //------------------------- input of design node dirichlet conditions
      input_designnode_dirich(designlines);
    
    
    
    
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
void input_designnode_dirich(RefCountPtr<CCADISCRETIZATION::DesignDiscretization> designdis)
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
    CCADISCRETIZATION::Condition condition(CCADISCRETIZATION::Condition::condition_Dirichlet);
    condition.Add("onoff",dirich_onoff);
    condition.Add("val",dirich_val);
    condition.Add("curve",dirich_curve);
    condition.Add("funct",dirich_funct);
    cout << condition;
    
    //------------------------------- get the dnode from the discretization
    CCADISCRETIZATION::Node* node = designdis->gNode(dnodeid);
    if (node==NULL) dserror("Cannot find design node");


    //-------------------------------------------------- read the next line
    frread();
  } // while(strncmp(allfiles.actplace,"------",6)!=0)
  return;
} // input_designnode_dirich


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
