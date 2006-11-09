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
 | input of design volumes to fe-node topology            m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_design_dvol_fenode_read(vector<vector<int> >& dvol_fenode,
                                   vector<int>& ndvol_fenode)
{
INT    i,ierr;
INT    counter;
INT    dvol;
DSTraceHelper dst("input_design_dvol_fenode_read");
/*----------------------------------------------------------------------*/
/*---------------------------------- count number of nodes on this vol */
if (frfind("--DVOL-NODE TOPOLOGY")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DVOL",&dvol,&ierr);
   dsassert(ierr==1,"Cannot read DVOL-NODE TOPOLOGY");
   ndvol_fenode[dvol-1]++;
   frread();
}
frrewind();
for (i=0; i<(int)ndvol_fenode.size(); i++)
   dvol_fenode[i].resize(ndvol_fenode[i]);
/*------------------------------- find fe-nodes belonging to this dvol */
for (i=0; i<(int)ndvol_fenode.size(); i++)
{
   frfind("--DVOL-NODE TOPOLOGY");
   frread();
   //design->dvol[i].Id=i;
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DVOL",&dvol,&ierr);
      dsassert(ierr==1,"Cannot read DVOL-NODE TOPOLOGY");
      if (dvol-1==i)
      {
         dsassert(counter<ndvol_fenode[i],"Cannot read DVOL-NODE TOPOLOGY");
         frint("NODE",&dvol_fenode[i][counter],&ierr);
         dsassert(ierr==1,"Cannot read DVOL-NODE TOPOLOGY");
         dvol_fenode[i][counter]--;
         counter++;
      }
      frread();
   }
}
/*----------------------------------------------------------------------*/

end:
return;
} /* end of input_design_dvol_fenode_read */


/*----------------------------------------------------------------------*
 | input of design surface to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
void input_design_dsurf_fenode_read(vector<vector<int> >& dsurf_fenode,
                                     vector<int>& ndsurf_fenode)
{
INT    i,ierr;
INT    counter;
INT    dsurf;
DSTraceHelper dst("input_design_dsurf_fenode_read");
/*----------------------------------------------------------------------*/
/*---------------------------------- count number of nodes on this surf */
if (frfind("--DSURF-NODE TOPOLOGY")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DSURF",&dsurf,&ierr);
   dsassert(ierr==1,"Cannot read DSURF-NODE TOPOLOGY");
   ndsurf_fenode[dsurf-1]++;
   frread();
}
frrewind();
for (i=0; i<(int)ndsurf_fenode.size(); i++)
   dsurf_fenode[i].resize(ndsurf_fenode[i]);
/*------------------------------- find fe-nodes belonging to this dsurf */
for (i=0; i<(int)ndsurf_fenode.size(); i++)
{
   frfind("--DSURF-NODE TOPOLOGY");
   frread();
   //design->dsurf[i].Id=i;
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DSURF",&dsurf,&ierr);
      dsassert(ierr==1,"Cannot read DSURF-NODE TOPOLOGY");
      if (dsurf-1==i)
      {
         dsassert(counter<ndsurf_fenode[i],"Cannot read DSURF-NODE TOPOLOGY");
         frint("NODE",&(dsurf_fenode[i][counter]),&ierr);
         dsassert(ierr==1,"Cannot read DSURF-NODE TOPOLOGY");
         dsurf_fenode[i][counter]--;
         counter++;
      }
      frread();
   }
}
/*----------------------------------------------------------------------*/

end:
return;
} /* end of input_design_dsurf_fenode_read */


/*----------------------------------------------------------------------*
 | input of design line  to fe-node topology             m.gee 11/06    |
 *----------------------------------------------------------------------*/
void input_design_dline_fenode_read(vector<vector<int> >& dline_fenode,
                                     vector<int>& ndline_fenode)
{
INT    i,ierr;
INT    counter;
INT    dline;
DSTraceHelper dst("input_design_dline_fenode_read");
/*----------------------------------------------------------------------*/
/*---------------------------------- count number of nodes on this line */
if (frfind("--DLINE-NODE TOPOLOGY")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("DLINE",&dline,&ierr);
   dsassert(ierr==1,"Cannot read DLINE-NODE TOPOLOGY");
   ndline_fenode[dline-1]++;
   frread();
}
frrewind();
for (i=0; i<(int)ndline_fenode.size(); i++)
   dline_fenode[i].resize(ndline_fenode[i]);
/*------------------------------- find fe-nodes belonging to this dline */
for (i=0; i<(int)ndline_fenode.size(); i++)
{
   frfind("--DLINE-NODE TOPOLOGY");
   frread();
   //design->dline[i].Id=i;
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DLINE",&dline,&ierr);
      dsassert(ierr==1,"Cannot read DLINE-NODE TOPOLOGY");
      if (dline-1==i)
      {
         dsassert(counter<ndline_fenode[i],"Cannot read DLINE-NODE TOPOLOGY");
         frint("NODE",&(dline_fenode[i][counter]),&ierr);
         dsassert(ierr==1,"Cannot read DLINE-NODE TOPOLOGY");
         dline_fenode[i][counter]--;
         counter++;
      }
      frread();
   }
}
/*----------------------------------------------------------------------*/

end:
return;
} /* end of input_design_dline_fenode_read */

/*----------------------------------------------------------------------*
 | input of design nodes  to fe-node topology             m.gee 11/06   |
 *----------------------------------------------------------------------*/
void input_design_dpoint_fenode_read(vector<vector<int> >& dnode_fenode,
                                     vector<int>& ndnode_fenode)
{
  DSTraceHelper dst("input_design_dpoint_fenode_read");
/*-------------------------------------------------------------- rewind */
frrewind();
/*------------------------------- find fe-nodes belonging to this dnode */
int ierr=0;
int found=0;
for (int i=0; i<(int)ndnode_fenode.size(); i++)
{
   //design->dnode[i].Id=i;
   if (frfind("--DNODE-NODE TOPOLOGY")==0) goto end;
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      int dnode=0;
      frint("DNODE",&dnode,&ierr);
      if (ierr==1)
      {
         found = 1; /*fount at least one DNODE-NODE TOPOLOGY*/
         if (dnode==i+1)
         {
            ndnode_fenode[i]=1;
            dnode_fenode[i].resize(ndnode_fenode[i]);
            frint("NODE",&(dnode_fenode[i][0]),&ierr);
            dnode_fenode[i][0]--;
            goto nextdnode;
         }
      }
      frread();
   }
   nextdnode:
   frrewind();
}
/*----------------------------------------------------------------------*/
if (found == 0) dserror("Cannot make DNODE-NODE topology");
/*----------------------------------------------------------------------*/

end:
return;
} /* end of input_design_dpoint_fenode_read */



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
