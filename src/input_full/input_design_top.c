/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
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
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;
/* global variable: flag for the creation of a second discretisation */
extern INT      create_dis;

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
 | prototypes of static routines in this file                      3/02 |
 *----------------------------------------------------------------------*/
static void inpdesign_dpoint_fenode_read(void);
static void inpdesign_dline_fenode_read(void);
static void inpdesign_dsurf_fenode_read(void);
static void inpdesign_dvol_fenode_read(void);

static void inpdesign_dpoint_fenode(void);
static void inpdesign_dline_feline(void);
static void inpdesign_dsurf_fesurf(void);
static void inpdesign_dvol_fevol(void);

/*----------------------------------------------------------------------*
 | global variables in this file in this file                      3/02 |
 *----------------------------------------------------------------------*/
static INT *ndnode_fenode;
static INT **dnode_fenode;

static INT *ndline_fenode;
static INT **dline_fenode;

static INT *ndsurf_fenode;
static INT **dsurf_fenode;

static INT *ndvol_fenode; 
static INT **dvol_fenode; 
  
/*----------------------------------------------------------------------*
 | create the connectivity of the design                           1/02 |
 *----------------------------------------------------------------------*/
void inpdesign_topology_design()
{
INT       i,j,k;
DNODE    *actdnode;
DLINE    *actdline;
DSURF    *actdsurf;
DVOL     *actdvol;
INT       nodeid;
INT       lineid;
INT       surfid;

#ifdef DEBUG 
dstrc_enter("inpdesign_topology_design");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------- init the topology of the design to zero */
for (i=0; i<design->ndnode; i++) 
{
   design->dnode[i].ndline = 0;
   design->dnode[i].dline  = NULL;
}
for (i=0; i<design->ndline; i++) 
{
   design->dline[i].ndsurf = 0;
   design->dline[i].dsurf  = NULL;
/*   design->dline[i].dnode  = NULL;*/
}
for (i=0; i<design->ndsurf; i++) 
{
   design->dsurf[i].ndvol  = 0;
   design->dsurf[i].dvol   = NULL;
   design->dsurf[i].dline  = NULL;
}
for (i=0; i<design->ndvol; i++)
{
   design->dvol[i].dsurf   = NULL;
}
/*----------------------------- make connectivity from dlines to dnodes */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*actdline->dnode = (DNODE**)CCACALLOC(2,sizeof(DNODE*));*/
   for (k=0; k<2; k++)
   {
      nodeid=actdline->my_dnodeId[k];
      /* j = nodeid; ??? */
      for (j=0; j<design->ndnode; j++)
      {
         if (design->dnode[j].Id==nodeid) break;
      }
      actdline->dnode[k] = &(design->dnode[j]);
      actdline->dnode[k]->ndline++;
   }
}
/*----------------------------- make connectivity from dnodes to dlines */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   for (j=0; j<2; j++)
   {
      actdnode = actdline->dnode[j];
      if (actdnode->dline==NULL)
      {
        /* when there are many nodes it would be better to allocate
         * just one junk of memory. */
         actdnode->dline = (DLINE**)CCACALLOC(actdnode->ndline,sizeof(DLINE*));
         actdnode->dline[0] = actdline;
      }
      else
      {
         k=0;
         while (k<actdnode->ndline && actdnode->dline[k]!=NULL) k++;
         if (k==actdnode->ndline-1 && actdnode->dline[k]!=NULL)
           dserror("Cannot make dnode to dline topology");
         actdnode->dline[k] = actdline;
      }
   }
}
/*-------------------------- make connectivity from dsurfaces to dlines */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   actdsurf->dline = (DLINE**)CCACALLOC(actdsurf->ndline,sizeof(DLINE*));
   for (k=0; k<actdsurf->my_dlineId.fdim; k++)
   {
      lineid = actdsurf->my_dlineId.a.ia[k][0];
      for (j=0; j<design->ndline; j++)
      {
         if (design->dline[j].Id==lineid) break;
      }
      actdsurf->dline[k] = &(design->dline[j]);
      actdsurf->dline[k]->ndsurf++;
   }
}
/*-------------------------- make connectivity from dlines to dsurfaces */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   for (j=0; j<actdsurf->ndline; j++)
   {
      actdline = actdsurf->dline[j];
      if (actdline->dsurf==NULL)
      {
         actdline->dsurf = (DSURF**)CCACALLOC(actdline->ndsurf,sizeof(DSURF*));
         actdline->dsurf[0] = actdsurf;
      }
      else
      {
         k=0;
         while (k<actdline->ndsurf && actdline->dsurf[k]!=NULL) k++;
         if (k==actdline->ndsurf-1 && actdline->dsurf[k]!=NULL)
           dserror("Cannot make dline to dsurf topology");
         actdline->dsurf[k] = actdsurf;
      }
   }
}
/*----------------------  make connectivity from dvolumes to dsurfaces */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   actdvol->dsurf = (DSURF**)CCACALLOC(actdvol->ndsurf,sizeof(DSURF**));
   for (k=0; k<actdvol->ndsurf; k++)
   {
      surfid = actdvol->my_dsurfId.a.ia[k][0];
      for (j=0; j<design->ndsurf; j++)
      {
         if (design->dsurf[j].Id==surfid) break;
      }
      actdvol->dsurf[k] = &(design->dsurf[j]);
      actdvol->dsurf[k]->ndvol++;
   }
}
/*----------------------- make connectivity from dsurfaces to dvolumes */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   for (j=0; j<actdvol->ndsurf; j++)
   {
      actdsurf = actdvol->dsurf[j];
      if (actdsurf->dvol==NULL)
      {
         actdsurf->dvol = (DVOL**)CCACALLOC(actdsurf->ndvol,sizeof(DVOL*));
         actdsurf->dvol[0] = actdvol;
      }
      else
      {
         k=0;
         while (k<actdsurf->ndvol && actdsurf->dvol[k]!=NULL) k++;
         if (k==actdsurf->ndvol-1 && actdsurf->dvol[k]!=NULL)
           dserror("Cannot make dsurf to dvol topology");
         actdsurf->dvol[k] = actdvol;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_topology_design */



/*----------------------------------------------------------------------*
 | create the connectivity between the design and the fields m.gee 4/01 |
 *----------------------------------------------------------------------*/
void inpdesign_topology_fe()
{
INT i,j;

#ifdef DEBUG 
dstrc_enter("inpdesign_topology_fe");
#endif

/*--------------------------------------------------read fe-design info */
ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
dnode_fenode  = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
dline_fenode  = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
dsurf_fenode  = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

ndvol_fenode = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
dvol_fenode  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;
/*----------------------------- read the fe-nodes on each design object */
  inpdesign_dpoint_fenode_read();
  inpdesign_dline_fenode_read();
  inpdesign_dsurf_fenode_read();
  inpdesign_dvol_fenode_read();
/*----------------------- make the topology between DESIGN and GEOMETRY */
  inpdesign_dpoint_fenode();
  inpdesign_dline_feline();
  inpdesign_dsurf_fesurf();
  inpdesign_dvol_fevol();
  /* make topolopy for gentities on the second discretisation */
  if (create_dis ==1)
  {
    for (i=0; i<design->ndnode; i++)
      for ( j=0; j<ndnode_fenode[i]; j++)
        dnode_fenode[i][j] += genprob.nodeshift;
    for (i=0; i<design->ndline; i++)
      for ( j=0; j<ndline_fenode[i]; j++)
        dline_fenode[i][j] += genprob.nodeshift;
    for (i=0; i<design->ndsurf; i++)
      for ( j=0; j<ndsurf_fenode[i]; j++)
        dsurf_fenode[i][j] += genprob.nodeshift;
    for (i=0; i<design->ndvol; i++)
      for ( j=0; j<ndvol_fenode[i]; j++)
        dvol_fenode[i][j] += genprob.nodeshift;
    inpdesign_dpoint_fenode();
    inpdesign_dline_feline();
    inpdesign_dsurf_fesurf();
    inpdesign_dvol_fevol();
  }
/*------------------------------------------------------------- tidy up */
for (i=0; i<design->ndnode; i++) CCAFREE(dnode_fenode[i]);
CCAFREE(dnode_fenode);
CCAFREE(ndnode_fenode);
for (i=0; i<design->ndline; i++) CCAFREE(dline_fenode[i]);
CCAFREE(dline_fenode);
CCAFREE(ndline_fenode);
for (i=0; i<design->ndsurf; i++) CCAFREE(dsurf_fenode[i]);
CCAFREE(dsurf_fenode);
CCAFREE(ndsurf_fenode);
for (i=0; i<design->ndvol; i++) CCAFREE(dvol_fenode[i]);
CCAFREE(dvol_fenode);
CCAFREE(ndvol_fenode);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_topology_fe */



/*----------------------------------------------------------------------*
 | input of design nodes  to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dpoint_fenode_read()
{
INT    i,ierr;
INT    dnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_dpoint_fenode_read");
#endif
/*-------------------------------------------------------------- rewind */
frrewind();
/*------------------------------- find fe-nodes belonging to this dnode */
for (i=0; i<design->ndnode; i++)
{
   design->dnode[i].Id=i;
   if (frfind("--DNODE-NODE TOPOLOGY")==0) goto end;
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DNODE",&dnode,&ierr);
      if (ierr==1)
      {
         if (dnode==i+1)
         {
            ndnode_fenode[i]=1;
            dnode_fenode[i] = (INT*)CCAMALLOC(ndnode_fenode[i]*sizeof(INT));
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

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dpoint_fenode_read */


/*----------------------------------------------------------------------*
 | input of design line  to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dline_fenode_read()
{
INT    i,ierr;
INT    counter;
INT    dline;
#ifdef DEBUG 
dstrc_enter("inpdesign_dline_fenode_read");
#endif
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
for (i=0; i<design->ndline; i++)
{
   dline_fenode[i] = (INT*)CCAMALLOC(ndline_fenode[i]*sizeof(INT));
}
/*------------------------------- find fe-nodes belonging to this dline */
for (i=0; i<design->ndline; i++)
{
   frfind("--DLINE-NODE TOPOLOGY");
   frread();
   design->dline[i].Id=i;
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DLINE",&dline,&ierr);
      dsassert(ierr==1,"Cannot read DLINE-NODE TOPOLOGY");
      if (dline-1==i)
      {
         dsassert(counter<ndline_fenode[i],"Cannot read DLINE-NODE TOPOLOGY");
         frint("NODE",&dline_fenode[i][counter],&ierr);
         dsassert(ierr==1,"Cannot read DLINE-NODE TOPOLOGY");
         dline_fenode[i][counter]--;
         counter++;
      }
      frread();
   }
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dline_fenode_read */


/*----------------------------------------------------------------------*
 | input of design surface to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dsurf_fenode_read()
{
INT    i,ierr;
INT    counter;
INT    dsurf;
#ifdef DEBUG 
dstrc_enter("inpdesign_dsurf_fenode_read");
#endif
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
for (i=0; i<design->ndsurf; i++)
{
   dsurf_fenode[i] = (INT*)CCAMALLOC(ndsurf_fenode[i]*sizeof(INT));
}
/*------------------------------- find fe-nodes belonging to this dsurf */
for (i=0; i<design->ndsurf; i++)
{
   frfind("--DSURF-NODE TOPOLOGY");
   frread();
   design->dsurf[i].Id=i;
   counter=0;
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("DSURF",&dsurf,&ierr);
      dsassert(ierr==1,"Cannot read DSURF-NODE TOPOLOGY");
      if (dsurf-1==i)
      {
         dsassert(counter<ndsurf_fenode[i],"Cannot read DSURF-NODE TOPOLOGY");
         frint("NODE",&dsurf_fenode[i][counter],&ierr);
         dsassert(ierr==1,"Cannot read DSURF-NODE TOPOLOGY");
         dsurf_fenode[i][counter]--;
         counter++;
      }
      frread();
   }
}
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dsurf_fenode_read */


/*----------------------------------------------------------------------*
 | input of design volumes to fe-node topology            m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dvol_fenode_read()
{
INT    i,ierr;
INT    counter;
INT    dvol;
#ifdef DEBUG 
dstrc_enter("inpdesign_dvol_fenode_read");
#endif
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
for (i=0; i<design->ndvol; i++)
{
   dvol_fenode[i] = (INT*)CCAMALLOC(ndvol_fenode[i]*sizeof(INT));
}
/*------------------------------- find fe-nodes belonging to this dvol */
for (i=0; i<design->ndvol; i++)
{
   frfind("--DVOL-NODE TOPOLOGY");
   frread();
   design->dvol[i].Id=i;
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
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dvol_fenode_read */





/*----------------------------------------------------------------------*
 | make topology DNODE <-> GNODE                          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dpoint_fenode()
{
INT           i;
INT           nodeId;
DNODE        *actdnode;
GNODE        *actgnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_dpoint_fenode");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<design->ndnode; i++)
{
   actdnode = &(design->dnode[i]);
   nodeId   = dnode_fenode[i][0];
   actgnode = genprob.gnodes[nodeId];
   if (actgnode == NULL) continue;
   actgnode->ondesigntyp = ondnode;
   actgnode->d.dnode = actdnode;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dpoint_fenode */



/*----------------------------------------------------------------------*
 | make topology DLINE <-> GLINE                          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dline_feline()
{
INT           i,j,k;
DLINE        *actdline;
GLINE        *actgline;
GNODE        *actgnode;
INT           nnodeonline;
INT           nodeId;
INT           firstnode,scndnode,thirdnode;
INT           firstmatch,scndmatch,thirdmatch;
INT           linematch;
#ifdef DEBUG 
dstrc_enter("inpdesign_dline_feline");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   nnodeonline = ndline_fenode[i];
   /* loop all GLINES and check whether there nodes are on this actdline */
   for (j=0; j<genprob.ngline; j++)
   {
      linematch=0;
      firstmatch=scndmatch=thirdmatch=0;
      actgline = genprob.glines[j];
      if (actgline->dline) continue;
      /*----------------------- check whether firstnode is on actdline */
      firstnode = actgline->gnode[0]->node->Id;
      for (k=0; k<nnodeonline; k++)
      {
         if (dline_fenode[i][k]==firstnode)
         {
            firstmatch=1;
            break;
         }
      }
      /*------------------------ if firstnode not on actdline continue */
      if (!firstmatch) continue;
      /*------------------------------- check for scndnode on actdline */
      scndnode = actgline->gnode[1]->node->Id;
      for (k=0; k<nnodeonline; k++)
      {
         if (dline_fenode[i][k]==scndnode)
         {
            scndmatch=1;
            break;
         }
      }
      /*------------------------- if scndnode not on actdline continue */
      if (!scndmatch) continue;
      /*----------------- check whether this actgline has 3 or 2 nodes */
      /*-------------------- a 2-node actgline does match the actdline */ 
      if (actgline->ngnode==2) linematch=1;
      /*----------------------- for a 3-node actgline check third node */
      else
      {
         thirdnode = actgline->gnode[2]->node->Id;
         for (k=0; k<nnodeonline; k++)
         {
            if (dline_fenode[i][k]==thirdnode)
            {
               thirdmatch=1;
               break;
            }
         }
         dsassert(thirdmatch==1,"Problems with quadratic GLINEs");
         if (!thirdmatch) continue;
         linematch=1;
      }
      dsassert(linematch==1,"Problems with GLINEs");
      actgline->dline = actdline;
   } /* end loop j over glines */    
}/* end loop i over dlines */
/* loop all dlines again and make pointers from all gnodes to the dlines*/
/*------------- if a gnode already has a pointer to a dnode, do nothing */

  for (i=0; i<design->ndline; i++)
  {
    actdline = &(design->dline[i]);
    nnodeonline = ndline_fenode[i];
    for (j=0; j<nnodeonline; j++)
    {
      nodeId = dline_fenode[i][j];
      actgnode = genprob.gnodes[nodeId];
      if (actgnode == NULL) continue;
      if (actgnode->ondesigntyp != ondnothing) continue;
      actgnode->ondesigntyp = ondline;
      actgnode->d.dline = actdline;
    }
  }/* end loop i over dlines */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dline_feline */



/*----------------------------------------------------------------------*
 | make topology DSURF <-> GSURF                          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dsurf_fesurf()
{
INT           i,j,k;
GSURF        *actgsurf;
DSURF        *actdsurf;
GNODE        *actgnode;
INT           nnodeonsurf;
INT           firstnode,scndnode,thirdnode;
INT           firstmatch,scndmatch,thirdmatch;
INT           surfmatch;
INT           nodeId;

#ifdef DEBUG 
dstrc_enter("inpdesign_dsurf_fesurf");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   nnodeonsurf = ndsurf_fenode[i];
   /* loop all GSURFS and check whether there nodes are on this actdsurf */
   for (j=0; j<genprob.ngsurf; j++)
   {
      surfmatch=0;
      firstmatch=scndmatch=thirdmatch=0;
      actgsurf = genprob.gsurfs[j];
      if (actgsurf->dsurf) continue;
      /*----------------------- check whether firstnode is on actdsurf */
      firstnode = actgsurf->gnode[0]->node->Id;
      for (k=0; k<nnodeonsurf; k++)
      {
         if (dsurf_fenode[i][k]==firstnode)
         {
            firstmatch=1;
            break;
         }
      }
      /*------------------------ if firstnode not on actdsurf continue */
      if (!firstmatch) continue;
      /*------------------------------- check for scndnode on actdsurf */
      scndnode = actgsurf->gnode[1]->node->Id;
      for (k=0; k<nnodeonsurf; k++)
      {
         if (dsurf_fenode[i][k]==scndnode)
         {
            scndmatch=1;
            break;
         }
      }
      /*------------------------ if scndnode not on actdsurf continue */
      if (!scndmatch) continue;
      /*----------------------------- check for thirdnode on actdsurf */
      thirdnode = actgsurf->gnode[2]->node->Id;
      for (k=0; k<nnodeonsurf; k++)
      {
         if (dsurf_fenode[i][k]==thirdnode)
         {
            thirdmatch=1;
            break;
         }
      }
      /*----------------------- if thirdnode not on actdsurf continue */
      if (!thirdmatch) continue;
      /* three nodes of actgsurf are on actdsurf, so actgsurf is in 
         actdsurf
      */
      actgsurf->dsurf = actdsurf;
   } /* end loop j over gsurfs */    
}/* end loop i over dsurfs */
/*- loop dsurfaces again and set pointers GNODE -> DSURF for all GNODEs */
/*------------- which do NOT already have a pointer to a DNODE or DLINE */

  for (i=0; i<design->ndsurf; i++)
  {
    actdsurf = &(design->dsurf[i]);
    nnodeonsurf = ndsurf_fenode[i];
    for (j=0; j<nnodeonsurf; j++)
    {
      nodeId = dsurf_fenode[i][j];
      actgnode = genprob.gnodes[nodeId];
      if (actgnode == NULL) continue;
      if (actgnode->ondesigntyp != ondnothing) continue;
      actgnode->ondesigntyp = ondsurf;
      actgnode->d.dsurf = actdsurf;
    }/* end loop j over GNODEs */
  } /* end loop i over DSURF */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dsurf_fesurf */



/*----------------------------------------------------------------------*
 | make topology DSURF <-> GSURF                          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dvol_fevol()
{
INT           i,j,k;
DVOL         *actdvol;
GVOL         *actgvol;
GNODE        *actgnode;
INT           ngnode;
ELEMENT      *actele;
INT           firstnode  = 0;
INT           scndnode   = 0;
INT           thirdnode  = 0;
INT           fourthnode = 0;
INT           firstmatch,scndmatch,thirdmatch,fourthmatch;
INT           nodeId;

#ifdef DEBUG 
dstrc_enter("inpdesign_dvol_fevol");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   for (j=0; j<genprob.ngvol; j++)
   {
      actgvol = genprob.gvols[j];
      actele  = actgvol->element;
      firstmatch=scndmatch=thirdmatch=fourthmatch=0;
      /*-------- find four nodes of element that are not in one surface */
      switch (actele->distyp)
      {
      case hex8:
      case hex20: 
      case hex27: 
         firstnode  = actele->node[0]->Id;
         scndnode   = actele->node[2]->Id;
         thirdnode  = actele->node[4]->Id;
         fourthnode = actele->node[6]->Id;
      break;
      case tet4:  
      case tet10: 
         firstnode  = actele->node[0]->Id;
         scndnode   = actele->node[1]->Id;
         thirdnode  = actele->node[2]->Id;
         fourthnode = actele->node[3]->Id;
      break;
      default: dserror("Unknown type of discretization"); break;
      }
      /*------------------------ check whether first node is in actdvol */
      for (k=0; k<ndvol_fenode[i]; k++)
      {
         if (firstnode==dvol_fenode[i][k])
         {
            firstmatch=1;
            break;
         }
      }
      /*--------------------------- if firstnode not in volume continue */
      if (!firstmatch) continue;
      /*------------------------- check whether scnd node is in actdvol */
      for (k=0; k<ndvol_fenode[i]; k++)
      {
         if (scndnode==dvol_fenode[i][k])
         {
            scndmatch=1;
            break;
         }
      }
      /*---------------------------- if scndnode not in volume continue */
      if (!scndmatch) continue;
      /*------------------------ check whether third node is in actdvol */
      for (k=0; k<ndvol_fenode[i]; k++)
      {
         if (thirdnode==dvol_fenode[i][k])
         {
            thirdmatch=1;
            break;
         }
      }
      /*---------------------------- if thirnode not in volume continue */
      if (!thirdmatch) continue;
      /*------------------------ check whether fourth node is in actdvol */
      for (k=0; k<ndvol_fenode[i]; k++)
      {
         if (fourthnode==dvol_fenode[i][k])
         {
            fourthmatch=1;
            break;
         }
      }
      /*---------------------------- if fourthnode not in volume continue */
      if (!fourthmatch) continue;
      /*---------------------------------- Yes, this acgvol is in actdvol */
      actgvol->dvol = actdvol;
   }/* end loop j over gvols */
}/* end loop i over dvols */

/* loop all design volumes again and make pointers GNODE -> DVOL in all */
/* GNODEs which do not already have a pointer to a DSURF,DLINE or DNODE */

  for (i=0; i<design->ndvol; i++)
  {
    actdvol = &(design->dvol[i]);
    ngnode  = ndvol_fenode[i];
    for (j=0; j<ngnode; j++)
    {
      nodeId = dvol_fenode[i][j];
      actgnode = genprob.gnodes[nodeId];
      if (actgnode == NULL) continue;
      if (actgnode->ondesigntyp != ondnothing) continue;
      actgnode->ondesigntyp = ondvol;
      actgnode->d.dvol = actdvol;
    }/* end loop j over GNODEs */
  } /* end loop i over DVOL */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_dvol_fevol */


