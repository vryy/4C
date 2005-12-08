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
static void inpdesign_dpoint_fenode_read(
    INT        **dnode_fenode,
    INT         *ndnode_fenode
    );

static void inpdesign_dline_fenode_read(
    INT        **dline_fenode,
    INT         *ndline_fenode
    );

static void inpdesign_dsurf_fenode_read(
    INT         **dsurf_fenode,
    INT          *ndsurf_fenode
    );

static void inpdesign_dvol_fenode_read(
    INT           **dvol_fenode,
    INT            *ndvol_fenode
    );


static void inpdesign_dpoint_fenode(
    DISCRET     *actdis,
    INT        **act_dnode_fenode,
    INT         *ndnode_fenode
    );

static void inpdesign_dline_feline(
    DISCRET     *actdis,
    INT        **act_dline_fenode,
    INT         *ndline_fenode
    );


static void inpdesign_dsurf_fesurf(
    DISCRET      *actdis,
    INT         **act_dsurf_fenode,
    INT          *ndsurf_fenode
    );


static void inpdesign_dvol_fevol(
    DISCRET        *actdis,
    INT           **act_dvol_fenode,
    INT            *ndvol_fenode
    );

/*----------------------------------------------------------------------*
 | global variables in this file in this file                      3/02 |
 *----------------------------------------------------------------------*/
#if 0
static INT *ndnode_fenode;
static INT **dnode_fenode;
static INT **dnode_fenode2;

static INT *ndline_fenode;
static INT **dline_fenode;
static INT **dline_fenode2;

static INT *ndsurf_fenode;
static INT **dsurf_fenode;
static INT **dsurf_fenode2;

static INT *ndvol_fenode;
static INT **dvol_fenode;
static INT **dvol_fenode2;
#endif

static INT  *gnode_ind;

static INT   init = 0;

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
void inpdesign_topology_fe(
    DISCRET        *actdis
    )
{

  INT i,j,k,l,id;


#ifdef DEBUG
  dstrc_enter("inpdesign_topology_fe");
#endif


  if (init != 1)
  {

    /* allocate arrays for the fe-design info */
    design->ndnode_fenode = (INT*)CCACALLOC(design->ndnode,sizeof(INT));
    design->dnode_fenode  = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));

    design->ndline_fenode = (INT*)CCACALLOC(design->ndline,sizeof(INT));;
    design->dline_fenode  = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;

    design->ndsurf_fenode = (INT*)CCACALLOC(design->ndsurf,sizeof(INT));;
    design->dsurf_fenode  = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;

    design->ndvol_fenode = (INT*)CCACALLOC(design->ndvol,sizeof(INT));;
    design->dvol_fenode  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;


    /* read the fe-nodes on each design object */
    inpdesign_dpoint_fenode_read(design->dnode_fenode,design->ndnode_fenode);
    inpdesign_dline_fenode_read(design->dline_fenode,design->ndline_fenode);
    inpdesign_dsurf_fenode_read(design->dsurf_fenode,design->ndsurf_fenode);
    inpdesign_dvol_fenode_read(design->dvol_fenode,design->ndvol_fenode);



    /* create topology for artifical dis */
    if (genprob.create_dis == 1 || genprob.create_ale == 1)
    {
      design->dnode_fenode2 = (INT**)CCACALLOC(design->ndnode,sizeof(INT*));
      design->dline_fenode2 = (INT**)CCACALLOC(design->ndline,sizeof(INT*));;
      design->dsurf_fenode2 = (INT**)CCACALLOC(design->ndsurf,sizeof(INT*));;
      design->dvol_fenode2  = (INT**)CCACALLOC(design->ndvol,sizeof(INT*));;

      for (i=0; i<design->ndnode; i++)
      {
        design->dnode_fenode2[i]= (INT*)CCAMALLOC(design->ndnode_fenode[i]*sizeof(INT));
        for ( j=0; j<design->ndnode_fenode[i]; j++)
          design->dnode_fenode2[i][j] = design->dnode_fenode[i][j] + genprob.nodeshift;
      }

      for (i=0; i<design->ndline; i++)
      {
        design->dline_fenode2[i]= (INT*)CCAMALLOC(design->ndline_fenode[i]*sizeof(INT));
        for ( j=0; j<design->ndline_fenode[i]; j++)
          design->dline_fenode2[i][j] = design->dline_fenode[i][j] + genprob.nodeshift;
      }

      for (i=0; i<design->ndsurf; i++)
      {
        design->dsurf_fenode2[i]= (INT*)CCAMALLOC(design->ndsurf_fenode[i]*sizeof(INT));
        for ( j=0; j<design->ndsurf_fenode[i]; j++)
          design->dsurf_fenode2[i][j] = design->dsurf_fenode[i][j] + genprob.nodeshift;
      }

      for (i=0; i<design->ndvol; i++)
      {
        design->dvol_fenode2[i] = (INT*)CCAMALLOC(design->ndvol_fenode[i]*sizeof(INT));
        for ( j=0; j<design->ndvol_fenode[i]; j++)
          design->dvol_fenode2[i][j]  = design->dvol_fenode[i][j]  + genprob.nodeshift;
      }
    }



    init = 1;
  }


  /* make the topology between DESIGN and GEOMETRY */
  /*===============================================*/

  gnode_ind = (INT*)CCACALLOC(genprob.maxnode,sizeof(INT));


  /* initialize gnode_ind */
  for(l=0;l<genprob.maxnode;l++)
    gnode_ind[l] = -1;


  /* loop the gnodes in this dis */
  for (k=0; k<actdis->ngnode; k++)
  {
    /* fill gnode_ind to find the nodes */
    id = actdis->gnode[k].node->Id;
    dsassert(id < genprob.maxnode,"Zu wenig KNOTEN");

    gnode_ind[id] = k;
  }


  switch (actdis->disclass)
  {
    /* the normal case */
    case dc_normal:
    case dc_subdiv_io:
      inpdesign_dpoint_fenode(actdis,design->dnode_fenode,design->ndnode_fenode);
      inpdesign_dline_feline(actdis,design->dline_fenode,design->ndline_fenode);
      inpdesign_dsurf_fesurf(actdis,design->dsurf_fenode,design->ndsurf_fenode);
      inpdesign_dvol_fevol(actdis,design->dvol_fenode,design->ndvol_fenode);
      break;

    /* create ale field */
    case dc_created_ale:
    case dc_subdiv_io_created_ale:
    case dc_created_f2p:
    case dc_created_tu:
      inpdesign_dpoint_fenode(actdis,design->dnode_fenode2,design->ndnode_fenode);
      inpdesign_dline_feline(actdis,design->dline_fenode2,design->ndline_fenode);
      inpdesign_dsurf_fesurf(actdis,design->dsurf_fenode2,design->ndsurf_fenode);
      inpdesign_dvol_fevol(actdis,design->dvol_fenode2,design->ndvol_fenode);
      break;

    /* create discretization by subdivision */
    case dc_subdiv_calc:
      inpdesign_dpoint_fenode(actdis,design->dnode_fenode2,design->ndnode_fenode);
      inpdesign_dline_feline(actdis,design->dline_fenode2,design->ndline_fenode);
      inpdesign_dsurf_fesurf(actdis,design->dsurf_fenode2,design->ndsurf_fenode);
      inpdesign_dvol_fevol(actdis,design->dvol_fenode2,design->ndvol_fenode);
      break;

    default:
      dserror("Unknown disclass!!");
      break;
  }






  CCAFREE(gnode_ind);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of inpdesign_topology_fe */







/*----------------------------------------------------------------------*
 | input of design nodes  to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dpoint_fenode_read(
    INT        **dnode_fenode,
    INT         *ndnode_fenode
    )
{
INT    i,ierr;
INT    dnode;
INT    found = 0;
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
         found = 1; /*fount at least one DNODE-NODE TOPOLOGY*/
         if (dnode==i+1)
         {
            ndnode_fenode[i]=1;
            dnode_fenode[i] = (INT*)CCAMALLOC(ndnode_fenode[i]*sizeof(INT));
            frint("NODE",&(design->dnode_fenode[i][0]),&ierr);
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
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpdesign_dpoint_fenode_read */


/*----------------------------------------------------------------------*
 | input of design line  to fe-node topology             m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dline_fenode_read(
    INT        **dline_fenode,
    INT         *ndline_fenode
    )
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
static void inpdesign_dsurf_fenode_read(
    INT         **dsurf_fenode,
    INT          *ndsurf_fenode
    )
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
static void inpdesign_dvol_fenode_read(
    INT           **dvol_fenode,
    INT            *ndvol_fenode
    )
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
static void inpdesign_dpoint_fenode(
    DISCRET     *actdis,
    INT        **act_dnode_fenode,
    INT         *ndnode_fenode
    )
{

  INT           i;
  INT           nodeId;
  DNODE        *actdnode;
  GNODE        *actgnode;


#ifdef DEBUG
  dstrc_enter("inpdesign_dpoint_fenode");
#endif


  for (i=0; i<design->ndnode; i++)
  {
    actdnode = &(design->dnode[i]);

    if (act_dnode_fenode[i] == NULL)
      dserror("DNODE without FE node. Uncollapsed nodes in GiD?");

    if (design->ndnode_fenode[i] == 0 )
      continue;
      /*dserror("Only one fe_node possible on each dnode !!");*/

      nodeId   = gnode_ind[act_dnode_fenode[i][0]];

    if (nodeId == -1 )
      continue;
      /*dserror("MIXUP GNODE dnode: %d; fenode; %d; nodeid: %d"
          ,design->dnode[i].Id,act_dnode_fenode[i][0],nodeId);*/


    actgnode = &(actdis->gnode[nodeId]);
    actgnode->ondesigntyp = ondnode;
    actgnode->d.dnode = actdnode;
  }


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of inpdesign_dpoint_fenode */



/*----------------------------------------------------------------------*
 | make topology DLINE <-> GLINE                          m.gee 3/02    |
 *----------------------------------------------------------------------*/
static void inpdesign_dline_feline(
    DISCRET     *actdis,
    INT        **act_dline_fenode,
    INT         *ndline_fenode
    )
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
   for (j=0; j<actdis->ngline; j++)
   {
      linematch=0;
      firstmatch=scndmatch=thirdmatch=0;
      actgline = &(actdis->gline[j]);
      if (actgline->dline) continue;
      /*----------------------- check whether firstnode is on actdline */
      firstnode = actgline->gnode[0]->node->Id;
      for (k=0; k<nnodeonline; k++)
      {
         if (act_dline_fenode[i][k]==firstnode)
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
         if (act_dline_fenode[i][k]==scndnode)
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
            if (act_dline_fenode[i][k]==thirdnode)
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
      nodeId = gnode_ind[act_dline_fenode[i][j]];
      if (nodeId == -1 )
        continue;
        /*dserror("MIXUP GLINE");*/

      actgnode = &(actdis->gnode[nodeId]);
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
static void inpdesign_dsurf_fesurf(
    DISCRET      *actdis,
    INT         **act_dsurf_fenode,
    INT          *ndsurf_fenode
    )
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
   for (j=0; j<actdis->ngsurf; j++)
   {
      surfmatch=0;
      firstmatch=scndmatch=thirdmatch=0;
      actgsurf = &(actdis->gsurf[j]);
      if (actgsurf->dsurf) continue;
      /*----------------------- check whether firstnode is on actdsurf */
      firstnode = actgsurf->gnode[0]->node->Id;
      for (k=0; k<nnodeonsurf; k++)
      {
         if (act_dsurf_fenode[i][k]==firstnode)
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
         if (act_dsurf_fenode[i][k]==scndnode)
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
         if (act_dsurf_fenode[i][k]==thirdnode)
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
      nodeId = gnode_ind[act_dsurf_fenode[i][j]];

      if (nodeId == -1 )
        continue;
        /*dserror("MIXUP GSURF");*/

      actgnode = &(actdis->gnode[nodeId]);
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
static void inpdesign_dvol_fevol(
    DISCRET        *actdis,
    INT           **act_dvol_fenode,
    INT            *ndvol_fenode
    )
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
   for (j=0; j<actdis->ngvol; j++)
   {
      actgvol = &(actdis->gvol[j]);
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
         if (firstnode==act_dvol_fenode[i][k])
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
         if (scndnode==act_dvol_fenode[i][k])
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
         if (thirdnode==act_dvol_fenode[i][k])
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
         if (fourthnode==act_dvol_fenode[i][k])
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
      nodeId = gnode_ind[act_dvol_fenode[i][j]];

      if (nodeId == -1 )
        continue;
        /*dserror("MIXUP GVOL");*/

      actgnode = &(actdis->gnode[nodeId]);
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


