/*!----------------------------------------------------------------------
\file
\brief contains history relevant routines of the shell contact

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef S8CONTACT
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "s8contact.h"
#include "shell8.h"



/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*!
\addtogroup CONTACTS8
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
extern struct _SHELLCONTACT shellcontact;



/*!---------------------------------------------------------------------
\brief write restart for contact variables

<pre>                                                        m.gee 3/03
</pre>
\param  actintra        *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_restartwrite(INTRA *actintra,INT step)
{
INT              i,l;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
INT              ierr;
FILE            *out;
long int         handle;
RESTARTCONTACT   res;
SHELLNODERES     resnode;
char             resname[100];
#ifdef DEBUG
dstrc_enter("s8_contact_restartwrite");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
res.numnp   = numnp;
out         = allfiles.out_pss;
/*----------------------------------------------------------------------*/
res.handles = (long int*)CCACALLOC(numnp,sizeof(long int));
/*------------------------ loop nodes and write the shellnode structure */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   for (l=0; l<6; l++)
   resnode.xc_his[l]  = actcnode->xc_his[l];
   resnode.topflag    = actcnode->topflag;
   resnode.topproj    = actcnode->topproj;
   resnode.histopflag = actcnode->histopflag;
   resnode.histopproj = actcnode->histopproj;
   resnode.top_ln     = actcnode->top_ln;
   resnode.top_lt[0]  = actcnode->top_lt[0];
   resnode.top_lt[1]  = actcnode->top_lt[1];

   resnode.botflag    = actcnode->botflag;
   resnode.botproj    = actcnode->botproj;
   resnode.hisbotflag = actcnode->hisbotflag;
   resnode.hisbotproj = actcnode->hisbotproj;
   resnode.bot_ln     = actcnode->bot_ln;
   resnode.bot_lt[0]  = actcnode->bot_lt[0];
   resnode.bot_lt[1]  = actcnode->bot_lt[1];

   pss_write("cnode",1,1,sizeof(SHELLNODERES),&resnode,&(res.handles[i]),out,&ierr);
}
/*----------------------------------------------- now write the handles */
pss_write("mainhandle",numnp,1,sizeof(long int),res.handles,&(res.mainhandle),out,&ierr);
/*---------------------------- and finally the restart structure itself */
sprintf(resname,"cs8_%d",step);
pss_write(resname,1,1,sizeof(RESTARTCONTACT),&res,&handle,out,&ierr);
/*-------------------------------------------- free the allocated space */
CCAFREE(res.handles);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_restartwrite */



/*!---------------------------------------------------------------------
\brief read restart for contact variables

<pre>                                                        m.gee 3/03
</pre>
\param  actintra        *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_restartread(INTRA *actintra,INT step)
{
INT              i,j,k,l;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
SHELLNODERES     tmpnode;
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
INT              ierr;
FILE            *in;
long int         handle;
RESTARTCONTACT   res;
char             resname[100];
#ifdef DEBUG
dstrc_enter("s8_contact_restartread");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
in          = allfiles.in_pss;
/*---------------------------------- find the contact restart structure */
sprintf(resname,"cs8_%d",step);
pss_chck(resname,&handle,in,&ierr);
if (ierr != 1) dserror("Cannot restart contact, step doesn't exist");
/*---------------------------------- read the contact restart structure */
pss_read_name_handle(resname,&i,&i,&j,&res,&handle,in,&ierr);
/*------------------------------------------ read the vector of handles */
/*----------------------------------------------------------------------*/
res.handles = (long int*)CCACALLOC(res.numnp,sizeof(long int));
/*---------------------------------------------- read vector of handles */
i = 1;
j = sizeof(long int);
pss_read_name_handle("mainhandle",&(res.numnp),&i,&j,res.handles,&(res.mainhandle),in,&ierr);
/*-------------------------------------------- loop nodes and read them */
j = 1;
k = sizeof(SHELLNODERES);
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   pss_read_name_handle("cnode",&j,&j,&k,&tmpnode,&(res.handles[i]),in,&ierr);
   if (ierr != 1) dserror("Failed to read contact nodes");
   /* copy data to the node */
   for (l=0; l<6; l++)
   actcnode->xc_his[l]  = tmpnode.xc_his[l];

   actcnode->topflag    = tmpnode.topflag;
   actcnode->topproj    = tmpnode.topproj;
   actcnode->histopflag = tmpnode.histopflag;
   actcnode->histopproj = tmpnode.histopproj;
   actcnode->top_ln     = tmpnode.top_ln;
   actcnode->top_lt[0]  = tmpnode.top_lt[0];
   actcnode->top_lt[1]  = tmpnode.top_lt[1];

   actcnode->botflag    = tmpnode.botflag;
   actcnode->botproj    = tmpnode.botproj;
   actcnode->hisbotflag = tmpnode.hisbotflag;
   actcnode->hisbotproj = tmpnode.hisbotproj;
   actcnode->bot_ln     = tmpnode.bot_ln;
   actcnode->bot_lt[0]  = tmpnode.bot_lt[0];
   actcnode->bot_lt[1]  = tmpnode.bot_lt[1];
}
/*-------------------------------------------- free the allocated space */
CCAFREE(res.handles);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_restartread */

/*!---------------------------------------------------------------------
\brief set lagrange multipliers for contact to zero

<pre>                                                        m.gee 3/03
</pre>
\param  actfield       *FIELD        (i) the field
\param  actpart        *PARTITION    (i) mypartition
\param  actintra       *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_setlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra)
{
INT              i;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
#ifdef DEBUG
dstrc_enter("s8_contact_setlagr");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
/*--------------------------------------------------- loop all my nodes */
for (i=0; i<numnp; i++)
{
   if (cnode[i].node->proc != myrank) continue;
   actcnode = &(cnode[i]);
   if (actcnode->topflag == s8_c_off)
   {
      actcnode->top_ln    = 0.0;
      actcnode->top_lt[0] = 0.0;
      actcnode->top_lt[1] = 0.0;
   }
   if (actcnode->botflag == s8_c_off)
   {
      actcnode->bot_ln    = 0.0;
      actcnode->bot_lt[0] = 0.0;
      actcnode->bot_lt[1] = 0.0;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_setlagr */



/*!---------------------------------------------------------------------
\brief update lagrange multipliers for contact

<pre>                                                        m.gee 3/03
</pre>
\param  actfield       *FIELD        (i) the field
\param  actpart        *PARTITION    (i) mypartition
\param  actintra       *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_updlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra)
{
INT              i;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
#ifdef DEBUG
dstrc_enter("s8_contact_updlagr");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
/*--------------------------------------------------- loop all my nodes */
for (i=0; i<numnp; i++)
{
   if (cnode[i].node->proc != myrank) continue;
   actcnode = &(cnode[i]);
   /*----------------------------------- check for top contact to be on */
   if (actcnode->topflag == s8_c_on)
   {
      /*------- update the lagrange multiplier for normal contact force */
      actcnode->top_ln    = actcnode->top_tn;
      actcnode->top_lt[0] = actcnode->top_tT[0];
      actcnode->top_lt[1] = actcnode->top_tT[1];
   }
   /*----------------------------------- check for bot contact to be on */
   if (actcnode->botflag == s8_c_on)
   {
      /*------- update the lagrange multiplier for normal contact force */
      actcnode->bot_ln    = actcnode->bot_tn;
      actcnode->bot_lt[0] = actcnode->bot_tT[0];
      actcnode->bot_lt[1] = actcnode->bot_tT[1];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_updlagr */
/*!---------------------------------------------------------------------
\brief store converged contact variables

<pre>                                                        m.gee 2/03
</pre>
\param  actintra        *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_history(INTRA *actintra)
{
INT              i,j;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
#ifdef DEBUG
dstrc_enter("s8_contact_history");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
/*--------------------------------------- loop nodes which belong to me */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   if (actcnode->node->proc != myrank) continue;
   /* this is converged state, so store values */
   /* top */
   actcnode->histopflag    = actcnode->topflag;
   actcnode->histopproj    = actcnode->topproj;
   actcnode->histopele     = actcnode->topele;
   actcnode->oldhisxitop[0]= actcnode->hisxitop[0];
   actcnode->oldhisxitop[1]= actcnode->hisxitop[1];
   actcnode->hisxitop[0]   = actcnode->xitop[0];
   actcnode->hisxitop[1]   = actcnode->xitop[1];
   actcnode->histopgap     = actcnode->topgap;
   actcnode->his_top_ln    = actcnode->top_ln;
   actcnode->his_top_lt[0] = actcnode->top_lt[0];
   actcnode->his_top_lt[1] = actcnode->top_lt[1];
   /* bot */
   actcnode->hisbotflag    = actcnode->botflag;
   actcnode->hisbotproj    = actcnode->botproj;
   actcnode->hisbotele     = actcnode->botele;
   actcnode->oldhisxibot[0]= actcnode->hisxibot[0];
   actcnode->oldhisxibot[1]= actcnode->hisxibot[1];
   actcnode->hisxibot[0]   = actcnode->xibot[0];
   actcnode->hisxibot[1]   = actcnode->xibot[1];
   actcnode->hisbotgap     = actcnode->botgap;
   actcnode->his_bot_ln    = actcnode->bot_ln;
   actcnode->his_bot_lt[0] = actcnode->bot_lt[0];
   actcnode->his_bot_lt[1] = actcnode->bot_lt[1];
   /*
    current confiuration of this converged step is needed,
    when the node is sliding from one element to another
   */
   for (j=0; j<6; j++)
   actcnode->xc_his[j] = actcnode->xr[j] + actcnode->node->sol.a.da[0][j];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_history */


/*!---------------------------------------------------------------------
\brief set contact variables back to last converged step

<pre>                                                        m.gee 2/03
</pre>
\param  actintra        *INTRA        (i) the intracommunicator
\return void

------------------------------------------------------------------------*/
void s8_contact_historyback(INTRA *actintra)
{
INT              i;
SHELLNODE       *actcnode,*cnode;                  /* the active contact node */
INT              myrank,nproc;                     /* parallel stuff */
INT              numnp;                            /* number of slave nodes (ususally all nodes) */
#ifdef DEBUG
dstrc_enter("s8_contact_historyback");
#endif
/*----------------------------------------------------------------------*/
numnp       = shellcontact.numnp;
cnode       = shellcontact.cnode;
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
/*--------------------------------------- loop nodes which belong to me */
for (i=0; i<numnp; i++)
{
   actcnode = &(cnode[i]);
   if (actcnode->node->proc != myrank) continue;
   /* top */
   actcnode->topflag     =  actcnode->histopflag     ;
   actcnode->topproj     =  actcnode->histopproj     ;
   actcnode->topele      =  actcnode->histopele      ;
   actcnode->xitop[0]    =  actcnode->hisxitop[0]    ;
   actcnode->xitop[1]    =  actcnode->hisxitop[1]    ;
   actcnode->hisxitop[0] =  actcnode->oldhisxitop[0] ;
   actcnode->hisxitop[1] =  actcnode->oldhisxitop[1] ;
   actcnode->topgap      =  actcnode->histopgap      ;
   actcnode->top_ln      =  actcnode->his_top_ln     ;
   actcnode->top_lt[0]   =  actcnode->his_top_lt[0]  ;
   actcnode->top_lt[1]   =  actcnode->his_top_lt[1]  ;
   /* bot */
   actcnode->botflag     =  actcnode->hisbotflag     ;
   actcnode->botproj     =  actcnode->hisbotproj     ;
   actcnode->botele      =  actcnode->hisbotele      ;
   actcnode->xibot[0]    =  actcnode->hisxibot[0]    ;
   actcnode->xibot[1]    =  actcnode->hisxibot[1]    ;
   actcnode->hisxibot[0] =  actcnode->oldhisxibot[0] ;
   actcnode->hisxibot[1] =  actcnode->oldhisxibot[1] ;
   actcnode->botgap      =  actcnode->hisbotgap      ;
   actcnode->bot_ln      =  actcnode->his_bot_ln     ;
   actcnode->bot_lt[0]   =  actcnode->his_bot_lt[0]  ;
   actcnode->bot_lt[1]   =  actcnode->his_bot_lt[1]  ;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_contact_historyback */



/*! @} (documentation module close)*/
#endif



#endif
