/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                        genk 04/04    |
 | number of local co-ordinate systems                                  |
 | vector of structures of local co-ordinate systems                    |
 | defined in input_locsys.c                                          |
 *----------------------------------------------------------------------*/
extern INT            numlocsys;
extern struct _LOCSYS *locsys;
/*!---------------------------------------------------------------------
\brief calculate normal of fluid nodes

<pre>                                                         genk 04/02

   init = 0    just calculate the normal
   init = 1    initialise and calculate the normal and old normal
   init = 2    copy first actn to oldn and then calculate new normal

</pre>
\return void

------------------------------------------------------------------------*/
void fluid_cal_normal(FIELD *actfield, INT init,
                      CALC_ACTION *action         )
{
INT         i,j,l;       /* simply some counters */
INT         foundit;       /* found flag */
INT         ilocsys;       /* index acutal locsys */
static INT  numele_total;  /* total number of elements in actual dis */
static INT  numnp_total;   /* total number of nodes in actual dis */
static INT  numnorm;       /* number of nodes normals have to calculated for */
DOUBLE      norm;          /* length of normal */
ELEMENT    *actele;        /* the actual element */
NODE       *actnode;       /* the actual node */
#ifdef D_FSI
INT         ngnode,k;        /* number of gnodes */
INT         ngline;        /* number of glines */
GLINE      *actgline;      /* the actual gline */
#endif
LOCSYS     *actlocsys;     /* the actual local co-system xzy* */
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("fluid_cal_normal");
#endif

fdyn = alldyn[genprob.numff].fdyn;
numele_total=actfield->dis[0].numele;
numnp_total=actfield->dis[0].numnp;

if (init==1) /* initialise and calculate initial normal */
{
   numnorm=0;
   numele_total=actfield->dis[0].numele;
   numnp_total=actfield->dis[0].numnp;
   if (numlocsys==0 && fdyn->freesurf<6) goto end;
   /* loop elements and allocate normals at the nodes */
   for (i=0;i<numele_total;i++)
   {
      actele=&(actfield->dis[0].element[i]);
      dsassert(actele->eltyp==el_fluid2,
      "normals only for fluid2 element implemented!\n");
      /* we need normals for
         local co-sys and for a generalised free surface */
      /* check for local co-sys */
      if (actele->locsys==locsys_yes) /* locsys */
      {
         for (j=0;j<actele->numnp;j++)
         {
            actnode=actele->node[j];
            if (actnode->actn!=NULL) continue;
            ilocsys=actnode->locsysId-1;
            if (ilocsys==-1) continue;
            actlocsys=&(locsys[ilocsys]);
            if (actlocsys->locsystyp!=locsys_fmc) continue;
            actnode->actn=(DOUBLE*)CCACALLOC(3,sizeof(DOUBLE));
            for(l=0;l<3;l++) actnode->actn[l]=ZERO;
            numnorm++;
         }
      }
#ifdef D_FSI
#ifdef D_FLUID2
      if (actele->e.f2->fs_on==6) /* generalised free surface */
      {
         ngline=actele->g.gsurf->ngline;
         for (j=0; j<ngline; j++)
         {
            actgline = actele->g.gsurf->gline[j];
            if (actgline->freesurf==NULL) continue;
            ngnode = actgline->ngnode;
            for (k=0;k<ngnode;k++)
            {
               actnode=actgline->gnode[k]->node;
               if (actnode->actn!=NULL) continue;
               actnode->actn=(DOUBLE*)CCACALLOC(3,sizeof(DOUBLE));
               actnode->oldn=(DOUBLE*)CCACALLOC(3,sizeof(DOUBLE));
               for(l=0;l<3;l++) actnode->actn[l]=ZERO;
               numnorm++;
            }
         }
      }
#endif
#endif
   }
}


if (numnorm>0) /* calculate nodal normals */
{
   if (init==2) /* copy actn to oldn */
   {
      for (i=0;i<numnp_total;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
         if (actnode->oldn==NULL) continue;
         for(l=0;l<3;l++)
            actnode->oldn[l]=actnode->actn[l];
      }
   }
   *action=calc_fluid_normal;
   for (i=0;i<numele_total;i++)
   {
      actele=&(actfield->dis[0].element[i]);
      dsassert(actele->eltyp==el_fluid2,
      "normals only for fluid2 element implemented!\n");
      foundit=0;
      for (j=0;j<actele->numnp;j++)
      {
         actnode=actele->node[j];
         if (actnode->actn!=NULL)
         {
            foundit++;
            break;
         }
      }
      if (foundit>0)
      {
#ifdef D_FLUID2
         fluid2(NULL,NULL,actele,NULL,NULL,NULL,NULL,NULL,NULL,
                action,NULL,NULL,NULL);
#else
         dserror("Fluid2 not compiled in!\n");
#endif
      }
   }
   /*------------------------------normalise actn and copy actn to oldn */
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[0].node[i]);
      if (actnode->actn==NULL) continue;
      norm=sqrt(DSQR(actnode->actn[0])+DSQR(actnode->actn[1])
               +DSQR(actnode->actn[2]));
      for (l=0;l<3;l++)
         actnode->actn[l]/=norm;
      /*------------------------------------------ output to the screen */
#ifdef DEUBG
      if (init==2)
      printf("Id = %5d   n1 = % 10.8lf   n2 = % 10.8lf \n",
         actnode->Id, actnode->actn[0],actnode->actn[1]);
#endif
   }
   if (init==1) /* copy actn to oldn */
   {
      for (i=0;i<numnp_total;i++)
      {
         actnode=&(actfield->dis[0].node[i]);
         if (actnode->oldn==NULL) continue;
         for(l=0;l<3;l++)
            actnode->oldn[l]=actnode->actn[l];
      }
   }

}


/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_cal_normal */

#endif
/*! @} (documentation module close)*/
