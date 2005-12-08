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
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 |                                                        genk 04/04    |
 | number of local co-ordinate systems                                  |
 | vector of structures of local co-ordinate systems                    |
 | defined in input_locsys.c                                          |
 *----------------------------------------------------------------------*/
extern INT            numlocsys;
extern struct _LOCSYS *locsys;
/*!---------------------------------------------------------------------
\brief define locsys for slip BC along curved DLINEs

<pre>                                                         genk 04/02

</pre>
\return void

------------------------------------------------------------------------*/
void fluid_locsys(
    FIELD              *actfield,
    INT                 disnum,
    FLUID_DYNAMIC      *fdyn
    )
{
INT i,j;
INT numnp_total;
INT ilocsys;
DOUBLE   *xloc, *yloc, *zloc;
DOUBLE    norm;
NODE     *actnode;
LOCSYS   *actlocsys;


#ifdef DEBUG
dstrc_enter("fluid_locsys");
#endif

numnp_total=actfield->dis[disnum].numnp;

for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[disnum].node[i]);
   /*------------------------------------------- check for local co-sys */
   ilocsys=actnode->locsysId-1;
   if (ilocsys==-1) continue;
   actlocsys=&(locsys[ilocsys]);
   if (actlocsys->locsystyp!=locsys_fmc) continue;
   dsassert(actnode->actn!=NULL,"locsys_fmc but no normal defined at node!\n");
   /*--------------- node has individual local co-sys based on the mass
      consistent normal
   */
   numlocsys++;
   locsys=(LOCSYS*)CCAREALLOC(locsys,numlocsys*sizeof(LOCSYS));
   actlocsys=&(locsys[numlocsys-1]);
   actlocsys->Id=numlocsys;
   actnode->locsysId=numlocsys;
   /*-------------------------------------------- allocate base vectors */
   xloc=amdef("xloc",&actlocsys->xloc,3,1,"DV");
   yloc=amdef("yloc",&actlocsys->yloc,3,1,"DV");
   zloc=amdef("zloc",&actlocsys->zloc,3,1,"DV");
   actlocsys->locsystyp = locsys_fmc;
   if (fdyn->numdf==3)
   {
      /* 2D-case is straight forward: normal = local y-axis */
      yloc[0]=actnode->actn[0];
      yloc[1]=actnode->actn[1];
      yloc[2]=actnode->actn[2];
      dsassert(FABS(yloc[2])<EPS15,"3D base vector y in 2D case!\n");
      /* define xloc orthogonal in xy-plane to yloc */
      xloc[0] = -yloc[1];
      xloc[1] =  yloc[0];
      xloc[2] = ZERO;
      /* vector product to check zloc */
      zloc[0] = xloc[1]*yloc[2] - xloc[2]*yloc[1];
      zloc[1] = xloc[2]*yloc[0] - xloc[0]*yloc[2];
      zloc[2] = xloc[0]*yloc[1] - xloc[1]*yloc[0];
      if (zloc[2]<ZERO) /* change sign of yloc */
      {
         xloc[0] *=-ONE;
         xloc[1] *=-ONE;
         zloc[2] *=-ONE;
      }
   }
   else if (fdyn->numdf==4) /* 3D how to choose the tangents??? */
      dserror("fluid locsys for 3D case not implemented yet!\n");
   else
      dserror("fluid locsys for 4D and higer dim. cases not implemented yet!\n");

   /*---------------------------------------- normalise the base vectors */
   norm = sqrt(DSQR(xloc[0])+DSQR(xloc[1])+DSQR(xloc[2]));
   for (j=0;j<3;j++) xloc[j] /=norm;
   norm = sqrt(DSQR(yloc[0])+DSQR(yloc[1])+DSQR(yloc[2]));
   for (j=0;j<3;j++) yloc[j] /=norm;
   norm = sqrt(DSQR(zloc[0])+DSQR(zloc[1])+DSQR(zloc[2]));
   for (j=0;j<3;j++) zloc[j] /=norm;

   /*----------------------------------------------- plausibility checks */
   /* vectors have to be orthogonal: check scalar products */
   if (FABS(xloc[0]*yloc[0]+xloc[1]*yloc[1]+xloc[2]*yloc[2])>EPS8)
      dserror("locsys base vectors are not orthogonal!\n");
   if (FABS(xloc[0]*zloc[0]+xloc[1]*zloc[1]+xloc[2]*zloc[2])>EPS8)
      dserror("locsys base vectors are not orthogonal!\n");
   if (FABS(zloc[0]*yloc[0]+zloc[1]*yloc[1]+zloc[2]*yloc[2])>EPS8)
      dserror("locsys base vectors are not orthogonal!\n");

   /* we need a right hand system: check vector product */
   if (FABS(zloc[0]-(xloc[1]*yloc[2] - xloc[2]*yloc[1]))
      +FABS(zloc[1]-(xloc[2]*yloc[0] - xloc[0]*yloc[2]))
      +FABS(zloc[2]-(xloc[0]*yloc[1] - xloc[1]*yloc[0]))>EPS8)
       dserror("chosen local cosys no right hand system!\n");

   /*--------------------------------------- store the direction cosini */
   actlocsys->lXx=xloc[0];
   actlocsys->lXy=yloc[0];
   actlocsys->lXz=zloc[0];
   actlocsys->lYx=xloc[1];
   actlocsys->lYy=yloc[1];
   actlocsys->lYz=zloc[1];
   actlocsys->lZx=xloc[2];
   actlocsys->lZy=yloc[2];
   actlocsys->lZz=zloc[2];
}

#if 0
if (fdyn->freesurf==6)
{
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[disnum].node[i]);
      actgnode=actnode->gnode;
      /*---------------------------------------- check for free surface */
      if (actgnode->freesurf==NULL) continue;
      /*---------------------------------------- check for local co-sys */
      ilocsys=actnode->locsysId-1;
      if (ilocsys>=0)
         dserror("locsys at generalised free surface not possible yet!\n");
      /*------------ node has individual local co-sys based on the mass
         consistent normal
      */
      numlocsys++;
      locsys=(LOCSYS*)CCAREALLOC(locsys,numlocsys*sizeof(LOCSYS));
      actlocsys=&(locsys[numlocsys-1]);
      actlocsys->Id=numlocsys;
      actnode->locsysId=numlocsys;
      actlocsys->locsystyp = locsys_fsnd;
      for (l=0;l<actnode->numele;l++)
      {
         actele=actnode->element[l];
         actele->locsys=locsys_yes;
      } /* end loop over nodes */
      /*----------------------------------------- allocate base vectors */
      xloc=amdef("xloc",&actlocsys->xloc,3,1,"DV");
      yloc=amdef("yloc",&actlocsys->yloc,3,1,"DV");
      zloc=amdef("zloc",&actlocsys->zloc,3,1,"DV");
      if (fdyn->numdf==3)
      {
         /* 2D-case is straight forward: normal = local y-axis */
         yloc[0]=actnode->actn[0];
         yloc[1]=actnode->actn[1];
         yloc[2]=actnode->actn[2];
         dsassert(FABS(yloc[2])<EPS15,"3D base vector y in 2D case!\n");
         /* define xloc orthogonal in xy-plane to yloc */
         xloc[0] = -yloc[1];
         xloc[1] =  yloc[0];
         xloc[2] = ZERO;
         /* vector product to check zloc */
         zloc[0] = xloc[1]*yloc[2] - xloc[2]*yloc[1];
         zloc[1] = xloc[2]*yloc[0] - xloc[0]*yloc[2];
         zloc[2] = xloc[0]*yloc[1] - xloc[1]*yloc[0];
         if (zloc[2]<ZERO) /* change sign of yloc */
         {
            xloc[0] *=-ONE;
            xloc[1] *=-ONE;
            zloc[2] *=-ONE;
         }
      }
      else if (fdyn->numdf==4) /* 3D how to choose the tangents??? */
         dserror("fluid locsys for 3D case not implemented yet!\n");
      else
         dserror("fluid locsys for 4D and higer dim. cases not implemented yet!\n");
      /*------------------------------------- normalise the base vectors */
      norm = sqrt(DSQR(xloc[0])+DSQR(xloc[1])+DSQR(xloc[2]));
      for (j=0;j<3;j++) xloc[j] /=norm;
      norm = sqrt(DSQR(yloc[0])+DSQR(yloc[1])+DSQR(yloc[2]));
      for (j=0;j<3;j++) yloc[j] /=norm;
      norm = sqrt(DSQR(zloc[0])+DSQR(zloc[1])+DSQR(zloc[2]));
      for (j=0;j<3;j++) zloc[j] /=norm;

      /*-------------------------------------------- plausibility checks */
      /* vectors have to be orthogonal: check scalar products */
      if (FABS(xloc[0]*yloc[0]+xloc[1]*yloc[1]+xloc[2]*yloc[2])>EPS8)
         dserror("locsys base vectors are not orthogonal!\n");
      if (FABS(xloc[0]*zloc[0]+xloc[1]*zloc[1]+xloc[2]*zloc[2])>EPS8)
         dserror("locsys base vectors are not orthogonal!\n");
      if (FABS(zloc[0]*yloc[0]+zloc[1]*yloc[1]+zloc[2]*yloc[2])>EPS8)
         dserror("locsys base vectors are not orthogonal!\n");

      /* we need a right hand system: check vector product */
      if (FABS(zloc[0]-(xloc[1]*yloc[2] - xloc[2]*yloc[1]))
         +FABS(zloc[1]-(xloc[2]*yloc[0] - xloc[0]*yloc[2]))
         +FABS(zloc[2]-(xloc[0]*yloc[1] - xloc[1]*yloc[0]))>EPS8)
         dserror("chosen local cosys no right hand system!\n");

      /*------------------------------------ store the direction cosini */
      actlocsys->lXx=xloc[0];
      actlocsys->lXy=yloc[0];
      actlocsys->lXz=zloc[0];
      actlocsys->lYx=xloc[1];
      actlocsys->lYy=yloc[1];
      actlocsys->lYz=zloc[1];
      actlocsys->lZx=xloc[2];
      actlocsys->lZy=yloc[2];
      actlocsys->lZz=zloc[2];
   }
}
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_locsys */

#endif
/*! @} (documentation module close)*/
