/*!---------------------------------------------------------------------
\file controlling lift&drag calculation
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

---------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;  
/*!--------------------------------------------------------------------- 
\brief controlling lift and drag calculation

<pre>                                                         genk 12/03

</pre>


\return void                                            

------------------------------------------------------------------------*/
void fluid_liftdrag(INT            init,
                    CALC_ACTION   *action,
		    CONTAINER     *container,
		    FIELD         *actfield,
		    SOLVAR        *actsolv,
		    PARTITION     *actpart,
		    INTRA         *actintra,
		    FLUID_DYNAMIC *fdyn)
{
INT           i,j,k;
INT           actcurve,numline;
FIELD        *structfield;
DLINE        *actdline,*aledline;
NODE         *actnode;
DOUBLE	      liftdrag[(FLUID_NUM_LD+1)*12];	/* array with lift & drag coeff.*/
DOUBLE        acttimefac;
DOUBLE        timefac[5];
DOUBLE        initval;

#ifdef PARALLEL
DOUBLE        recv[(FLUID_NUM_LD+1)*12];          /*  receive buffer              */
#endif

#ifdef DEBUG 
dstrc_enter("fluid_liftdrag");
#endif


/*****************************************************************
 *
 * W A R N I N G ! ! !
 *
 * lift and drag for multifield (pseudo fsi and real fsi)
 * only for 2 D  !!!
 *
 *****************************************************************/


switch (init)
{
case 0: /* init multifield */
structfield = &(field[genprob.numsf]);
for (i=0; i<design->ndline; i++) /* loop over dlines */
{
   actdline = &(design->dline[i]);
   if (actdline->liftdrag == NULL) continue;
   for (k=0;k<structfield->dis[0].numnp;k++) /* loop over nodes */
   {
      actnode = &(structfield->dis[0].node[k]);
      if (actnode->x[0] - actdline->liftdrag->ld_center[0] < EPS8)
      if (actnode->x[1] - actdline->liftdrag->ld_center[1] < EPS8)
      if (actnode->x[2] - actdline->liftdrag->ld_center[2] < EPS8)
      {
         actdline->liftdrag->alenode = k;
         break;
      }
   } /* end loop over nodes */
}/* end loop over dlines */
break;

case 1: /* evaluate for fluid problem */
   /* do nothing at the moment */
break;

case 2: /* evaluate for pseudo fsi problem */
   for (i=0; i<design->ndline; i++) /* loop over dlines */
   {
      actdline = &(design->dline[i]);
      if (actdline->liftdrag == NULL) continue;
      numline = actdline->liftdrag->aledline;
      aledline = &(design->dline[numline]);
      for (j=0;j<2;j++)
      {
         actcurve = aledline->dirich->curve.a.iv[j]-1; 
         if (actcurve<0) acttimefac = ONE;
         else acttimefac = timefac[actcurve];
         initval  = aledline->dirich->dirich_val.a.dv[j];               
         actdline->liftdrag->ld_center[j] += initval*acttimefac;  
      }
   }
break;

case 3: /* evaluate for real fsi problem */
   structfield = &(field[genprob.numsf]);
   /*------------------------------- update center point coordinates ---*/
   for (i=0; i<design->ndline; i++) /* loop over dlines */
   {
      actdline = &(design->dline[i]);
      if (actdline->liftdrag == NULL) continue;
      actnode = &(structfield->dis[0].node[actdline->liftdrag->alenode]);
      actdline->liftdrag->ld_center[0] = actnode->x[0]+actnode->sol_mf.a.da[0][0];
      actdline->liftdrag->ld_center[1] = actnode->x[1]+actnode->sol_mf.a.da[0][1];
      actdline->liftdrag->ld_center[2] = actnode->x[1]+actnode->sol_mf.a.da[0][2];
   }	/* end loop over dlines */
break;
}

if (init>0)
{
  /*---------------- initialise lift- and drag- coefficient to zero ---*/
  for (i=0; i<(FLUID_NUM_LD+1)*12; i++)
  {
   liftdrag[i] = ZERO;
  }
   
   /*------------------------- get fluid stresses and integrate them ---*/
   container->liftdrag = liftdrag;
   container->nii= 0;
   container->nim= 0;
   container->nif= 0;
   calelm(actfield,actsolv,actpart,actintra,0,-1,
          container,action);

   /*----------- distribute lift- and drag- coefficient to all procs ---*/
#ifdef PARALLEL
   MPI_Reduce(liftdrag,recv,12*(FLUID_NUM_LD+1),MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
   for (i=0; i<(FLUID_NUM_LD+1)*12; i++)
       liftdrag[i] = recv[i];
#endif

   /* viscous */
   /* calculate the sums of all liftdrag conditions */ 
   for (j=0; j<6; j++) 
     for (i=0; i<FLUID_NUM_LD; i++)
       liftdrag[FLUID_NUM_LD*6+j] += liftdrag[i*6+j];
   
   /* only pressure */
   /* calculate the sums of all liftdrag conditions */ 
   for (j=0; j<6; j++) 
     for (i=0; i<FLUID_NUM_LD; i++)
       liftdrag[(2*FLUID_NUM_LD+1)*6+j] += liftdrag[(FLUID_NUM_LD+1)*6+i*6+j];
   
   /*-------------------------------------------------------- output ---*/
   if(par.myrank == 0)
   {
     /* viscous */
     printf("F_x = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+0]);
     printf("\n");

     printf("F_y = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+1]);
     printf("\n");

     printf("F_z = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+2]);
     printf("\n");

     printf("M_x = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+3]);
     printf("\n");

     printf("M_y = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+4]);
     printf("\n");

     printf("M_z = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[i*6+5]);
     printf("\n\n");

     /* only pressure */
     printf("F_x = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+0]);
     printf("\n");

     printf("F_y = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+1]);
     printf("\n");

     printf("F_z = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+2]);
     printf("\n");

     printf("M_x = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+3]);
     printf("\n");

     printf("M_y = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+4]);
     printf("\n");

     printf("M_z = ");
     for (i=0; i<FLUID_NUM_LD+1; i++)
       printf(" %12.4E ",  liftdrag[(FLUID_NUM_LD +1+i)*6+5]);
     printf("\n\n");
   }
   if (par.myrank == 0)
     plot_liftdrag(fdyn->time, liftdrag);
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_liftdrag*/
#endif
