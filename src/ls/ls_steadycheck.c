#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "ls_prototypes.h"



extern struct _PAR   par;





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
INT ls_steadycheck(
  LS_DYNAMIC	*lsdyn,
  FIELD		*actfield,
  INT		 numeq_total
  )
{
  INT        steady=0;
  DOUBLE     lrat;

#ifdef DEBUG
  dstrc_enter("ls_steadycheck");
#endif
/*----------------------------------------------------------------------*/

  /* determine the conv. ratios */
  ls_norm(lsdyn,actfield,numeq_total,&lrat);
  /* output to the screen */
  if (par.myrank==0)
  {
    switch (lsdyn->stchk)
    {
        case 1:
          printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_in] \n",
                 lsdyn->sttol);
          break;
        case 2:
          printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_1 ] \n",
                 lsdyn->sttol);
          break;
        case 3:
          printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_2 ] \n",
                 lsdyn->sttol);
          break;
        default:
          dserror("Norm for steady state check unknwon!\n");
    }
    printf("         levelset: %10.3E  \n",lrat);
  }
  /* check if the ratios are smaller than the given tolerance and set flag */
  if (lrat < lsdyn->sttol)
  {
    steady=1;
    if (par.myrank==0)
    {
      printf("\n");
      printf("    >>>>>> STEADY STATE REACHED <<<<<< \n");
      printf("\n");
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return (steady);
} /* end of ls_steadycheck*/



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls_norm(
  LS_DYNAMIC      *lsdyn,
  FIELD           *actfield,
  INT              numeq_total,
  DOUBLE          *lrat
  )
{
  INT       i,j;
  INT       dof;
  NODE     *actnode;
/***********************************************************************/
  DOUBLE   dlnorm  =ZERO;
  DOUBLE    lnorm  =ZERO;
  DOUBLE   dlnorm00=ZERO;
  DOUBLE    lnorm00=ZERO;
  DOUBLE   dlnorm01=ZERO; /* variables for norm calculation */
  DOUBLE    lnorm01=ZERO;
  DOUBLE   dlnorm02=ZERO;
  DOUBLE    lnorm02=ZERO;
/************************************************************************/

#ifdef DEBUG
  dstrc_enter("ls_norm");
#endif

  /* loop over all nodes */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode=&(actfield->dis[0].node[i]);
    for (j=0; j<actnode->numdf; j++)
    {
      dof = actnode->dof[j];
      if (dof >= numeq_total) continue;
      /* compute L-infinity norm */
      dlnorm00  = DMAX(dlnorm, FABS(actnode->sol_increment.a.da[1][j]- \
                                    actnode->sol_increment.a.da[0][j]));
      lnorm00  = DMAX( lnorm, FABS(actnode->sol_increment.a.da[0][j]));

      /* compute L-1        norm */
      dlnorm01 += FABS(actnode->sol_increment.a.da[1][j]- \
                       actnode->sol_increment.a.da[0][j]);
      lnorm01 += FABS(actnode->sol_increment.a.da[0][j]);

      /* compute L-2        norm */
      dlnorm02 += pow(actnode->sol_increment.a.da[1][j]- \
                      actnode->sol_increment.a.da[0][j],2);
      lnorm02 += pow(actnode->sol_increment.a.da[0][j],2);
    } /* end of loop over dofs */
  } /* end of loop over nodes */


  switch (lsdyn->stchk) /* steady state check */
  {
      case 0:
	/* no convergence check */
	goto end;
	break;
      case 1:
	/* convergence check by L-infinity norm */
	dlnorm = dlnorm00;
        lnorm =  lnorm00;
	break;
      case 2:
	/* convergence check by L-1        norm */
	dlnorm = dlnorm01;
        lnorm =  lnorm01;
	break;
      case 3:
	/* convergence check by L-2        norm */
	dlnorm = sqrt(dlnorm02);
        lnorm = sqrt( lnorm02);
	break;
      default:
        dserror("unknown norm for stedy state check!\n");
  }

  /* check for "ZERO-field" */
  if (lnorm<EPS5)
  {
    lnorm = ONE;
#ifdef DEBUT
    printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
  }
  /* set final convergence ratios */
  *lrat = dlnorm/lnorm;

 end:
/*-------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of  ls_norm */
#endif
