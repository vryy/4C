/*!----------------------------------------------------------------------
\file
\brief ls_convcheck.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;



/*!----------------------------------------------------------------------
\brief subroutine for convergence check of level set problem

<pre>                                                            irhan 05/04
subroutine for convergence check of level set problem
</pre>

*----------------------------------------------------------------------*/
INT ls_convcheck(
  LS_DYNAMIC	*lsdyn,
  DOUBLE       lrat,
  INT          itnum,
  DOUBLE       te,
  DOUBLE       ts
  )
{
  INT     converged=0;  /* flag for convergence check */

#ifdef DEBUG
  dstrc_enter("ls_convcheck");
#endif
/*----------------------------------------------------------------------*/

  if (lsdyn->itchk!=0)
  {
    if (par.myrank==0) /* output to the screen */
    {
      switch(lsdyn->itchk)
      {
          case 1: /* infinity norm */
            printf("|  %3d/%3d   | %10.3E[L_in]   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                   itnum,lsdyn->itemax,lsdyn->ittol,lrat,te,ts);
            break;
          case 2: /* L_1 norm */
            printf("|  %3d/%3d   | %10.3E[L_in]   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                   itnum,lsdyn->itemax,lsdyn->ittol,lrat,te,ts);
            break;
          case 3: /* L_2 norm */
            printf("|  %3d/%3d   | %10.3E[L_in]   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n",
                   itnum,lsdyn->itemax,lsdyn->ittol,lrat,te,ts);
            break;
          default:
            dserror("Norm for nonlin. convergence check unknown!!\n");
      }
    }
    /* convergence check */
    if (lrat<lsdyn->ittol)
      converged=2;
    if (itnum==lsdyn->itemax)
      converged++;
    if (converged==1 && par.myrank==0)
    {
      printf("---------------------------------------------------------------- \n");
      printf("|          >>>>>> not converged in itemax steps!               | \n");
    }
    if (converged>0 && par.myrank==0)
    {
      printf("---------------------------------------------------------------- \n");
      printf("\n");
    }
  }
  else
  {
    if (par.myrank==0)
    {
      printf("      iteration step: %3d / 3d \n",
             itnum, lsdyn->itemax);
      if(itnum==lsdyn->itemax)
      {
       printf("----------------------------------------------------------------- \n");
       printf("\n");
      }
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return (converged);
} /* end of ls_convcheck*/
/*! @} (documentation module close)*/
#endif
