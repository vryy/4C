/*!----------------------------------------------------------------------
\file
\brief ls_result_incre.c

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
\brief map solution obtained back to nodes and compute convergence ratio

<pre>                                                            irhan 05/04
map solution obtained back to nodes and compute convergence ratio
</pre>

*----------------------------------------------------------------------*/
void ls_result_incre(
  FIELD             *actfield,
  INTRA             *actintra,
  DIST_VECTOR       *sol,
  INT                place,
  SPARSE_ARRAY      *sysarray,
  SPARSE_TYP        *sysarray_typ,
  DOUBLE            *lrat,
  LS_DYNAMIC       *lsdyn
  )
{
  INT         i,j;
  INT         max;
  INT         diff;
  INT         dof;
  INT         numeq_total;
  NODE       *actnode;
  ARRAY       result_a;
  DOUBLE     *result; /* redundant result vector */
/************************************************************************/
  DOUBLE      lnorm  =ZERO;
  DOUBLE      lnorm00=ZERO;
  DOUBLE      lnorm01=ZERO;
  DOUBLE      lnorm02=ZERO;
/************************************************************************/

#ifdef DEBUG
  dstrc_enter("ls_result_incre");
#endif
/*----------------------------------------------------------------------*/

  numeq_total = sol->numeq_total;
  /* allocate space to allreduce the DIST_VECTOR */
  result = amdef("result",&result_a,numeq_total,1,"DV");
  amzero(&result_a);
  /* copy distributed result to redundant result vector */
  solserv_reddistvec(
    sol,
    sysarray,
    sysarray_typ,
    result,
    sol->numeq_total,
    actintra
    );
  /* loop nodes and put the result back to the node structure */
  for (i=0; i<actfield->dis[0].numnp; i++)
  {
    actnode = &(actfield->dis[0].node[i]);
    /* enlarge sol_increment, if necessary */
    if (place >= actnode->sol_increment.fdim)
    {
      diff = place - actnode->sol_increment.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
              actnode->sol_increment.sdim,"DA");
    }
    for (j=0; j<actnode->numdf; j++)
    {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
/*****************************BE CAREFUL*********************************/
      if (actnode->sol_increment.a.da[place][j]==0.0) continue;
/*****************************BE CAREFUL*********************************/
      /* compute L-inf     norm */
      lnorm00  = DMAX(lnorm, FABS(result[dof])/actnode->sol_increment.a.da[place][j]);
      /* compute L-1       norm */
      lnorm01 += FABS(result[dof]/actnode->sol_increment.a.da[place][j]);
      /* compute L-2       norm */
      lnorm02 += pow(result[dof]/actnode->sol_increment.a.da[place][j],2);
      /* put result to the node */
      actnode->sol_increment.a.da[place][j] += result[dof];
    }
  }
  switch (lsdyn->itchk) /* iteration convergence check */
  {
      case 0:
	/* no convergence check */
	goto end;
	break;
      case 1:
	/* convergence check by L-infinity norm */
        lnorm =  lnorm00;
	break;
      case 2:
	/* convergence check by L-1        norm */
        lnorm =  lnorm01;
	break;
      case 3:
	/* convergence check by L-2        norm */
        lnorm = sqrt( lnorm02);
	break;
      default:
        dserror("unknown norm for convergence check!\n");
  }
  /* set final convergence ratios */
  *lrat = lnorm;

 end:
  amdel(&result_a);

/*-------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of  ls_result_incre */
/*! @} (documentation module close)*/
#endif
