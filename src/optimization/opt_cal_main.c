/*!----------------------------------------------------------------------
\file
\brief contains the routine 'caloptmain',
       controlling execution of optimization

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/optimization.h"
#include "opt_prototypes.h"
/*!
\addtogroup OPTIMIZATION
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
\brief the optimization main structure
<pre>                                                            al 06/01
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 | control execution of optimization                        al 05/01    |
 *----------------------------------------------------------------------*/
void caloptmain()
{
/*----------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("caloptmain");
  #endif
/*----------------------------------------------------------------------*/
  if(opt->objective==oj_frequency&&par.nprocs>1)
  {
    printf(" no optimization in parallel implemented");
    printf(" for eigen frequencies\n");
    goto endoptmain;
  }
/*-------------------------------------- initialize graphical output ---*/
  opt_g_out(gr_init);
/*----------------------- initialize execution stage of optimization ---*/
  opcini();
/*------------------------- sequential quadratic programming (nlpql) ---*/
  if (opt->strategy==os_nlp)
  {
  }
/*-------------------------------------------- fully stressed design ---*/
  if (opt->strategy==os_fsd)
  {
    optfsd();
  }
/*----------------------------------------------------------------------*/
  endoptmain:;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_exit();
  #endif
return;
} /* end of caloptmain */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
