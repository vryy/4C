/*!----------------------------------------------------------------------
\file
\brief contains the routine 'caloptmain',
       controlling execution of optimization

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
  int  i;            /* a counter */
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("caloptmain");
  #endif
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
