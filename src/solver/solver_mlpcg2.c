/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;




/*!---------------------------------------------------------------------
\brief init multilevel preconditioner                                              

<pre>                                                        m.gee 6/02 

</pre>
\param actsolv    SOLVAR*      (i)   general structure of solver informations                   
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_init(MLPCGVARS *mlpcgvars)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_init");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_init */































/*! @} (documentation module close)*/
