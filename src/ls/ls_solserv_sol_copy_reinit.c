#include "../headers/standardtypes.h"
#include "../solver/solver.h"



/*----------------------------------------------------------------------*
 |                                                            m.gee 4/03|
 | copy values from array arrayfrom in place from      to               |
 |                  array arrayto   in place to                         |
 |                                                                      |
 | FIELD *actfield  (i) active field                                    |
 | INT    disnum    (i) indize of the discretization in actfield to be used|
 | INT    arrayfrom (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 | INT    arrayto   (i) indize of the array, 0 = sol                    |
 |                                           1 = sol_increment          |
 |                                           2 = sol_residual           |
 | INT    from    (i) row in ARRAY sol to be set to zero                |
 *----------------------------------------------------------------------*/
void solserv_sol_copy_reinit(
  FIELD *actfield, INT disnum, INT arrayfrom, INT arrayto,
  INT from, INT to
  )
{
INT               i,j;
INT               diff,max;
ARRAY            *arrayf,*arrayt;
NODE             *actnode;
DISCRET          *actdis;
#ifdef DEBUG
dstrc_enter("solserv_sol_copy_reinit");
#endif
/*----------------------------------------------------------------------*/

printf("\n**WARNING** COPYING SOLUTION IN REINITIALIZATION!\n");

actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
/*********************BE CAREFUL*****************************************/
/*********************BE CAREFUL*****************************************/
/*********************BE CAREFUL*****************************************/
   if (actnode->gnode->is_node_active==1) continue;
/*********************BE CAREFUL*****************************************/
/*********************BE CAREFUL*****************************************/
/*********************BE CAREFUL*****************************************/
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayfrom)
   {
   case 0:
      arrayf = &(actnode->sol);
   break;
   case 1:
      arrayf = &(actnode->sol_increment);
   break;
   case 2:
      arrayf = &(actnode->sol_residual);
   break;
   case 3:
      arrayf = &(actnode->sol_mf);
   break;
   default:
      arrayf = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayto)
   {
   case 0:
      arrayt = &(actnode->sol);
   break;
   case 1:
      arrayt = &(actnode->sol_increment);
   break;
   case 2:
      arrayt = &(actnode->sol_residual);
   break;
   case 3:
      arrayt = &(actnode->sol_mf);
   break;
   default:
      arrayt = NULL;
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of arrayf */
   if (from >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");
   /* check the size of arrayt */
   if (to >= arrayt->fdim)
   {
      diff = to - arrayt->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(arrayt,arrayt->fdim+max+1,arrayt->sdim,"DA");
   }
   for (j=0; j<arrayt->sdim; j++)
      arrayt->a.da[to][j] = arrayf->a.da[from][j];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of solserv_sol_copy_reinit */
