#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{

#ifdef DEBUG 
dstrc_enter("ntacal");
#endif
/*----------------------------------------------------------------------*/
/*------------------------do initial partitioning of nodes and elements */
part_fields();
/*------------------------------------------------ assign dofs to nodes */
assign_dof();
/*--------make the procs know their own nodes and elements a bit better */
part_assignfield();
/*-------------------calculate system matrices parallel storage formats */
mask_global_matrices();
/*------------------------------------------------ write general output */
out_general();
/*--------------------------------------------------- write mesh to gid */
if (par.myrank==0) 
{
   out_gid_sol_init();
   out_gid_msh();
}
/*----------------- call controll programs of static or dynamic control */
if (genprob.timetyp==time_static)
{
     calsta();
}
if (genprob.timetyp==time_dynamic)
{
     caldyn(); 
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ntacal */
