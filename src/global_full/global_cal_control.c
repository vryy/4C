#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
int i;
FIELD *actfield;

#ifdef DEBUG 
dstrc_enter("ntacal");
#endif
/*----------------------------------------------------------------------*/
/*------------------------do initial partitioning of nodes and elements */
part_fields();

/*------------------------------------------------ assign dofs to nodes */
for(i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]); 
   if (actfield->ndis==1) assign_dof(actfield);
   if (actfield->ndis>1) assign_dof_ndis(actfield);
}
/*--------make the procs know their own nodes and elements a bit better */
part_assignfield();
 
/*-------------------calculate system matrices parallel storage formats */
mask_global_matrices();
/*------------------------------------------------ write general output */
out_general();
/*--------------------------------------------------- write mesh to gid */
if (par.myrank==0) 
{
   if (ioflags.struct_disp_gid||ioflags.struct_stress_gid||ioflags.fluid_sol_gid) 
   {
      out_gid_sol_init();
      out_gid_msh();
   }
}
/*------------------ call control programs of static or dynamic control */
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
