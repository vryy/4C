/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
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
INT i;
FIELD *actfield;

#ifdef DEBUG 
dstrc_enter("ntacal");
#endif
/*----------------------------------------------------------------------*/
/*------------------------do initial partitioning of nodes and elements */
part_fields();

#ifdef D_FLUID
/*---------------------------------- set dofs for implicit free surface */
if (genprob.numff>=0) fluid_freesurf_setdofs();
/*----------------------------- modify coordinates for special problems */
if (genprob.numff>=0) fluid_modcoor();
#endif

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
   if (ioflags.struct_disp_gid||ioflags.struct_stress_gid||ioflags.fluid_sol_gid
       ||ioflags.ale_disp_gid) 
   {
      out_gid_sol_init();
      if (field[genprob.numsf].fieldtyp != structure)
         out_gid_msh();
   }
}
/*------------------------ program to control execution of optimization */
#ifdef D_OPTIM                   /* include optimization code to ccarat */
if(genprob.probtyp==prb_opt)
{
  caloptmain();
  goto end;
}   
#endif
/*------------------ call control programs of static or dynamic control */
if (genprob.timetyp==time_static)
{
     calsta();
}
if (genprob.timetyp==time_dynamic)
{
     caldyn(); 
}
/*------------------------------------------------------- check results */
#ifdef RESULTTEST
global_result_test();
#endif

/*----------------------------------------------------------------------*/
end:;
/*--------------------------------------------------- write warnings ---*/
dswarning(2,0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ntacal */
