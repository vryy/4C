/*!----------------------------------------------------------------------
\file
\brief write and read element restart data of FLUID3

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID3
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;

/*!---------------------------------------------------------------------
\brief write element restart data of fluid3 element

<pre>                                                         genk 09/03

   at the moment only the curvature at the free surface has to be
   written to the binary file

</pre>
\param  *actele	   ELEMENT	   (i)    actual element
\param   nhandle   INT             (i)    number of handles
\param  *handles   long int        (i)    handles

\return void

------------------------------------------------------------------------*/
void f3_write_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
INT counter=1;
FLUID_DYNAMIC   *fdyn;
FILE *out;
#ifdef DEBUG
dstrc_enter("f3_write_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   fluid3 element has to write the curvature data at the free surface
   and the stab parameter of the last time step.
   
*/
/*----------------------------------------------------------------------*/
fdyn= alldyn[genprob.numff].fdyn;
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
handles[0]+=1;
if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
/*--------------------------------------------------------- write kappa */
#if 0 /* no KAPPA at the moment for fluid3 element! */
if (fdyn->surftens>0 && actele->e.f3->fs_on>0)
{
   pss_write_array(&(actele->e.f3->kappa_ND),&(handles[counter]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   counter++;
}
#endif
/*----------------------------------------------------- write stab-par */
pss_write_array(&(actele->e.f3->tau_old),&(handles[counter]),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
counter++;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_write_restart */

/*!---------------------------------------------------------------------
\brief read element restart data of fluid3 element

<pre>                                                         genk 09/03

   at the moment only the curvature at the free surface has to be
   read from the binary file

</pre>
\param  *actele	   ELEMENT	   (i)    actual element
\param   nhandle   INT             (i)    number of handles
\param  *handles   long int        (i)    handles

\return void

------------------------------------------------------------------------*/
void f3_read_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
INT dims[3];
INT          counter=1;
FLUID_DYNAMIC   *fdyn;
FILE *in;
#ifdef DEBUG
dstrc_enter("f3_read_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   fluid3 element has to read the curvature data at the free surface
   and the stab parameter of the last time step.
   
*/
/*----------------------------------------------------------------------*/
fdyn= alldyn[genprob.numff].fdyn;
in  = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
#if 0 /* no KAPPA at the moment for fluid3 element! */
/*----------------------------------------------------------------------*/
if (fdyn->surftens!=0 && actele->e.f3->fs_on>0)
{
   /*-------------------------------------------- check dimensions kappa*/
   pss_getdims_name_handle(actele->e.f3->kappa_ND.name,&dims[0],&dims[1],
                           &dims[2],&handles[counter],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.f3->kappa_ND.fdim != dims[0] ||
       actele->e.f3->kappa_ND.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*------------------------------------------------------- read kappa */
   pss_read_array_name_handle(actele->e.f3->kappa_ND.name,
      &(actele->e.f3->kappa_ND),&handles[counter],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   counter++;
}
#endif
/*-------------------------------------------- check dimensions tau_old */
pss_getdims_name_handle(actele->e.f3->tau_old.name,&dims[0],&dims[1],
                        &dims[2],&handles[counter],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
if (actele->e.f3->tau_old.fdim != dims[0] ||
    actele->e.f3->tau_old.sdim != dims[1])
    dserror("Mismatch in reading element restart data");
/*-------------------------------------------------------- read tau_old */
pss_read_array_name_handle(actele->e.f3->tau_old.name,
   &(actele->e.f3->tau_old),&handles[counter],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
counter++;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_read_restart */

#endif
