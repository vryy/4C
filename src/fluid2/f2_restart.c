/*!----------------------------------------------------------------------
\file
\brief write and read element restart data of FLUID2

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
\brief write element restart data of fluid2 element

<pre>                                                         genk 09/03

   at the moment only the curvature at the free surface has to be 
   written to the binary file
		     
</pre>
\param  *actele	   ELEMENT	   (i)    actual element
\param   nhandle   INT             (i)    number of handles
\param  *handles   long int        (i)    handles

\return void                                                                       

------------------------------------------------------------------------*/
void f2_write_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
FILE *out;
#ifdef DEBUG 
dstrc_enter("f2_write_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the fluid2 element has to write are the curvature
   kappa at the free surface
*/
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
if (actele->e.f2->fs_on>0)
{
   handles[0]+=1;
   if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
   /*------------------------------------------------------ write kappa */
   pss_write_array(&(actele->e.f2->kappa_ND),&(handles[1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_write_restart */ 

/*!---------------------------------------------------------------------
\brief read element restart data of fluid2 element

<pre>                                                         genk 09/03

   at the moment only the curvature at the free surface has to be 
   read from the binary file
		     
</pre>
\param  *actele	   ELEMENT	   (i)    actual element
\param   nhandle   INT             (i)    number of handles
\param  *handles   long int        (i)    handles

\return void                                                                       

------------------------------------------------------------------------*/
void f2_read_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
INT dims[3];
FILE *in;
#ifdef DEBUG 
dstrc_enter("f2_read_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the fluid2 element has to read are the curvature
   kappa at the free surface
*/
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
/*----------------------------------------------------------------------*/
if (actele->e.f2->fs_on>0)
{
   /*-------------------------------------------- check dimensions kappa*/
   pss_getdims_name_handle(actele->e.f2->kappa_ND.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.f2->kappa_ND.fdim != dims[0] ||
       actele->e.f2->kappa_ND.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*-------------------------------------------------------- read alfa */
   pss_read_array_name_handle(actele->e.f2->kappa_ND.name,&(actele->e.f2->kappa_ND),&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_read_restart */

#endif
