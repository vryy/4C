/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
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
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/02    |
 | write the data needed to restart this element                        |
 *----------------------------------------------------------------------*/
void s8_write_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
FILE *out;
#ifdef DEBUG 
dstrc_enter("s8_write_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the shell element has to write are the eas parameters
   at the moment
*/
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
if (actele->e.s8->nhyb)
{
   handles[0]+=4;
   if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
   /*------------------------------------------------------- write alfa */
   pss_write_array(&(actele->e.s8->alfa),&(handles[1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*--------------------------------------------------- write Dtildinv */
   pss_write_array(&(actele->e.s8->Dtildinv),&(handles[2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*--------------------------------------------------------- write Lt */
   pss_write_array(&(actele->e.s8->Lt),&(handles[3]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*----------------------------------------------------- write Rtilde */
   pss_write_array(&(actele->e.s8->Rtilde),&(handles[4]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_write_restart */
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/02    |
 | read the data needed to restart this element                         |
 *----------------------------------------------------------------------*/
void s8_read_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
INT dims[3];
FILE *in;
#ifdef DEBUG 
dstrc_enter("s8_read_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the shell element has to read are the eas parameters
   at the moment
*/
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
/*----------------------------------------------------------------------*/
if (actele->e.s8->nhyb)
{
   /*------------------------------------------------- check dimensions alfa*/
   pss_getdims_name_handle(actele->e.s8->alfa.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->alfa.fdim != dims[0] ||
       actele->e.s8->alfa.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*-------------------------------------------------------- read alfa */
   pss_read_array_name_handle(actele->e.s8->alfa.name,&(actele->e.s8->alfa),&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Dtildinv*/
   pss_getdims_name_handle(actele->e.s8->Dtildinv.name,&dims[0],&dims[1],&dims[2],&handles[2],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Dtildinv.fdim != dims[0] ||
       actele->e.s8->Dtildinv.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------- read Dtildinv */
   pss_read_array_name_handle(actele->e.s8->Dtildinv.name,&(actele->e.s8->Dtildinv),&handles[2],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Lt*/
   pss_getdims_name_handle(actele->e.s8->Lt.name,&dims[0],&dims[1],&dims[2],&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Lt.fdim != dims[0] ||
       actele->e.s8->Lt.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------------- read Lt */
   pss_read_array_name_handle(actele->e.s8->Lt.name,&(actele->e.s8->Lt),&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Rtilde*/
   pss_getdims_name_handle(actele->e.s8->Rtilde.name,&dims[0],&dims[1],&dims[2],&handles[4],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Rtilde.fdim != dims[0] ||
       actele->e.s8->Rtilde.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*------------------------------------------------------ read Rtilde */
   pss_read_array_name_handle(actele->e.s8->Rtilde.name,&(actele->e.s8->Rtilde),&handles[4],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_read_restart */
#endif
