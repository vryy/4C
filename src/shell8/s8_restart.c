#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/02    |
 | write the data needed to restart this element                        |
 *----------------------------------------------------------------------*/
void s8_write_restart(ELEMENT *actele, int nhandle, int *handles)
{
int ierr;

#ifdef DEBUG 
dstrc_enter("s8_write_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the shell element has to write are the eas parameters
   at the moment
*/
/*----------------------------------------------------------------------*/
handles[0]=0;
if (actele->e.s8->nhyb)
{
   handles[0]+=4;
   if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
   /*------------------------------------------------------- write alfa */
   pss_write_array(&(actele->e.s8->alfa),&(handles[1]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*--------------------------------------------------- write Dtildinv */
   pss_write_array(&(actele->e.s8->Dtildinv),&(handles[2]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*--------------------------------------------------------- write Lt */
   pss_write_array(&(actele->e.s8->Lt),&(handles[3]),&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*----------------------------------------------------- write Rtilde */
   pss_write_array(&(actele->e.s8->Rtilde),&(handles[4]),&ierr);
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
void s8_read_restart(ELEMENT *actele, int nhandle, int *handles)
{
int ierr;
int dims[3];
#ifdef DEBUG 
dstrc_enter("s8_read_restart");
#endif
/*----------------------------------------------------------------------*/
/*
   the only thing the shell element has to read are the eas parameters
   at the moment
*/
/*----------------------------------------------------------------------*/
if (actele->e.s8->nhyb)
{
   /*------------------------------------------------- check dimensions alfa*/
   pss_getdims_name_handle(actele->e.s8->alfa.name,&dims[0],&dims[1],&dims[2],&handles[1],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->alfa.fdim != dims[0] ||
       actele->e.s8->alfa.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*-------------------------------------------------------- read alfa */
   pss_read_array_name_handle(actele->e.s8->alfa.name,&(actele->e.s8->alfa),&handles[1],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Dtildinv*/
   pss_getdims_name_handle(actele->e.s8->Dtildinv.name,&dims[0],&dims[1],&dims[2],&handles[2],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Dtildinv.fdim != dims[0] ||
       actele->e.s8->Dtildinv.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------- read Dtildinv */
   pss_read_array_name_handle(actele->e.s8->Dtildinv.name,&(actele->e.s8->Dtildinv),&handles[2],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Lt*/
   pss_getdims_name_handle(actele->e.s8->Lt.name,&dims[0],&dims[1],&dims[2],&handles[3],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Lt.fdim != dims[0] ||
       actele->e.s8->Lt.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------------- read Lt */
   pss_read_array_name_handle(actele->e.s8->Lt.name,&(actele->e.s8->Lt),&handles[3],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------------------------------- check dimensions Rtilde*/
   pss_getdims_name_handle(actele->e.s8->Rtilde.name,&dims[0],&dims[1],&dims[2],&handles[4],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s8->Rtilde.fdim != dims[0] ||
       actele->e.s8->Rtilde.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*------------------------------------------------------ read Rtilde */
   pss_read_array_name_handle(actele->e.s8->Rtilde.name,&(actele->e.s8->Rtilde),&handles[4],&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_read_restart */
#endif
