/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_write_restart:  writes all the element data that is needed for a restart
 - s9_read_restart:   reads all the element data that is needed for a restart

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

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


/*!----------------------------------------------------------------------
\brief write the data needed to restart this element                                       

<pre>                     m.gee 5/02              modified by    sh 02/03
This routine writes the data needed to restart the shell9 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9_write_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT i,ierr;
INT num_klay = actele->e.s9->num_klay;
FILE *out;
#ifdef DEBUG 
dstrc_enter("s9_write_restart");
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
if (actele->e.s9->nhyb)
{
   handles[0]+=4;
   if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
   /*------------------------------------------------------- write alfa */
   pss_write_array(&(actele->e.s9->alfa),&(handles[1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*------------------------ write EAS-Arrays for each kinematic layer */
   for (i=0; i<num_klay; i++)
   {
      /*------------------------------------------------ write Dtildinv */
      pss_write_array(&(actele->e.s9->Dtildinv[i]),&(handles[2]),out,&ierr);
      if (ierr != 1) dserror("Error writing restart data");
      /*------------------------------------------------------ write Lt */
      pss_write_array(&(actele->e.s9->Lt[i]),&(handles[3]),out,&ierr);
      if (ierr != 1) dserror("Error writing restart data");
      /*-------------------------------------------------- write Rtilde */
      pss_write_array(&(actele->e.s9->Rtilde[i]),&(handles[4]),out,&ierr);
      if (ierr != 1) dserror("Error writing restart data");
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_write_restart */

/*!----------------------------------------------------------------------
\brief read the data needed to restart this element                                       

<pre>                     m.gee 5/02              modified by    sh 02/03
This routine reads the data needed to restart the shell9 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/02    |
 | read the data needed to restart this element                         |
 | modified from shell8                                    sh 02/03     |
 *----------------------------------------------------------------------*/
void s9_read_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT i,ierr;
INT num_klay = actele->e.s9->num_klay;
INT dims[3];
FILE *in;
#ifdef DEBUG 
dstrc_enter("s9_read_restart");
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
if (actele->e.s9->nhyb)
{
   /*------------------------------------------------- check dimensions alfa*/
   pss_getdims_name_handle(actele->e.s9->alfa.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->alfa.fdim != dims[0] ||
       actele->e.s9->alfa.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*-------------------------------------------------------- read alfa */
   pss_read_array_name_handle(actele->e.s9->alfa.name,&(actele->e.s9->alfa),&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
   /*------------------------- read EAS-Arrays for each kinematic layer */
   for (i=0; i<num_klay; i++)
   {
      /*------------------------------------------------- check dimensions Dtildinv*/
      pss_getdims_name_handle(actele->e.s9->Dtildinv[i].name,&dims[0],&dims[1],&dims[2],&handles[2],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");
      if (actele->e.s9->Dtildinv[i].fdim != dims[0] ||
          actele->e.s9->Dtildinv[i].sdim != dims[1])
          dserror("Mismatch in reading element restart data");
      /*---------------------------------------------------- read Dtildinv */
      pss_read_array_name_handle(actele->e.s9->Dtildinv[i].name,&(actele->e.s9->Dtildinv[i]),&handles[2],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");  
      /*------------------------------------------------- check dimensions Lt*/
      pss_getdims_name_handle(actele->e.s9->Lt[i].name,&dims[0],&dims[1],&dims[2],&handles[3],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");
      if (actele->e.s9->Lt[i].fdim != dims[0] ||
          actele->e.s9->Lt[i].sdim != dims[1])
          dserror("Mismatch in reading element restart data");
      /*---------------------------------------------------------- read Lt */
      pss_read_array_name_handle(actele->e.s9->Lt[i].name,&(actele->e.s9->Lt[i]),&handles[3],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");  
      /*------------------------------------------------- check dimensions Rtilde*/
      pss_getdims_name_handle(actele->e.s9->Rtilde[i].name,&dims[0],&dims[1],&dims[2],&handles[4],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");
      if (actele->e.s9->Rtilde[i].fdim != dims[0] ||
          actele->e.s9->Rtilde[i].sdim != dims[1])
          dserror("Mismatch in reading element restart data");
      /*------------------------------------------------------ read Rtilde */
      pss_read_array_name_handle(actele->e.s9->Rtilde[i].name,&(actele->e.s9->Rtilde[i]),&handles[4],in,&ierr);
      if (ierr != 1) dserror("Cannot read restart data");  
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_read_restart */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
