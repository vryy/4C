#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

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

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief writes the data needed to restart the beam element

<pre>                                                              fh 04/03
This routine writes all the data needed to restart the beam element

</pre>
\param *actele         ELEMENT    (i/o) actual element
\param nhandle         INT         (i)  number of handles
\param *handles        LONGINT     (i)  handles
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_write_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
FILE *out;
#ifdef DEBUG 
dstrc_enter("b3_write_restart");
#endif
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
handles[0]+=1;
if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
/*------------------------------------------------------ write elewa */
pss_write_array(&(actele->e.b3->elewa),&(handles[1]),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_write_restart */

/*!----------------------------------------------------------------------
\brief reads the data needed to restart the beam element

<pre>                                                              fh 04/03
This routine reads all the data needed to restart the beam element

</pre>
\param *actele         ELEMENT    (i/o) actual element
\param nhandle         INT         (i)  number of handles
\param *handles        LONGINT     (i)  handles
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/

void b3_read_restart(ELEMENT *actele, INT nhandle, long int *handles)
{
INT ierr;
INT dims[3];
FILE *in;
#ifdef DEBUG 
dstrc_enter("b3_read_restart");
#endif
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
/*----------------------------------------------------------------------*/
pss_getdims_name_handle(actele->e.b3->elewa.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);  
if (ierr != 1) dserror("Cannot read restart data");						   
if (actele->e.b3->elewa.fdim != dims[0] ||
    actele->e.b3->elewa.sdim != dims[1])
    dserror("Mismatch in reading element restart data");
pss_read_array_name_handle(actele->e.b3->elewa.name,&(actele->e.b3->elewa),&handles[1],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");                
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_read_restart */
#endif
/*! @} (documentation module close)*/
