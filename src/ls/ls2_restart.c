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
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"
#include "ls.h"
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


extern struct _LS_POLY_DATA_ARRAY    lspolyarr;


/*!---------------------------------------------------------------------
\brief write element restart data of ls2 element

<pre>                                                         genk 09/03

   at the moment only the curvature at the free surface has to be
   written to the binary file

</pre>
\param  *actele	   ELEMENT	   (i)    actual element
\param   nhandle   INT             (i)    number of handles
\param  *handles   long int        (i)    handles

\return void

------------------------------------------------------------------------*/
void ls2_write_restart(
  ELEMENT *actele,
  INT nhandle,
  long int *handles
  )
{
  INT ierr;
  FILE *out;
  ARRAY            arr1_;
  INT             *arr1;

#ifdef DEBUG
  dstrc_enter("ls2_write_restart");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  arr1 = amdef("...", &arr1_,4,1, "IV");

  /* */
  out = allfiles.out_pss;
  /* */
  handles[0]=0;
  handles[0]+=1;
  if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
  /* write nlayer */
  arr1[0] = actele->e.ls2->nlayer;
  if (actele->e.ls2->my_fluid!=NULL)
  {
    arr1[1] = actele->e.ls2->my_fluid->mat;
  }
  else
  {
    arr1[1] = 0;
  }
  arr1[2] = actele->e.ls2->is_el_first;
  arr1[3] = actele->e.ls2->is_el_last;
  pss_write_array(&(arr1_),&(handles[1]),out,&ierr);
  if (ierr != 1) dserror("Error writing restart data");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of ls2_write_restart */

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
void ls2_read_restart(
  ELEMENT *actele,
  INT nhandle,
  long int *handles
  )
{
  INT ierr;
  INT dims[3];
  FILE *in;
  ARRAY            arr1_;
  INT             *arr1;

#ifdef DEBUG
  dstrc_enter("ls2_read_restart");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  arr1 = amdef("...", &arr1_,4,1, "IV");

  /* */
  in = allfiles.in_pss;
  /* */
  if (!in) dserror("There is no restart input file open");

  /* check dimensions nlayer */
  pss_getdims_name_handle(arr1_.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
  if (ierr != 1) dserror("Cannot read restart data");
  if (arr1_.fdim != dims[0] || arr1_.sdim != dims[1])
    dserror("Mismatch in reading element restart data");
  /* read nlayer */
  pss_read_array_name_handle(arr1_.name,&(arr1_),&handles[1],in,&ierr);
  if (ierr != 1) dserror("Cannot read restart data");
  /* set nlayer */
  actele->e.ls2->nlayer = arr1[0];
  if (actele->e.ls2->my_fluid!=NULL)
  {
    /* set material index of associated fluid element */
    actele->e.ls2->my_fluid->mat = arr1[1];
  }
  else
  {
    /* */
  }
  /* set is_el_first and is_el_last flags */
  actele->e.ls2->is_el_first = arr1[2];
  actele->e.ls2->is_el_last = arr1[3];

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f2_read_restart */

#endif
