/*-----------------------------------------------------------------------*/
/*!
\file
\brief write and read element restart data of FLUID3_FAST


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "f3f_prototypes.h"
#include "../fluid3/fluid3.h"

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



/*-----------------------------------------------------------------------*/
/*!
  \brief write element restart data of fast fluid3 element

  at the moment only the stab parameter of the last time step has to be
  written to the binary file

  \param actele    *ELEMENT  (i)    actual element
  \param nhandle    INT      (i)    number of handles
  \param handles   *long int (i)    handles

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3f_write_restart(
    ELEMENT   *actele,
    INT        nhandle,
    long int  *handles
    )
{

  INT ierr;
  INT counter=1;
  FLUID_DYNAMIC   *fdyn;
  FILE *out;

#ifdef DEBUG
  dstrc_enter("f3f_write_restart");
#endif


  /* fluid3_fast element has to write the stab parameter
     of the last time step.  */


  fdyn= alldyn[genprob.numff].fdyn;
  out = allfiles.out_pss;

  handles[0]=0;
  handles[0]+=1;

  if (handles[0]+1 > nhandle)
    dserror("Handle range too small for element");


  /* write stab-par */
  pss_write_array(&(actele->e.f3->tau_old),&(handles[counter]),out,&ierr);
  if (ierr != 1)
    dserror("Error writing restart data");
  counter++;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3f_write_restart */




/*-----------------------------------------------------------------------*/
/*!
  \brief read element restart data of fast fluid3 element

  at the moment only the stab parameter of the last time step has to be
  read to the binary file

  \param actele    *ELEMENT  (i)    actual element
  \param nhandle    INT      (i)    number of handles
  \param handles   *long int (i)    handles

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3f_read_restart(
    ELEMENT   *actele,
    INT        nhandle,
    long int  *handles
    )
{

  INT ierr;
  INT dims[3];
  INT          counter=1;
  FLUID_DYNAMIC   *fdyn;
  FILE *in;

#ifdef DEBUG
  dstrc_enter("f3f_read_restart");
#endif


  /* fluid3 element has to read the stab parameter
     of the last time step. */


  fdyn= alldyn[genprob.numff].fdyn;
  in  = allfiles.in_pss;
  if (!in)
    dserror("There is no restart input file open");


  /* check dimensions tau_old */
  pss_getdims_name_handle(actele->e.f3->tau_old.name,&dims[0],&dims[1],
      &dims[2],&handles[counter],in,&ierr);

  if (ierr != 1)
    dserror("Cannot read restart data");

  if (actele->e.f3->tau_old.fdim != dims[0] ||
      actele->e.f3->tau_old.sdim != dims[1])
    dserror("Mismatch in reading element restart data");


  /* read tau_old */
  pss_read_array_name_handle(actele->e.f3->tau_old.name,
      &(actele->e.f3->tau_old),&handles[counter],in,&ierr);
  if (ierr != 1)
    dserror("Cannot read restart data");
  counter++;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3f_read_restart */




#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


