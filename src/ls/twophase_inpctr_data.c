#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"





/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES                 allfiles;
struct _TWOPHASE_DATA               *twophase_data;





/*---------------------------------------------------------------------
\brief input of the TWO PHASE FLOW DATA block in the input-file

In this routine the data in the TWO PHASE DATA block of the input file
are read and stored in xfem_data
------------------------------------------------------------------------*/
void twophase_inpctr_data()
{
  INT      ierr;
  CHAR     buffer[50];

#ifdef DEBUG
  dstrc_enter("twophase_inpctr_data");
#endif
/*----------------------------------------------------------------------*/

  /* allocate memory */
  twophase_data = (TWOPHASE_DATA*)CCACALLOC(1,sizeof(TWOPHASE_DATA));

  /* read data */
  if (frfind("-TWO PHASE FLOW GENERAL DATA")==0) dserror("frfind: TWO PHASE FLOW GENERAL DATA not in input file");
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frchar("SOLN_METHOD"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"lsxfem",6)==0)
        twophase_data->soln_method=tp_lsxfem;
      else if (strncmp(buffer,"smeared",7)==0)
        twophase_data->soln_method=tp_smeared;
      else
        twophase_data->soln_method=tp_none;
    }
    frread();
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of twophase_inpctr_data */
#endif /*  D_LS */
