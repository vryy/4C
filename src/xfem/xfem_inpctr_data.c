#include "../headers/standardtypes.h"
#include "xfem_prototypes.h"





/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES                 allfiles;
struct _XFEM_DATA                    xfem_data;





/*---------------------------------------------------------------------
\brief input of the XFEM DATA block in the input-file

In this routine the data in the XFEM DATA block of the input file
are read and stored in xfem_data
------------------------------------------------------------------------*/
void xfem_inpctr_data()
{
  INT      ierr;
  CHAR     buffer[50];

#ifdef DEBUG
  dstrc_enter("xfem_inpctr_data");
#endif
/*----------------------------------------------------------------------*/

  /* read data */
  if (frfind("-XFEM GENERAL DATA")==0) dserror("frfind: XFEM GENERAL DATA not in input file");
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read character */
    frchar("XFEM_ONOFF"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"on",2)==0)
        xfem_data.xfem_on_off=1;
      else
        xfem_data.xfem_on_off=0;
    }
    frchar("XFEM_OPTIMIZE"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"yes",3)==0)
        xfem_data.xfem_optimize=1;
      else
        xfem_data.xfem_optimize=0;
    }
    frchar("XFEM_AREACHECK"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"yes",3)==0)
        xfem_data.xfem_area_check=1;
      else
        xfem_data.xfem_area_check=0;
    }
    frread();
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of xfem_inpctr_data */
