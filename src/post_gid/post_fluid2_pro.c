
#ifdef D_FLUID2_PRO

#include "post_fluid2_pro.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void fluid2_pro_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("fluid2_pro_write_gauss");
#endif

  if (gid->is_fluid2_pro_22)
  {
    GiD_BeginGaussPoint(gid->fluid2_pro_22_name, GiD_Quadrilateral, gid->fluid2_pro_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_fluid2_pro_33)
  {
    GiD_BeginGaussPoint(gid->fluid2_pro_33_name, GiD_Quadrilateral, gid->fluid2_pro_33_name, 9, 0, 1);
    GiD_EndGaussPoint();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void fluid2_pro_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("fluid2_pro_write_mesh");
#endif

  if (gid->is_fluid2_pro_22)
  {
    dserror("fluid2_pro_22 mesh not supported");
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid2_pro_33)
  {
    dserror("fluid2_pro_33 mesh not supported");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
