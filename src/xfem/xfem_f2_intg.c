#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "xfem_prototypes.h"





/************************************************************************
------------------------------------------ last checked by Irhan 15.04.04
************************************************************************/
void xfem_f2_intg(
  FLUID_DATA    *data
  )
{
  DOUBLE     DUMMY = 0.0;
  
#ifdef DEBUG 
  dstrc_enter("xfem_f2_intg");
#endif
/*----------------------------------------------------------------------*/

  /* QUARILATERALS */
  /* weights */                                 
  data->qwgt[0][0] = 2.0000000000000;
  data->qwgt[0][1] = 1.0000000000000;
  data->qwgt[1][1] = 1.0000000000000;  
  data->qwgt[0][2] = 0.5555555555555;
  data->qwgt[1][2] = 0.8888888888888;
  data->qwgt[2][2] = 0.5555555555555;  

  /* coordinates for  gauss points */
  data->qxg[0][0] =  0.0000000000000;
  data->qxg[0][1] = -0.5773502691896;                                            
  data->qxg[1][1] =  0.5773502691896;
  data->qxg[0][2] = -0.7745966692414;
  data->qxg[1][2] =  0.0000000000000;
  data->qxg[2][2] =  0.7745966692414;

  /* TRIANGLES */
  /* weights */
  data->twgt[0][0] = 0.5000000000000;
  data->twgt[0][1] = DUMMY;
  data->twgt[1][1] = DUMMY;  
  data->twgt[0][2] = 0.1666666666667;
  data->twgt[1][2] = 0.1666666666667;
  data->twgt[2][2] = 0.1666666666667;
  /* Gauss points */
  data->txgr[0][0] = 0.3333333333333;
  data->txgr[0][1] = DUMMY;
  data->txgr[1][1] = DUMMY;
  data->txgr[0][2] = 0.6666666666667;
  data->txgr[1][2] = 0.1666666666667;
  data->txgr[2][2] = 0.1666666666667;

  data->txgs[0][0] = 0.3333333333333;
  data->txgs[0][1] = DUMMY;
  data->txgs[1][1] = DUMMY;
  data->txgs[0][2] = 0.1666666666667;
  data->txgs[1][2] = 0.6666666666667;
  data->txgs[2][2] = 0.1666666666667;
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of xfem_f2_intg */
#endif
