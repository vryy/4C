#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"
/*----------------------------------------------------------------------*
 | integration points                                        mn 06/02   |
 -----------------------------------------------------------------------|
 | QUAD-ELEMENT                                                          |
 | COORDINATES AND WEIGHTING FACTORS OF GAUSS-INTEGRATION-POINTS FOR    |
 | NUMERICAL INTEGRATION                                                |
 *----------------------------------------------------------------------*/
void ale2_intg(const ELEMENT   *ele,
              ALE2_DATA        *data,
              int              option)
{
int i, k;

double zero  = 0.0;
double one   = 1.0;
double two   = 2.0;
double three = 3.0;

double  q12, q13, q16, q23;
double  xgr[13][8],xgs[13][8],wgtt[13][8];
static double xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("ale2_intg");
#endif
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    Q U A D                 ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
switch(ele->e.ale2->nGP[0])/* direction r */
{
case 2:
  data->xgpr[0] = -0.5773502691896;
  data->xgpr[1] =  0.5773502691896;
  
  data->wgtr[0] = 1.0            ;
  data->wgtr[1] = 1.0            ;
break;
case 1:
  data->xgpr[0] = 0.0;
  
  data->wgtr[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale2_intg");
break;
}
switch(ele->e.ale2->nGP[1])/* direction s */
{
case 2:
  data->xgps[0] = -0.5773502691896;
  data->xgps[1] =  0.5773502691896;
  
  data->wgts[0] = 1.0            ;
  data->wgts[1] = 1.0            ;
break;
case 1:
  data->xgps[0] = 0.0;
  
  data->wgts[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale2_intg");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2_intg */
/*----------------------------------------------------------------------*/
#endif
