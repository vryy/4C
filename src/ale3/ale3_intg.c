#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"
/*----------------------------------------------------------------------*
 | integration points                                        mn 06/02   |
 -----------------------------------------------------------------------|
 | HEX-ELEMENT                                                          |
 | COORDINATES AND WEIGHTING FACTORS OF GAUSS-INTEGRATION-POINTS FOR    |
 | NUMERICAL INTEGRATION                                                |
 *----------------------------------------------------------------------*/
void ale_intg(const ELEMENT   *ele,
              ALE3_DATA        *data,
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
dstrc_enter("ale_intg");
#endif
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L     ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
switch(ele->e.ale3->nGP[0])/* direction r */
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
   dserror("unknown number of gaussian points in ale_intg");
break;
}
switch(ele->e.ale3->nGP[1])/* direction s */
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
   dserror("unknown number of gaussian points in ale_intg");
break;
}
switch(ele->e.ale3->nGP[2])/* direction t */
{
case 2:
  data->xgpt[0] = -0.5773502691896;
  data->xgpt[1] =  0.5773502691896;
  
  data->wgtt[0] = 1.0            ;
  data->wgtt[1] = 1.0            ;
break;
case 1:
  data->xgpt[0] = 0.0;
  
  data->wgtt[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale_intg");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale_intg */
/*----------------------------------------------------------------------*/
#endif
