/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | integration points                                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8intg(const ELEMENT   *ele,
               S8_DATA         *data,
               INT              option)
{
DOUBLE b,wgt,wgt0;
#ifdef DEBUG
dstrc_enter("s8intg");
#endif
/*----------------------------------------------------------------------*/
if (option==0)
{
       switch(ele->e.s8->nGP[2])/*---------------- thickness direction t */
       {
          case 2:
             b   = 1.0/3.0;
             b   = sqrt(b);
             wgt = 1.0;
             data->xgpt[0] = -b;
             data->xgpt[1] =  b;
             data->xgpt[2] =  0.0;
             data->wgtt[0] =  wgt;
             data->wgtt[1] =  wgt;
             data->wgtt[2] =  0.0;
          break;
          default:
              dserror("unknown number of GP points");
          break;
       }
   if (ele->distyp == quad4 || /*---------------- quadrilateral elements */
       ele->distyp == quad8 ||
       ele->distyp == quad9 )
   {
       switch(ele->e.s8->nGP[0])/* direction r */
       {
       case 1:
          data->xgpr[0] = 0.0;
          data->xgpr[1] = 0.0;
          data->xgpr[2] = 0.0;
          data->wgtr[0] = 2.0;
          data->wgtr[1] = 0.0;
          data->wgtr[2] = 0.0;
       break;
       case 2:
          b   = 1.0/3.0;
          b   = sqrt(b);
          wgt = 1.0;
          data->xgpr[0] = -b;
          data->xgpr[1] =  b;
          data->xgpr[2] =  0.0;
          data->wgtr[0] =  wgt;
          data->wgtr[1] =  wgt;
          data->wgtr[2] =  0.0;
       break;
       case 3:
          b    = 3.0/5.0;
          b    = sqrt(b);
          wgt  = 5.0/9.0;
          wgt0 = 8.0/9.0;
          data->xgpr[0] = -b;
          data->xgpr[1] =  0.0;
          data->xgpr[2] =  b;
          data->wgtr[0] =  wgt;
          data->wgtr[1] =  wgt0;
          data->wgtr[2] =  wgt;
       break;
       default:
          dserror("unknown number of GP points");
       break;
       }
       switch(ele->e.s8->nGP[1])/* direction s */
       {
       case 1:
          data->xgps[0] = 0.0;
          data->xgps[1] = 0.0;
          data->xgps[2] = 0.0;
          data->wgts[0] = 2.0;
          data->wgts[1] = 0.0;
          data->wgts[2] = 0.0;
       break;
       case 2:
          b   = 1.0/3.0;
          b   = sqrt(b);
          wgt = 1.0;
          data->xgps[0] = -b;
          data->xgps[1] =  b;
          data->xgps[2] =  0.0;
          data->wgts[0] =  wgt;
          data->wgts[1] =  wgt;
          data->wgts[2] =  0.0;
       break;
       case 3:
          b    = 3.0/5.0;
          b    = sqrt(b);
          wgt  = 5.0/9.0;
          wgt0 = 8.0/9.0;
          data->xgps[0] = -b;
          data->xgps[1] =  0.0;
          data->xgps[2] =  b;
          data->wgts[0] =  wgt;
          data->wgts[1] =  wgt0;
          data->wgts[2] =  wgt;
       break;
       default:
          dserror("unknown number of GP points");
       break;
       }
    }
    else if (ele->distyp == tri3 || /*-------------- triangular elements */
             ele->distyp == tri6 )
    {
       switch(ele->e.s8->nGP_tri)
       {
       case 1:
          b   = 1.0/3,0;
          wgt = 1.0/2.0;
          data->xgpr[0] =  b;
          data->xgpr[1] =  0.0;
          data->xgpr[2] =  0.0;
          data->xgps[0] =  b;
          data->xgps[1] =  0.0;
          data->xgps[2] =  0.0;
          data->wgtr[0] =  wgt;
          data->wgtr[1] =  0.0;
          data->wgtr[2] =  0.0;
          data->wgts[0] =  wgt;
          data->wgts[1] =  0.0;
          data->wgts[2] =  0.0;
       break;
       case 3:
          b   = 1.0/2.0;
          wgt = 1.0/6.0;
          data->xgpr[0] =  b;
          data->xgpr[1] =  b;
          data->xgpr[2] =  0.0;
          data->xgps[0] =  0.0;
          data->xgps[1] =  b;
          data->xgps[2] =  b;
          data->wgtr[0] =  wgt;
          data->wgtr[1] =  wgt;
          data->wgtr[2] =  wgt;
          data->wgts[0] =  wgt;
          data->wgts[1] =  wgt;
          data->wgts[2] =  wgt;
       break;
       default:
             dserror("unknown number of GP points");
       break;
       }
    }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8intg */
#endif




