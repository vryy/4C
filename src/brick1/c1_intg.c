/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1intg' which calclate data of 
       integration points for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief compute coordinates and weighting factors of gauss-integration-points

<pre>                                                              al 06/02
This routine calcuates coordinates and weighting factors of
gauss-integration-points 
for an 3D-hex-element.

</pre>
\param  *ele     ELEMENT (i)   the element
\param *data     C1_DATA (i)   structure containing gaussian point and weight

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1intg(ELEMENT         *ele,
            C1_DATA         *data)
{

DOUBLE zero  = 0.0;
DOUBLE one   = 1.0;
DOUBLE two   = 2.0;
DOUBLE three = 3.0;

static DOUBLE xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("c1intg");
#endif
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L     ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
switch(ele->e.c1->nGP[0])/* direction r */
{
case 2:
  data->xgrr[0] = -0.5773502691896;
  data->xgrr[1] =  0.5773502691896;
  
  data->wgtr[0] = 1.0            ;
  data->wgtr[1] = 1.0            ;
break;
case 3:
  data->xgrr[0] = -0.7745966692415;
  data->xgrr[1] =  0.0;
  data->xgrr[2] =  0.7745966692415;
  
  data->wgtr[0] =  0.5555555555556;
  data->wgtr[1] =  0.8888888888889;
  data->wgtr[2] =  0.5555555555556;
break;
case 4:
  data->xgrr[0] =  -0.8611363115941;
  data->xgrr[1] =  -0.3399810435849;
  data->xgrr[2] =   0.3399810435849;
  data->xgrr[3] =   0.8611363115941;
  
  data->wgtr[0] =  0.3478548451375;
  data->wgtr[1] =  0.6521451548625;
  data->wgtr[2] =  0.6521451548625;
  data->wgtr[3] =  0.3478548451375;
break;
case 5:
  data->xgrr[0] =  -0.9061798459387;
  data->xgrr[1] =  -0.5384693101057;
  data->xgrr[2] =   0.0;
  data->xgrr[3] =   0.5384693101057;
  data->xgrr[4] =   0.9061798459387;
  
  data->wgtr[0] =  0.2369268850562;
  data->wgtr[1] =  0.4786286704994;
  data->wgtr[2] =  0.5688888888889;
  data->wgtr[3] =  0.4786286704994;
  data->wgtr[4] =  0.2369268850562;
break;
case 6:
  data->xgrr[0] =  -0.9324695142032;
  data->xgrr[1] =  -0.6612093864663;
  data->xgrr[2] =  -0.2386191860832;
  data->xgrr[3] =   0.2386191860832;
  data->xgrr[4] =   0.6612093864663;
  data->xgrr[5] =   0.9324695142032;
  
  data->wgtr[0] =  0.1713244923792;
  data->wgtr[1] =  0.3607615730481;
  data->wgtr[2] =  0.4679139345727;
  data->wgtr[3] =  0.4679139345727;
  data->wgtr[4] =  0.3607615730481;
  data->wgtr[5] =  0.1713244923792;
break;
default:
   dserror("unknown number of gaussian points in c1intg");
break;
}
switch(ele->e.c1->nGP[1])/* direction s */
{
case 2:
  data->xgss[0] = -0.5773502691896;
  data->xgss[1] =  0.5773502691896;
  
  data->wgts[0] = 1.0            ;
  data->wgts[1] = 1.0            ;
break;
case 3:
  data->xgss[0] = -0.7745966692415;
  data->xgss[1] =  0.0;
  data->xgss[2] =  0.7745966692415;
  
  data->wgts[0] =  0.5555555555556;
  data->wgts[1] =  0.8888888888889;
  data->wgts[2] =  0.5555555555556;
break;
case 4:
  data->xgss[0] =  -0.8611363115941;
  data->xgss[1] =  -0.3399810435849;
  data->xgss[2] =   0.3399810435849;
  data->xgss[3] =   0.8611363115941;
  
  data->wgts[0] =  0.3478548451375;
  data->wgts[1] =  0.6521451548625;
  data->wgts[2] =  0.6521451548625;
  data->wgts[3] =  0.3478548451375;
break;
case 5:
  data->xgss[0] = -0.9061798459387;
  data->xgss[1] = -0.5384693101057;
  data->xgss[2] =  0.0;
  data->xgss[3] =  0.5384693101057;
  data->xgss[4] =  0.9061798459387;
  
  data->wgts[0] =  0.2369268850562;
  data->wgts[1] =  0.4786286704994;
  data->wgts[2] =  0.5688888888889;
  data->wgts[3] =  0.4786286704994;
  data->wgts[4] =  0.2369268850562;
break;
case 6:
  data->xgss[0] = -0.9324695142032;
  data->xgss[1] = -0.6612093864663;
  data->xgss[2] = -0.2386191860832;
  data->xgss[3] =  0.2386191860832;
  data->xgss[4] =  0.6612093864663;
  data->xgss[5] =  0.9324695142032;
  
  data->wgts[0] =  0.1713244923792;
  data->wgts[1] =  0.3607615730481;
  data->wgts[2] =  0.4679139345727;
  data->wgts[3] =  0.4679139345727;
  data->wgts[4] =  0.3607615730481;
  data->wgts[5] =  0.1713244923792;
break;
default:
   dserror("unknown number of gaussian points in c1intg");
break;
}
switch(ele->e.c1->nGP[2])/* direction t */
{
case 2:
  data->xgtt[0] = -0.5773502691896;
  data->xgtt[1] =  0.5773502691896;
  
  data->wgtt[0] = 1.0            ;
  data->wgtt[1] = 1.0            ;
break;
case 3:
  data->xgtt[0] = -0.7745966692415;
  data->xgtt[1] =  0.0;
  data->xgtt[2] =  0.7745966692415;
  
  data->wgtt[0] =  0.5555555555556;
  data->wgtt[1] =  0.8888888888889;
  data->wgtt[2] =  0.5555555555556;
break;
case 4:
  data->xgtt[0] =  -0.8611363115941;
  data->xgtt[1] =  -0.3399810435849;
  data->xgtt[2] =   0.3399810435849;
  data->xgtt[3] =   0.8611363115941;
  
  data->wgtt[0] =  0.3478548451375;
  data->wgtt[1] =  0.6521451548625;
  data->wgtt[2] =  0.6521451548625;
  data->wgtt[3] =  0.3478548451375;
break;
case 5:
  data->xgtt[0] = -0.9061798459387;
  data->xgtt[1] = -0.5384693101057;
  data->xgtt[2] =  0.0;
  data->xgtt[3] =  0.5384693101057;
  data->xgtt[4] =  0.9061798459387;
  
  data->wgtt[0] =  0.2369268850562;
  data->wgtt[1] =  0.4786286704994;
  data->wgtt[2] =  0.5688888888889;
  data->wgtt[3] =  0.4786286704994;
  data->wgtt[4] =  0.2369268850562;
break;
case 6:
  data->xgtt[0] = -0.9324695142032;
  data->xgtt[1] = -0.6612093864663;
  data->xgtt[2] = -0.2386191860832;
  data->xgtt[3] =  0.2386191860832;
  data->xgtt[4] =  0.6612093864663;
  data->xgtt[5] =  0.9324695142032;
  
  data->wgtt[0] =  0.1713244923792;
  data->wgtt[1] =  0.3607615730481;
  data->wgtt[2] =  0.4679139345727;
  data->wgtt[3] =  0.4679139345727;
  data->wgtt[4] =  0.3607615730481;
  data->wgtt[5] =  0.1713244923792;
break;
default:
   dserror("unknown number of gaussian points in c1intg");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1intg */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
