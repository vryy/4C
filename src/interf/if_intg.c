/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'ifintg' which gives the coordinates and
weight factors for numerical integration of a interface element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h" 

/*! 
\addtogroup INTERF
*/
/*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief coordinates and weight factors for numerical integration 

<pre>                                                              ah 05/03 
This routine  gives the coordinates and weight factors for numerical
integration of a interface shell element. Gauss-Integration with 2 GP 
is used here.

</pre>
\param *data  IF_DATA  (I)   coordinates and weighting factors at gaussian points
\param *ele   ELEMENT  (I)   the element

\warning There is nothing special to this routine
\return void                                               
\sa calling:   ---; 
    called by: ifinit;

*----------------------------------------------------------------------*/
void ifintg(ELEMENT       *ele,
            INTERF_DATA   *data)
{
#ifdef DEBUG 
dstrc_enter("ifintg");
#endif


switch(ele->e.interf->nGP)/* direction s */
{
case 1:
   data->xgr[0]  = 0.0;  
   data->wgtr[1] = 2.0; 
break;
/*----------------------------------------------------------------------*/
case 2:
   data->xgr[0]  = - 1/sqrt(3.0);  
   data->xgr[1]  =   1/sqrt(3.0);   
   data->wgtr[0] = 1.0; 
   data->wgtr[1] = 1.0; 
break;
/*----------------------------------------------------------------------*/
case 3:
   data->xgr[0]  = - sqrt(3.0/5.0);  
   data->xgr[1]  =   0.0;   
   data->xgr[2]  =   sqrt(3.0/5.0);   
   data->wgtr[0] = 5.0/9.0; 
   data->wgtr[1] = 8.0/9.0; 
   data->wgtr[2] = 5.0/9.0; 
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ifintg */
/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
