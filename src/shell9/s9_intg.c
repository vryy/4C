/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9intg: which gets the natural coordinates of the integration points
           and their weights


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief integration points and weights

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine gets the natural coordinates of the integration points and
their weights for a numerical integration
</pre>
\param  const ELEMENT *ele    (i) element array of actual element
\param  S9_DATA       *data   (o) coordinates and weights at GP
\param  INT            option (i) ?

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9eleload()     [s9_load1.c]
                             s9static_keug() [s9_static_keug.c]
                             s9_stress()     [s9_stress.c]

*----------------------------------------------------------------------*/
void s9intg(const ELEMENT   *ele,
            S9_DATA         *data,
            INT              option)
{
DOUBLE b,wgt,wgt0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9intg");
#endif
/*----------------------------------------------------------------------*/
if (option==0)
{
       switch(ele->e.s9->nGP[2])/*---------------- thickness direction t */
       {
       case 1:
          data->xgpt[0] = 0.0;
          data->xgpt[1] = 0.0;
          data->xgpt[2] = 0.0;
          data->wgtt[0] = 2.0;
          data->wgtt[1] = 0.0;
          data->wgtt[2] = 0.0;
       break;
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
           dserror("'s9_intg.c': unknown number of GP points thickness direction t");
       break;
       }
   if (ele->distyp == quad4 || /*---------------- quadrilateral elements */
       ele->distyp == quad8 ||
       ele->distyp == quad9 )
   {
       switch(ele->e.s9->nGP[0])/* direction r */
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
          dserror("'s9_intg.c': unknown number of GP points direction r");
       break;
       }
       switch(ele->e.s9->nGP[1])/* direction s */
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
          dserror("'s9_intg.c': unknown number of GP points direction s");
       break;
       }
    }
    else if (ele->distyp == tri3 || /*-------------- triangular elements */
             ele->distyp == tri6 )
    {
       switch(ele->e.s9->nGP_tri)
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
             dserror("'s9_intg.c': unknown number of GP for triangular elements");
       break;
       }
    }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9intg */


/*!----------------------------------------------------------------------
\brief integration points for stress extrapolation

<pre>                                                            sh 07/03
This routine gets the natural coordinates of the integration points and
their weights for a numerical integration  -> for stress extrapolation
</pre>
\param  S9_DATA       *data   (o) coordinates and weights at GP
\param  INT            option (i) 4 (Q4) 9 (Q8/9)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_stress()     [s9_stress.c]

*----------------------------------------------------------------------*/
void s9intg_str(S9_DATA         *data,
                INT              option)
{
DOUBLE b;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9intg_str");
#endif
/*----------------------------------------------------------------------*/

/*---------------- thickness direction t is always 2 */
b   = 1.0/3.0;
b   = sqrt(b);
data->xgpt[0] = -b;
data->xgpt[1] =  b;
data->xgpt[2] =  0.0;

if (option == 4) /*2x2*/
{
   /* direction r */
   b   = 1.0/3.0;
   b   = sqrt(b);
   data->xgpr[0] = -b;
   data->xgpr[1] =  b;
   data->xgpr[2] =  0.0;

   /* direction s */
   b   = 1.0/3.0;
   b   = sqrt(b);
   data->xgps[0] = -b;
   data->xgps[1] =  b;
   data->xgps[2] =  0.0;
}
else if (option == 9) /*3x3*/
{
   /* direction r */
   b    = 3.0/5.0;
   b    = sqrt(b);
   data->xgpr[0] = -b;
   data->xgpr[1] =  0.0;
   data->xgpr[2] =  b;

   /* direction s */
   b    = 3.0/5.0;
   b    = sqrt(b);
   data->xgps[0] = -b;
   data->xgps[1] =  0.0;
   data->xgps[2] =  b;
}
else dserror("wrong option for s9_intg_str");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9intg_str */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/




