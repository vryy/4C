/*!----------------------------------------------------------------------
\file
\brief contains routines to evaluate element quality measures for 2-D
ale elements

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  minimum corner angle of 4 noded quad ale element

<pre>                                                             ck 06/03
This routine evaluates the maximal corner angle distortion 
of a 4 noded quad element, the resulting distortion parameter is normalised
such that 1 - element with four rectangular corners
          0 - failure element with one corner with 180 degree angle

</pre>
\param **xyz  DOUBLE    (i)   elemental coordinates

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---;
             called by: ale2_statik_ke(), ale2_statik_ke_test2()

*----------------------------------------------------------------------*/
DOUBLE ale2_corner_angle(DOUBLE **xyz)
{

DOUBLE edge[4];     /* lengths of the element edges */
DOUBLE delta_x[4];  /* x- and y-difference of nodes*/
DOUBLE delta_y[4];
DOUBLE win;         /* corner angle */
DOUBLE ca;          /* corner angle criterion */
DOUBLE max;         /* maximal angular distortion */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_corner_angle");
#endif
/*--------------------------------------------- x- and y-differences ---*/
delta_x[0] = xyz[0][0]-xyz[1][0];   /* line 0-1 */
delta_y[0] = xyz[0][1]-xyz[1][1];
delta_x[1] = xyz[1][0]-xyz[2][0];   /* line 1-2 */
delta_y[1] = xyz[1][1]-xyz[2][1];
delta_x[2] = xyz[2][0]-xyz[3][0];   /* line 2-3 */
delta_y[2] = xyz[2][1]-xyz[3][1];
delta_x[3] = xyz[3][0]-xyz[0][0];   /* line 3-0 */
delta_y[3] = xyz[3][1]-xyz[0][1];
/*------------------------------------------- evaluate element edges ---*/
edge[0] = sqrt( delta_x[0] * delta_x[0] 
              + delta_y[0] * delta_y[0] ); /*line 0-1*/
edge[1] = sqrt( delta_x[1] * delta_x[1] 
              + delta_y[1] * delta_y[1] ); /*line 1-2*/
edge[2] = sqrt( delta_x[2] * delta_x[2] 
              + delta_y[2] * delta_y[2] ); /*line 2-3*/
edge[3] = sqrt( delta_x[3] * delta_x[3] 
              + delta_y[3] * delta_y[3] ); /*line 3-0*/
/*----------------------------------------------- angular distortion ---*/
max = 0.0;
win = acos( (delta_x[3]*delta_x[0] + delta_y[3]*delta_y[0])/edge[3]/edge[0] );
ca = fabs(win-1.570796327) / 1.570796327;
max = ( ca > max ) ? ca : max;  /* angle at node 0 checked */
win = acos( (delta_x[0]*delta_x[1] + delta_y[0]*delta_y[1])/edge[0]/edge[1] );
ca = fabs(win-1.570796327) / 1.570796327;
max = ( ca > max ) ? ca : max;  /* angle at node 1 checked */
win = acos( (delta_x[1]*delta_x[2] + delta_y[1]*delta_y[2])/edge[1]/edge[2] );
ca = fabs(win-1.570796327) / 1.570796327;
max = ( ca > max ) ? ca : max;  /* angle at node 2 checked */
win = acos( (delta_x[2]*delta_x[3] + delta_y[2]*delta_y[3])/edge[2]/edge[3] );
ca = fabs(win-1.570796327) / 1.570796327;
max = ( ca > max ) ? ca : max;  /* angle at node 3 checked */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return (1.0-max);
}



/*!----------------------------------------------------------------------
\brief  distortion measure by means of aspect ratio

<pre>                                                             ck 06/03
This routine evaluates the minimal aspect ratio of a 4 noded quad element,
the ratio is evaluated between the largest inner and smallest outer circle
of each of the four participating triangles. Normalisation such that 1 stands
for a square element and 0 represents an element that is degenerated to a
triangle.

</pre>
\param **xyz  DOUBLE    (i)   elemental coordinates

\warning There is nothing special to this routine
\return void                                               
\sa calling: 
             called by: ale2_statik_ke(), ale2_statik_ke_test2()

*----------------------------------------------------------------------*/
DOUBLE ale2_aspect_ratio(DOUBLE **xyz)
{
INT i;             /* a counter */

DOUBLE ratio;      /* acutal aspect ratio */
DOUBLE min=2.0;    /* minimum aspect ratio */
DOUBLE factor = 1.0/(sqrt(2)-1.0)/2.0;

DOUBLE sqed[4];    /* edges squared */
DOUBLE sqdg[2];    /* diagonals squared */
DOUBLE edge[4];    /* lengths of edges */
DOUBLE diag[2];    /* lengths of diagonals*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_aspect_ratio");
#endif
/*-------------------------------------------------- squared lengths ---*/
sqed[0] = (xyz[0][0]-xyz[1][0])*(xyz[0][0]-xyz[1][0])
         +(xyz[0][1]-xyz[1][1])*(xyz[0][1]-xyz[1][1]); /* line 0-1 squared */
sqed[1] = (xyz[1][0]-xyz[2][0])*(xyz[1][0]-xyz[2][0])
         +(xyz[1][1]-xyz[2][1])*(xyz[1][1]-xyz[2][1]); /* line 1-2 squared */
sqed[2] = (xyz[2][0]-xyz[3][0])*(xyz[2][0]-xyz[3][0])
         +(xyz[2][1]-xyz[3][1])*(xyz[2][1]-xyz[3][1]); /* line 2-3 squared */
sqed[3] = (xyz[3][0]-xyz[0][0])*(xyz[3][0]-xyz[0][0])
         +(xyz[3][1]-xyz[0][1])*(xyz[3][1]-xyz[0][1]); /* line 3-0 squared */
sqdg[0] = (xyz[2][0]-xyz[0][0])*(xyz[2][0]-xyz[0][0])
         +(xyz[2][1]-xyz[0][1])*(xyz[2][1]-xyz[0][1]); /* diag 2-0 squared */
sqdg[1] = (xyz[1][0]-xyz[3][0])*(xyz[1][0]-xyz[3][0])
         +(xyz[1][1]-xyz[3][1])*(xyz[1][1]-xyz[3][1]); /* diag 1-3 squared */
/*---------------------------------------------------------- lengths ---*/
for (i=0; i<4; i++) edge[i] = sqrt(sqed[i]);
for (i=0; i<2; i++) diag[i] = sqrt(sqdg[i]);
/*----------------------------------------------------- aspect ratio ---*/
ratio = (2.0*(sqed[0]*sqed[3] + sqed[3]*sqdg[1] + sqdg[1]*sqed[0])
            - sqed[0]*sqed[0] - sqed[3]*sqed[3] - sqdg[1]*sqdg[1] )
	 /( (edge[0]+edge[3]+diag[1])*edge[0]*edge[3]*diag[1] )*factor;
min = (ratio < min) ? ratio : min;
ratio = (2.0*(sqed[0]*sqed[1] + sqed[1]*sqdg[0] + sqdg[0]*sqed[0])
            - sqed[0]*sqed[0] - sqed[1]*sqed[1] - sqdg[0]*sqdg[0] )
	 /( (edge[0]+edge[1]+diag[0])*edge[0]*edge[1]*diag[0] )*factor;
min = (ratio < min) ? ratio : min;
ratio = (2.0*(sqed[1]*sqed[2] + sqed[2]*sqdg[1] + sqdg[1]*sqed[1])
            - sqed[1]*sqed[1] - sqed[2]*sqed[2] - sqdg[1]*sqdg[1] )
	 /( (edge[1]+edge[2]+diag[1])*edge[1]*edge[2]*diag[1] )*factor;
min = (ratio < min) ? ratio : min;
ratio = (2.0*(sqed[2]*sqed[3] + sqed[3]*sqdg[0] + sqdg[0]*sqed[2])
            - sqed[2]*sqed[2] - sqed[3]*sqed[3] - sqdg[0]*sqdg[0] )
	 /( (edge[2]+edge[3]+diag[0])*edge[2]*edge[3]*diag[0] )*factor;
min = (ratio < min) ? ratio : min;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return min;
}/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*---------------------------------------------------------------------*/



/*!----------------------------------------------------------------------
\brief  element quality statistics

<pre>                                                             ck 06/03
element quality statistics

</pre>
\param *actfield  FIELD    (i)   actual field
\param  step      INT      (i)   actual time step

\warning Works at the moment for ale2-elements only!
\return void                                               
\sa calling: 
             called by: dyn_ale();

*----------------------------------------------------------------------*/
void ale_quality(FIELD *field,INT step)
{
INT i;      /* a counter */
INT numel;  /* number of elements in this discretisation */

DOUBLE quality;     /* current element quality measure */
DOUBLE square = 0;       
DOUBLE min, max;    /* minimal and maximal quality */
DOUBLE stand_degr;  /* standard degression*/
DOUBLE average;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale_quality");
#endif
/*----------------------------------------------------------------------*/

numel = field->dis->numele;
average = 0.0;
min = 2.0;
max = -1.0;

for (i=0; i<numel; i++)  /* loop over all elements */
{
   quality = field->dis->element[i].e.ale2->quality;
   average += quality;
   min = (quality < min) ? quality : min;
   max = (quality > max) ? quality : max;
   square += quality * quality;
}
if (numel > 1)
   stand_degr = 1.0/(numel-1.0) * ( square - 1.0/numel*average*average );
else
   stand_degr = 0.0;
average = average/numel;

/*----------------------------------------------------- gnuplot file ---*/
fprintf(allfiles.gnu,"%i  %8.7f  %8.7f  %8.7f  %8.7f\n", 
        step-1, average, stand_degr, min, max);
fflush(allfiles.gnu);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
}

/*! @} (documentation module close)*/
#endif
