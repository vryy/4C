/*!----------------------------------------------------------------------
\file
\brief contains routines to evaluate element quality measures for 2-D
ale elements

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
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
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      


/*!----------------------------------------------------------------------
\brief  minimum corner angle of 4 noded quad ale element

<pre>                                                             ck 06/03
This routine evaluates the maximal corner angle distortion 
of a 4 noded quad element, the resulting distortion parameter is normalised
such that 1 - element with for rectangular corners
          0 - failure element with one corner with 180 degree angle

</pre>
\param **xyz  DOUBLE    (i)   elemental coordinates

\warning There is nothing special to this routine
\return DOUBLE corner angle
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
\brief  minimum corner angle of 3 noded triangle ale element

<pre>                                                             ck 07/03
This routine evaluates the maximal corner angle distortion 
of a 3 noded triangle element, the resulting distortion parameter is 
normalised
such that 1 - element with for corners of 60 degree angle
          0 - failure element with one corner with 180 degree angle

</pre>
\param **xyz  DOUBLE    (i)   elemental coordinates

\warning There is nothing special to this routine
\return DOUBLE corner angle
\sa calling: ---;
             called by: ale2_statik_ke(), ale2_statik_ke_test2()

*----------------------------------------------------------------------*/
DOUBLE ale2_corner_angle_tria(DOUBLE **xyz)
{

DOUBLE edge[3];     /* lengths of the element edges */
DOUBLE delta_x[3];  /* x- and y-difference of nodes*/
DOUBLE delta_y[3];
DOUBLE win;         /* corner angle */
DOUBLE ca;          /* corner angle criterion */
DOUBLE max;         /* maximal angular distortion */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_corner_angle_tria");
#endif
/*--------------------------------------------- x- and y-differences ---*/
delta_x[0] = xyz[0][0]-xyz[1][0];   /* line 0-1 */
delta_y[0] = xyz[0][1]-xyz[1][1];
delta_x[1] = xyz[1][0]-xyz[2][0];   /* line 1-2 */
delta_y[1] = xyz[1][1]-xyz[2][1];
delta_x[2] = xyz[2][0]-xyz[0][0];   /* line 2-0 */
delta_y[2] = xyz[2][1]-xyz[0][1];
/*------------------------------------------- evaluate element edges ---*/
edge[0] = sqrt( delta_x[0] * delta_x[0] 
              + delta_y[0] * delta_y[0] ); /*line 0-1*/
edge[1] = sqrt( delta_x[1] * delta_x[1] 
              + delta_y[1] * delta_y[1] ); /*line 1-2*/
edge[2] = sqrt( delta_x[2] * delta_x[2] 
              + delta_y[2] * delta_y[2] ); /*line 2-0*/
/*----------------------------------------------- angular distortion ---*/
max = 0.0;
win = acos( (delta_x[2]*delta_x[0] + delta_y[2]*delta_y[0])/edge[2]/edge[0] );
ca = fabs(win-1.047197551) / 2.094395102;
max = ( ca > max ) ? ca : max;  /* angle at node 0 checked */
win = acos( (delta_x[0]*delta_x[1] + delta_y[0]*delta_y[1])/edge[0]/edge[1] );
ca = fabs(win-1.047197551) / 2.094395102;
max = ( ca > max ) ? ca : max;  /* angle at node 1 checked */
win = acos( (delta_x[1]*delta_x[2] + delta_y[1]*delta_y[2])/edge[1]/edge[2] );
ca = fabs(win-1.047197551) / 2.094395102;
max = ( ca > max ) ? ca : max;  /* angle at node 2 checked */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return (1.0-max);
} /* end of ale2_corner_angle-tria */



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
}

/*!----------------------------------------------------------------------
\brief  distortion measure by means of aspect ratio

<pre>                                                             ck 07/03
This routine evaluates the minimal aspect ratio of a 3 noded triangular element,
the ratio is evaluated between the largest inner and smallest outer circle.
Normalisation such that 1 stands for a regular triangle element and 
0 represents an element that is degenerated to a line.

</pre>
\param **xyz  DOUBLE    (i)   elemental coordinates

\warning There is nothing special to this routine
\return void                                               
\sa calling: 
             called by: ale2_statik_ke(), ale2_statik_ke_test2()

*----------------------------------------------------------------------*/
DOUBLE ale2_aspect_ratio_tria(DOUBLE **xyz)
{
INT i;             /* a counter */

DOUBLE ratio;      /* aspect ratio */

DOUBLE sqed[3];    /* edges squared */
DOUBLE edge[3];    /* lengths of edges */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_aspect_ratio");
#endif
/*-------------------------------------------------- squared lengths ---*/
sqed[0] = (xyz[0][0]-xyz[1][0])*(xyz[0][0]-xyz[1][0])
         +(xyz[0][1]-xyz[1][1])*(xyz[0][1]-xyz[1][1]); /* line 0-1 squared */
sqed[1] = (xyz[1][0]-xyz[2][0])*(xyz[1][0]-xyz[2][0])
         +(xyz[1][1]-xyz[2][1])*(xyz[1][1]-xyz[2][1]); /* line 1-2 squared */
sqed[2] = (xyz[2][0]-xyz[0][0])*(xyz[2][0]-xyz[0][0])
         +(xyz[2][1]-xyz[0][1])*(xyz[2][1]-xyz[0][1]); /* diag 2-0 squared */
/*---------------------------------------------------------- lengths ---*/
for (i=0; i<3; i++) edge[i] = sqrt(sqed[i]);
/*----------------------------------------------------- aspect ratio ---*/
ratio = (2.0*(sqed[0]*sqed[1] + sqed[1]*sqed[2] + sqed[2]*sqed[0])
            - sqed[0]*sqed[0] - sqed[1]*sqed[1] - sqed[2]*sqed[2] )
	 /( (edge[0]+edge[1]+edge[2])*edge[0]*edge[1]*edge[2] ); 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return ratio;
} /* end of ale2_aspect_ratio_tria */


/*!----------------------------------------------------------------------
\brief  writes element quality to ele->e.ale2->quality

<pre>                                                             ck 07/03
This routine writes the appropriate element quality to the ale element in
order to prepare the element quality statistics

</pre>
\param  *ele       ELEMENT   (i)   actual element
\param   quality   INT       (i)   flag, which quality to write
\param **xyz       DOUBLE    (i)   elemental coordinates
\param   min_detF  DOUBLE    (i)   minimal Jacobian determinant

\warning There is nothing special to this routine
\return void                                               
\sa calling: ale2_aspect_ratio(), ale2_aspect_ratio_tria(), 
             ale2_corner_angle(), ale2_corner_angle_tria(),
	     ale2_el_area();
             called by: ale2_static_ke_stiff(), ale2_static_ke_prestress(), 
	     ale2_static_ke_step2(), ale2_static_ke_spring(),
	     ale2_static_ke_laplace()

*----------------------------------------------------------------------*/
void write_element_quality(ELEMENT  *ele, 
                           INT       quality, 
			   DOUBLE  **xyz, 
			   DOUBLE    min_detF)
{
DOUBLE el_area;            /* area of the actual element */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale_quality");
#endif
/*----------------------------------------------------------------------*/
if (quality == 0)       /* no quality monitoring */
   ele->e.ale2->quality = 0.0; 
else if (quality == 1)  /* case aspect ratio */
   switch (ele->distyp)
   {
      case quad4:
         ele->e.ale2->quality = ale2_aspect_ratio(xyz); 
	 break;
      case tri3:
         ele->e.ale2->quality = ale2_aspect_ratio_tria(xyz); 
	 break;
      default:
         dswarning(1,2);
	 break;
   }
else if (quality == 2)  /* case corner angle */
   switch (ele->distyp)
   {
      case quad4:
         ele->e.ale2->quality = ale2_corner_angle(xyz); 
	 break;
      case tri3:
         ele->e.ale2->quality = ale2_corner_angle_tria(xyz); 
	 break;
      default:
         dswarning(1,3);
	 break;
   }
else if (quality == 3)  /* case (normalised) min. Jacobian determinant */
{
   switch (ele->distyp)
   {
       case quad4:
          el_area = ale2_el_area(xyz);
	  ele->e.ale2->quality = min_detF * 4.0/el_area;
	  break;
       case tri3:
          dswarning(1,1);
	  break;
       default:
	  dswarning(1,4);
	  break;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of write_element_quality */


/*!----------------------------------------------------------------------
\brief  element quality statistics

<pre>                                                             ck 06/03
element quality statistics

</pre>
\param *actfield  FIELD     (i)   actual field
\param  step      INT       (i)   actual time step
\param *actintra  INTRA     (i)   intra-communicator
\param *actpart   PARTITION (i)   actual partition

\warning Works at the moment for ale2-elements only!
\return void                                               
\sa calling: 
             called by: dyn_ale();

*----------------------------------------------------------------------*/
void ale_quality(FIELD *field,INT step, 
                 INTRA  *actintra, PARTITION    *actpart)
{
INT i;      /* a counter */
INT numel;  /* number of elements in this discretisation */
INT numele_total;

DOUBLE quality;     /* current element quality measure */
DOUBLE square = 0;       
DOUBLE min, max;    /* minimal and maximal quality */
DOUBLE stand_degr;  /* standard degression*/
DOUBLE average;
DOUBLE recv=ZERO;
ELEMENT *actele;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale_quality");
#endif
/*----------------------------------------------------------------------*/

numel = actpart->pdis[0].numele;
numele_total=field->dis[0].numele;
average = 0.0;
min = 2.0;
max = -1.0;

for (i=0; i<numel; i++)  /* loop over all elements */
{
   actele = actpart->pdis[0].element[i];
   if (actele->proc!=par.myrank) continue;
   quality = actele->e.ale2->quality;   
   average += quality;
   min = DMIN(quality,min);
   max = DMAX(quality,max);
   square += quality * quality;
}

#ifdef PARALLEL
MPI_Reduce(&square,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
square=recv;
MPI_Reduce(&average,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
average=recv;
MPI_Reduce(&min,&recv,1,MPI_DOUBLE,MPI_MIN,0,actintra->MPI_INTRA_COMM);
min=recv;
MPI_Reduce(&max,&recv,1,MPI_DOUBLE,MPI_MAX,0,actintra->MPI_INTRA_COMM);
max=recv;
#endif

if (par.myrank==0)
{
   if (numele_total > 1)
      stand_degr = 1.0/(numele_total-1.0) * ( square - 1.0/numele_total*average*average );
   else
      stand_degr = 0.0;
   average = average/numele_total;
/*----------------------------------------------------- gnuplot file ---*/
   fprintf(allfiles.gnu,"%i  %8.7f  %8.7f  %8.7f  %8.7f\n", 
           step-1, average, stand_degr, min, max);
   fflush(allfiles.gnu);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
}

/*! @} (documentation module close)*/
#endif
