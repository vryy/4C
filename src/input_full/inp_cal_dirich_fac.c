/*!----------------------------------------------------------------------
\file
\brief contains the routine 'cal_dirich_fac' which calculates the value
       of a given spatial function for a given gnode

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/



#include "../headers/standardtypes.h"

/*!
  \addtogroup Ale
  *//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
  |                                                          mn 02/04    |
  | number of spatual load functions   numcurve                          |
  | vector of structures of functions                                    |
  | defined in input_curves.c                                            |
  | INT                 numfunct;                                        |
  | struct _FUNCT      *funct;                                           |
 *----------------------------------------------------------------------*/
extern INT                 numfunct;
extern struct _FUNCT      *funct;


/*!----------------------------------------------------------------------
  \brief calculates the value of a given spatial function for a given gnode

  <pre>                                                              mn 06/02
  This routine calculates the value of a given spatial function for a given
  gnode

  </pre>
  \param *gnode  GNODE  (i)   the gnode
  \param  index  INT    (i)   the index of dirich condition

  \warning There is nothing special to this routine
  \return void
  \sa calling: ---; called by: inherit_design_dis_dirichlet()

 *----------------------------------------------------------------------*/
void cal_dirich_fac(
    GNODE              *gnode,
    INT                 index
    )
{

  DOUBLE  xi;
  DOUBLE  length,length_1,length_2;

  INT     funct_num;
  DOUBLE  xp[3],x1[3],x2[3];
  DOUBLE  fac;

  DOUBLE    a,d;
  DOUBLE    h,um;

  FUNCT_LINE_LIN     *f_line_lin;
  FUNCT_RADIUS_LIN   *f_radius_lin;
  FUNCT_LINE_QUAD    *f_line_quad;
  FUNCT_RADIUS_QUAD  *f_radius_quad;
  FUNCT_CYL          *f_cyl;

#ifdef DEBUG
  dstrc_enter("cal_dirich_fac");
#endif


  funct_num = gnode->dirich->funct.a.iv[index] - 1;
  xp[0]     = gnode->node->x[0];
  xp[1]     = gnode->node->x[1];
  xp[2]     = gnode->node->x[2];

  /* switch to the correct funtion */
  switch (funct[funct_num].functtyp)
  {

    case funct_line_lin:  /* linear function */
      f_line_lin = funct[funct_num].typ.funct_line_lin;

      length = f_line_lin->length;

      /* x1: vector along the line */
      x1[0] = f_line_lin->x2[0] - f_line_lin->x1[0];
      x1[1] = f_line_lin->x2[1] - f_line_lin->x1[1];
      x1[2] = f_line_lin->x2[2] - f_line_lin->x1[2];

      /* x2: vector from the beginning of the line the the point */
      x2[0] = xp[0] - f_line_lin->x1[0];
      x2[1] = xp[1] - f_line_lin->x1[1];
      x2[2] = xp[2] - f_line_lin->x1[2];

      /* length_1 = projection of x2 onto x1 */
      length_1 = ( x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] ) / length;
      /* length_2 = length of the vector x2 */
      length_2 = sqrt( x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2] );

      /* check for a point not on the line */
      if ( FABS(length_1 - length_2) > 10e-6 )
        dswarning(1,6);

      /* calculate xi and check for a point outside the range of the funct */
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        dswarning(1,5);

      /* calculate function value at point p */
      fac = f_line_lin->b + xi * f_line_lin->m;
      break;




    case funct_radius_lin:  /* linear function */
      f_radius_lin = funct[funct_num].typ.funct_radius_lin;

      length = f_radius_lin->length;

      /* x2: vector from the beginning of the line the the point */
      x1[0] = xp[0] - f_radius_lin->x1[0];
      x1[1] = xp[1] - f_radius_lin->x1[1];
      x1[2] = xp[2] - f_radius_lin->x1[2];

      /* length_1 = length of the vector x1 */
      length_1 = sqrt( x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2] );

      /* calculate xi and check for a point outside the range of the funct */
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        dswarning(1,5);

      /* calculate function value at point p */
      fac = f_radius_lin->b + xi * f_radius_lin->m;
      break;





    case funct_line_quad: /*quadratic parabola */
      f_line_quad = funct[funct_num].typ.funct_line_quad;

      length = f_line_quad->length;

      /* x1: vector along the line */
      x1[0] = f_line_quad->x2[0] - f_line_quad->x1[0];
      x1[1] = f_line_quad->x2[1] - f_line_quad->x1[1];
      x1[2] = f_line_quad->x2[2] - f_line_quad->x1[2];

      /* x2: vector from the beginning of the line the the point */
      x2[0] = xp[0] - f_line_quad->x1[0];
      x2[1] = xp[1] - f_line_quad->x1[1];
      x2[2] = xp[2] - f_line_quad->x1[2];

      /* length_1 = projection of x2 onto x1 */
      length_1 = ( x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] ) / length;
      /* length_2 = length of the vector x2 */
      length_2 = sqrt( x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2] );

      /* check for a point not on the line */
      if ( FABS(length_1 - length_2) > 10e-6 )
        dswarning(1,6);

      /* calculate xi and check for a point outside the range of the funct */
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        dswarning(1,5);

      /* calculate function value at point p */
      fac = 1.0 - 4 * (xi - 1.0/2.0)*(xi - 1.0/2.0);
      break;






    case funct_radius_quad: /*quadratic parabola */
      f_radius_quad = funct[funct_num].typ.funct_radius_quad;

      length = f_radius_quad->length;

      /* x1: vector from the beginning of the line the the point */
      x1[0] = xp[0] - f_radius_quad->x1[0];
      x1[1] = xp[1] - f_radius_quad->x1[1];
      x1[2] = xp[2] - f_radius_quad->x1[2];

      /* length_1 = length of the vector x1 */
      length_1 = sqrt( x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2] );

      /* calculate xi and check for a point outside the range of the funct */
      xi =  length_1/length;
      if ( xi > 1.0)
        dswarning(1,5);

      /* calculate function value at point p */
      fac = 1.0 - xi * xi ;
      break;






    case funct_bel:  /* spatial function for beltrami flow */
      /* set some constants */
      a    = PI/4.0;
      d    = PI/2.0;
      /* calculate values */
      switch (index)
      {
        case 0:
          fac  = -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
              exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) );
          break;
        case 1:
          fac  = -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
              exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) );
          break;
        case 2:
          fac  = -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
              exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) );
          break;
        case 3:
          fac  = -a*a/2 * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
              + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
              + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
              + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1])));
          break;
        default:
          fac = 1.0;
          break;
      }
      break;


    case funct_kim:  /* spatial function for kim-moin flow */
      /* set some constants */
      a    = 2.0;
      /* calculate values */
      switch (index)
      {
        case 0:
          fac  = - cos(a*PI*xp[0]) * sin(a*PI*xp[1]);
          break;
        case 1:
          fac  = + sin(a*PI*xp[0]) * cos(a*PI*xp[1]);
          break;
        case 2:
          fac  = -1.0/4.0 * ( cos(2.0*a*PI*xp[0]) + cos(2.0*a*PI*xp[1]) );
          break;
        default:
          fac = 1.0;
          break;
      }
      break;


    case funct_cyl:
      f_cyl = funct[funct_num].typ.funct_cyl;

      /* set some constants */
      h    = 0.41;
      um = f_cyl->um;

      /* calculate values */
      fac = 16*um*xp[1]*xp[2]*(h-xp[1])*(h-xp[2]) / (h*h*h*h);
      break;


    default:  /* default: no function */
      fac = 1.0;
      break;
  }


  /* write the factor to the gnode */
  gnode->d_funct[index] = fac;


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of cal_dirich_fac */


/*! @} (documentation module close)*/




