/*!----------------------------------------------------------------------
\file
\brief contains routines to serve elemental geometry calculations

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"



/*!----------------------------------------------------------------------
\brief  area of linear 2D elements

<pre>                                                              ck 01/03
This routine evaluates the area of 2D linear elements (quad4 and tri3)
The area of higher order elements is calculated ignoring the mid-side nodes, 
i.e. assuming stright sides. 

</pre>
\param  *ele        ELEMENT   (i)   element pointer
\param   xyz[2][9]  DOUBLE    (i)   elemental coordinates 1st index: x or y
                                                          2dn index: node

\warning This routine gives an exact result for quad4 elements only. For 
         higher order quads it offers an estimation.
\return DOUBLE element area
\sa calling:
             called by: ale2_static_ke_prestress(),
	                ale2_static_ke_step2()

*----------------------------------------------------------------------*/
DOUBLE area_lin_2d(ELEMENT *ele, DOUBLE xyz[2][MAXNOD])
{
DOUBLE a, b, c;  /* geometrical values */
DOUBLE el_area;  /* element area */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("area_lin_2d");
#endif
switch(ele->distyp)
{
case tri3:
   a = (xyz[0][0]-xyz[0][1])*(xyz[0][0]-xyz[0][1])
      +(xyz[1][0]-xyz[1][1])*(xyz[1][0]-xyz[1][1]); /* line 0-1 squared */
   b = (xyz[0][1]-xyz[0][2])*(xyz[0][1]-xyz[0][2])
      +(xyz[1][1]-xyz[1][2])*(xyz[1][1]-xyz[1][2]); /* line 1-2 squared */
   c = (xyz[0][2]-xyz[0][0])*(xyz[0][2]-xyz[0][0])
      +(xyz[1][2]-xyz[1][0])*(xyz[1][2]-xyz[1][0]); /* diag 2-0 squared */
   el_area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
break;
case quad4:
case quad9:
   a = (xyz[0][0]-xyz[0][1])*(xyz[0][0]-xyz[0][1])
      +(xyz[1][0]-xyz[1][1])*(xyz[1][0]-xyz[1][1]); /* line 0-1 squared */
   b = (xyz[0][1]-xyz[0][2])*(xyz[0][1]-xyz[0][2])
      +(xyz[1][1]-xyz[1][2])*(xyz[1][1]-xyz[1][2]); /* line 1-2 squared */
   c = (xyz[0][2]-xyz[0][0])*(xyz[0][2]-xyz[0][0])
      +(xyz[1][2]-xyz[1][0])*(xyz[1][2]-xyz[1][0]); /* diag 2-0 squared */
   el_area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
   a = (xyz[0][2]-xyz[0][3])*(xyz[0][2]-xyz[0][3])
      +(xyz[1][2]-xyz[1][3])*(xyz[1][2]-xyz[1][3]); /* line 2-3 squared */
   b = (xyz[0][3]-xyz[0][0])*(xyz[0][3]-xyz[0][0])
      +(xyz[1][3]-xyz[1][0])*(xyz[1][3]-xyz[1][0]); /* line 3-0 squared */
   /*----------------------------------------- evaluate element area ---*/
   el_area += 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
break;
default:
 dserror("evaluation of area of actual element not implemented!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return el_area;
}
