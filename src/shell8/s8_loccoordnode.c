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
 | local coordinates of nodal points                      m.gee 6/01    |
 *----------------------------------------------------------------------*/
DOUBLE s8_local_coord_node(INT node, INT flag, enum _DIS_TYP typ)
{
DOUBLE     coord;
const DOUBLE node489[9][2] = {{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{1.0,-1.0},
                              {0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,0.0}};
#ifdef DEBUG
dstrc_enter("s8_local_coord_node");
#endif
/*----------------------------------------------------------------------*/
switch(typ)
{
case quad4:
   coord = node489[node][flag];
break;
case quad8:
   coord = node489[node][flag];
break;
case quad9:
   coord = node489[node][flag];
break;
case tri3:
break;
case tri6:
break;
default:
   dserror("unknown number of nodes");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return(coord);
} /* end of s8_local_coord_node */
#endif
