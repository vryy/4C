#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | local coordinates of nodal points                      m.gee 6/01    |
 *----------------------------------------------------------------------*/
double s8_local_coord_node(int node, int flag, ELEMENT_TYP typ)
{
double     coord;
const double node489[9][2] = {{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{1.0,-1.0},
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
