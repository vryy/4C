/*!----------------------------------------------------------------------
\file
\brief contains the function 
 - s9_local_coord_node: which returns the local coordinates of the nodal 
                        points within an element

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief returns the coordinates of nodal points                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This function returns the local coordinates of nodal points within an
element
</pre>
\param  INT node          (i) local node number within element
\param  INT flag          (i) direction (0:r-direction; 1:s-direction)
\param  enum _DIS_TYP typ (i) type (quad4/8/9)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9a3(); s9a3ref_extern()  [s9_a3.c]

*----------------------------------------------------------------------*/
DOUBLE s9_local_coord_node(INT node, INT flag, enum _DIS_TYP typ)
{
DOUBLE     coord;
const DOUBLE node489[9][2] = {{1.0,1.0},{-1.0,1.0},{-1.0,-1.0},{1.0,-1.0},
                              {0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,0.0}}; 
#ifdef DEBUG 
dstrc_enter("s9_local_coord_node");
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
} /* end of s9_local_coord_node */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
