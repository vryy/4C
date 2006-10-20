/*======================================================================*/
/*!
\file
\brief Routines cope with THERM2 element configurations

\author bborn
\date 03/06
*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*----------------------------------------------------------------------*/
/*!
\brief Locally globals

\author bborn
\date 03/06
*/
static INT iegq[4][4][2];
static INT iegt[4][4][2];

/*======================================================================*/
/*!
\brief Initialise local globals

\author bborn
\date 03/06
*/
void th2_cfg_init()
{
  INT i;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_cfg_init");
#endif

  /*--------------------------------------------------------------------*/
   /*-------------------------------------------- egde nodes for quad4 */
   iegq[0][0][0] = 0;
   iegq[1][0][0] = 1;
   iegq[0][1][0] = 1;
   iegq[1][1][0] = 2;
   iegq[0][2][0] = 2;
   iegq[1][2][0] = 3;
   iegq[0][3][0] = 3;
   iegq[1][3][0] = 0;
   /*----------------------------------- egde nodes for quad8 and quad9 */
   iegq[0][0][1] = 0;
   iegq[1][0][1] = 4;
   iegq[2][0][1] = 1;
   iegq[0][1][1] = 1;
   iegq[1][1][1] = 5;
   iegq[2][1][1] = 2;
   iegq[0][2][1] = 2;
   iegq[1][2][1] = 6;
   iegq[2][2][1] = 3;
   iegq[0][3][1] = 3;
   iegq[1][3][1] = 7;
   iegq[2][3][1] = 0;
   /*---------------------------------------------- egde nodes for tri3 */
   iegt[0][0][0] = 0;
   iegt[1][0][0] = 1;
   iegt[0][1][0] = 1;
   iegt[1][1][0] = 2;
   iegt[0][2][0] = 2;
   iegt[1][2][0] = 0;
   /*---------------------------------------------- egde nodes for tri6 */
   iegt[0][0][1] = 0;
   iegt[1][0][1] = 3;
   iegt[2][0][1] = 1;
   iegt[0][1][1] = 1;
   iegt[1][1][1] = 4;
   iegt[2][1][1] = 2;
   iegt[0][2][1] = 2;
   iegt[1][2][1] = 5;
   iegt[2][2][1] = 0;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_cfg_init */

/*======================================================================*/
/*!
\brief 

What the heck is this? Must be some sort of numbering of the nodes or
DOFs of the element...

\param  *iegnod     INT      (o)   element nodes (indices) on actual line
\param  *ele        ELEMENT  (i)   actual element
\param  line        INT      (i)   actual element line (edge)
\return void

\author bborn
\date 03/06
*/
void th2_cfg_iedg(INT *iegnod, 
                  ELEMENT *ele, 
                  INT line)
{
INT i;

#ifdef DEBUG
dstrc_enter("th2_cfg_iedg");
#endif

/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/

   switch(ele->distyp)
   {
   case quad4:
      for(i=0;i<2;i++) iegnod[i] = iegq[i][line][0];
   break;
   case quad8: case quad9:
      for(i=0;i<3;i++) iegnod[i] = iegq[i][line][1];
   break;
   case tri3:
      for(i=0;i<2;i++) iegnod[i] = iegt[i][line][0];
   break;
   case tri6:
      for(i=0;i<3;i++) iegnod[i] = iegt[i][line][1];
/*      dserror("iegnode for tri6 not tested yet\n"); */
   break;
   default:
      dserror("distyp unknown\n");
   }  /* end of switch(ele->distyp) */

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of th2_cfg_iedg */


/*======================================================================*/
/*!
\brief Natural ccordinates of element nodes

Obtain the rs-coordinates (ie natural coords, or the coordinates in 
parameter space)  of the element nodes (ie quads with 4,8,9, or tris
with 3,6)

\author bborn
\date 03/06
*/
void th2_cfg_noders(ELEMENT *ele,
                    INT inode,
                    DOUBLE *rs)
{
  static DOUBLE qrs489[9][NDIM_THERM2] = { {1.0,1.0}, {-1.0,1.0}, 
                                           {-1.0,-1.0}, {1.0,-1.0},
                                           {0.0,1.0}, {-1.0,0.0},
                                           {0.0,-1.0}, {1.0,0.0},
                                           {0.0,0.0} };
  static DOUBLE trs36[6][NDIM_THERM2] = { {0.0,0.0}, {1.0,0.0}, {0.0,1.0},
                                          {0.5,0.0}, {0.5,0.5}, {0.0,0.5} };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_cfg_noders");
#endif

  /*--------------------------------------------------------------------*/
  switch (ele->distyp)
  {
    case quad4: case quad8: case quad9:
      rs[0] = qrs489[inode][0];
      rs[1] = qrs489[inode][1];
      break;
    case tri3: case tri6:
      rs[0] = trs36[inode][0];
      rs[1] = trs36[inode][1];
      break;
    default:
      dserror("Unknown discretisation type!");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_cfg_noders */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM2 */
