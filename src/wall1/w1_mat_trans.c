#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*-----------------------------------------------------------------------|
|      topic: blowing up plane stress/strain conditions           sh 7/02|
|             to 3D --> 3D-Material Law                                  |
|-----------------------------------------------------------------------*/
void w1mat_trans_up (double **d,/*current material matrix d14-d44       */
                     ELEMENT   *ele,                                        
                     WALL_TYPE wtype,
                     double  *stress,  /*actuel stress [4]              */
                     double  *strain,  /*actual strain [4]              */
                     double  *qn,
                     double  **bop)
{
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_enter("w1mat_trans_up");
    #endif
/*---------------------------- condensed stress- or backstressvector ---*/
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1mat_trans_up */



/*-----------------------------------------------------------------------|
|      topic: kondense 3D conditions                              sh 7/02|
|             to plane stress/strain conditions                          |
|-----------------------------------------------------------------------*/
void w1mat_trans_down (double **d,/*current material matrix d14-d44     */
                     ELEMENT   *ele,                                        
                     WALL_TYPE wtype,
                     double  *stress,  /*actuel stress [4]              */
                     double  *strain,  /*actual strain [4]              */
                     double  *qn,
                     double  **bop)
{
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_enter("w1mat_trans_down");
    #endif
/*---------------------------- condensed stress- or backstressvector ---*/
/*----------------------------------------------------------------------*/
    #ifdef DEBUG 
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1mat_trans_down */
