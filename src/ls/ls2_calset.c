#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calset( 
  ELEMENT         *ele,     
  DOUBLE         **xyze,
  DOUBLE          *lset00,    
  DOUBLE          *lset01
  )
{
  INT      i;
  NODE    *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls2_calset");
#endif
/*----------------------------------------------------------------------*/
  
  /* set element coordinates */
  for(i=0; i<ele->numnp; i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }
  
  /*--------------------------------------------------------------------*
   | position of the different solutions:                               |
   | node->sol_incement: solution history used for calculations         |
   |       sol_increment[0][0]: solution at (n)                         |
   |	 sol_increment[0][1]: solution at (n+1)                         |
   |	 sol_increment[1][2]: solution at (n+1)^0                       |
   *--------------------------------------------------------------------*/

  for(i=0; i<ele->numnp; i++)
  {
    actnode=ele->node[i];
    /* nodal values of levelset function at time step (n  ) */      
    lset00[i]=actnode->sol_increment.a.da[0][0];
    /* nodal values of levelset function at time step (n+1) */      	
    lset01[i]=actnode->sol_increment.a.da[1][0];
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return; 
} /* end of ls2_calset */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calset1( 
  ELEMENT         *ele,
  INT              pos,
  DOUBLE          *lset
  )
{
  INT      i;
  NODE    *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls2_calset1");
#endif
/*----------------------------------------------------------------------*/
  
  /*--------------------------------------------------------------------*
   | position of the different solutions:                               |
   | node->sol_incement: solution history used for calculations         |
   |     sol_increment[0][0]: solution at (n)                           |
   |	 sol_increment[0][1]: solution at (n+1)                         |
   |	 sol_increment[1][2]: solution at (n+1)^0                       |
   *--------------------------------------------------------------------*/
  for(i=0; i<ele->numnp; i++)
  {
    actnode=ele->node[i];
    /* nodal values of levelset function at time step (n+1) */      	
    lset[i]=actnode->sol_increment.a.da[pos][0];
  }
  
/*---------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return; 
} /* end of ls2_calset1 */
#endif
