/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'if_dirichnode' which eleminates 2 of the
       nodes by dirichlet conditions if the element is quadratic

*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief eleminates 2 of the nodes by dirichlet conditions if the 
       element is quadratic

<pre>                                                              mn 05/03
This routine eleminates 2 of the nodes by dirichlet conditions if the 
       element is quadratic

</pre>
\param **s       DOUBLE    (o)  blablabla 
\param   dl      DOUBLE    (i)  blablabal

\warning There is nothing special to this routine
\return void                                               
\sa calling:  ---; 
    caled by: assign_dof();

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


void if_dirichnode(ELEMENT       *actele)
{
INT      cnode,i,mynode;              /* some loopers     */
DOUBLE   xrefe[4],yrefe[4];  /* reference coordinates of corner nodes */
DOUBLE   L_one, L_two;       /* lengh of element edges */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_dirichnode");
#endif
/*----------- check orientation of element (which is my xi direction)---*/
for (cnode=0; cnode<4; cnode++)
{/* coordinates of corner nodes */
  xrefe[cnode] = actele->node[cnode]->x[0];          
  yrefe[cnode] = actele->node[cnode]->x[1];                
}
L_one = sqrt( (xrefe[1] - xrefe[0]) * (xrefe[1] - xrefe[0])
      +       (yrefe[1] - yrefe[0]) * (yrefe[1] - yrefe[0]));
L_two = sqrt( (xrefe[2] - xrefe[1]) * (xrefe[2] - xrefe[1])
      +       (yrefe[2] - yrefe[1]) * (yrefe[2] - yrefe[1]));
/*----------------------------------------------------------------------*/
if (L_one>L_two)
 {
  actele->node[5]->gnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
  actele->node[7]->gnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
  amdef("onoff",&(actele->node[5]->gnode->dirich->dirich_onoff),6,1,"IV");
  amdef("val",&(actele->node[5]->gnode->dirich->dirich_val),6,1,"DV");
  amdef("curve",&(actele->node[5]->gnode->dirich->curve),6,1,"IV");
  amdef("onoff",&(actele->node[7]->gnode->dirich->dirich_onoff),6,1,"IV");
  amdef("val",&(actele->node[7]->gnode->dirich->dirich_val),6,1,"DV");
  amdef("curve",&(actele->node[7]->gnode->dirich->curve),6,1,"IV");

  actele->node[5]->gnode->dirich->dirich_onoff.a.iv[0]=1; 
  actele->node[5]->gnode->dirich->dirich_onoff.a.iv[1]=1; 
  for (i=2; i<6; i++)
  actele->node[5]->gnode->dirich->dirich_onoff.a.iv[i]=0;
   
  for (i=0; i<6; i++)
  {
    actele->node[5]->gnode->dirich->dirich_val.a.dv[i]=0.0; 
    actele->node[5]->gnode->dirich->curve.a.iv[i]=0; 
  }

  actele->node[7]->gnode->dirich->dirich_onoff.a.iv[0]=1; 
  actele->node[7]->gnode->dirich->dirich_onoff.a.iv[1]=1; 
  for (i=2; i<6; i++)
  actele->node[7]->gnode->dirich->dirich_onoff.a.iv[i]=0; 
  
  for (i=0; i<6; i++)
  {
    actele->node[7]->gnode->dirich->dirich_val.a.dv[i]=0.0;
    actele->node[7]->gnode->dirich->curve.a.iv[i]=0; 
  } 
 }
else if (L_two>L_one)
{
  actele->node[4]->gnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
  actele->node[6]->gnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
  amdef("onoff",&(actele->node[4]->gnode->dirich->dirich_onoff),6,1,"IV");
  amdef("val",&(actele->node[4]->gnode->dirich->dirich_val),6,1,"DV");
  amdef("curve",&(actele->node[4]->gnode->dirich->curve),6,1,"IV");
  amdef("onoff",&(actele->node[6]->gnode->dirich->dirich_onoff),6,1,"IV");
  amdef("val",&(actele->node[6]->gnode->dirich->dirich_val),6,1,"DV");
  amdef("curve",&(actele->node[6]->gnode->dirich->curve),6,1,"IV");

  actele->node[4]->gnode->dirich->dirich_onoff.a.iv[0]=1; 
  actele->node[4]->gnode->dirich->dirich_onoff.a.iv[1]=1;
  for (i=2; i<6; i++)
  actele->node[4]->gnode->dirich->dirich_onoff.a.iv[i]=0; 

  for (i=0; i<6; i++)
  {
    actele->node[4]->gnode->dirich->dirich_val.a.dv[i]=0.0; 
    actele->node[4]->gnode->dirich->curve.a.iv[i]=0; 
  }

  actele->node[6]->gnode->dirich->dirich_onoff.a.iv[0]=1; 
  actele->node[6]->gnode->dirich->dirich_onoff.a.iv[1]=1; 
  for (i=2; i<6; i++)
  actele->node[6]->gnode->dirich->dirich_onoff.a.iv[i]=0; 
  
  for (i=0; i<6; i++)
  {
    actele->node[6]->gnode->dirich->dirich_val.a.dv[i]=0.0; 
    actele->node[6]->gnode->dirich->curve.a.iv[i]=0; 
  }
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of if_dirichnode */



/*!----------------------------------------------------------------------
\brief routine 'if_permstiff' -> resort of stiffness stiff
       into the element stiffness matrix "estif"

*----------------------------------------------------------------------*/
void if_permstiff(DOUBLE **estif,
                  DOUBLE **Kdd,
                  INT      iele,
                  INT      ield)   /* "mixed" element stiffness   */      
{

INT            i,j;
INT            nodei,nodej;
INT            nodestarti,nodestartj;
INT            dofi,dofj;
INT            numdf;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_permstiff");
#endif
/*-------------------- (upper left part of estif) node 1-4 * node1-4 ---*/
numdf = 2* ield + iele;
for (nodei=0; nodei<iele; nodei++)
{
  for (nodej=0; nodej<iele; nodej++)
  {
    dofi       = 2 * nodei;
    dofj       = 2 * nodej;
    nodestarti = 3 * nodei;
    nodestartj = 3 * nodej;
    /*-----------------------------------------------------*/
    estif[nodestarti][nodestartj]     = Kdd[dofi][dofj];
    estif[nodestarti+1][nodestartj]   = Kdd[dofi+1][dofj];
    estif[nodestarti][nodestartj+1]   = Kdd[dofi][dofj+1];
    estif[nodestarti+1][nodestartj+1] = Kdd[dofi+1][dofj+1];
    /*-----------------------------------------------------*/
    estif[nodestarti][nodestartj+2]   = 0.0;
    estif[nodestarti+1][nodestartj+2] = 0.0;
    /*-----------------------------------------------------*/
    estif[nodestarti+2][nodestartj]   = 0.0;
    estif[nodestarti+2][nodestartj+1] = 0.0;
    /*-----------------------------------------------------*/
    estif[nodestarti+2][nodestartj+2] = 0.0;
  }
}
/*-------------------------------------------- (rest part of estif) ---*/
if(ield>iele)
{
/*----------------- (upper right part of estif) node1-4 * node4-8/9 ---*/
  for (nodei=0; nodei<iele; nodei++)
  {
    for (j=12; j<numdf; j++)
    {
      dofi       = 2 * nodei;
      dofj       = j-4;
      nodestarti = 3 * nodei;
      /*-----------------------------------------------------*/
      estif[nodestarti][j]     = Kdd[dofi][dofj];
      estif[nodestarti+1][j]   = Kdd[dofi+1][dofj];
      /*-----------------------------------------------------*/
      estif[nodestarti+2][j]   = 0.0;
    }
  }
/*------------------- (down left part of estif) node4-8/9 * node1-4 ---*/
  for (i=12; i<numdf; i++)
  {
    for (nodej=0; nodej<iele; nodej++)
    {
      dofi       = i-4;
      dofj       = 2 * nodej;
      nodestartj = 3 * nodej;
      /*-----------------------------------------------------*/
      estif[i][nodestartj]     = Kdd[dofi][dofj];
      estif[i][nodestartj+1]   = Kdd[dofi][dofj+1];
      /*-----------------------------------------------------*/
      estif[i][nodestartj+2]   = 0.0;
    }
  }
/*---------------- (down right part of estif) node4-8/9 * node4-8/9 ---*/
  for (i=12; i<numdf; i++)
  {
    for (j=12; j<numdf; j++)
    {
      dofi       = i-4;
      dofj       = j-4;
      /*-----------------------------------------------------*/
      estif[i][j] = Kdd[dofi][dofj];
    }
  }
} /*-----endif-*/

#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of if_permstiff */


/*!----------------------------------------------------------------------
\brief routine 'if_permforce' -> resort of internal force parts fint
       into the element internal force "force"

*----------------------------------------------------------------------*/
void if_permforce(DOUBLE    *force,   /* "mixed" element int. force   */
                  DOUBLE    *fintd,    /*  int. force           */
                   INT       iele,    /* num.of equiv.strain nodes   */
                   INT       ield)    /* num of displacement nodes  */      
{

INT            i;
INT            nodei;
INT            nodestarti;
INT            dofi;
INT            numdf;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_permforce");
#endif
/*----------------------------------- (upper part of force) node 1-4 ---*/
numdf = 2* ield + iele;
for (nodei=0; nodei<iele; nodei++)
{
    dofi       = 2 * nodei;
    nodestarti = 3 * nodei;
    /*---------------------------------------*/
    force[nodestarti]   = fintd[dofi];
    force[nodestarti+1] = fintd[dofi+1];
    /*---------------------------------------*/
    force[nodestarti+2] = 0.0;
}
/*--------------------------------------------- (rest part of force) ---*/
if(ield>iele)
{
  for (i=12; i<numdf; i++)
  {
    dofi = i - 4;
    force[i] = fintd[dofi];
  }
} /*-----endif-*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of if_permforce */


/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
