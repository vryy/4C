/*!----------------------------------------------------------------------
\file
\brief  contains the routine 'wge_disd' which calculates displacement 
        derivatives for gradient enhanced wall element

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*! 
\addtogroup WALLGE 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  contains the routine 'wge_disd' which calculates displacement 
        derivatives for gradient enhanced wall element

*----------------------------------------------------------------------*/
void wge_strain(ELEMENT   *ele,     /* actual element                  */
                DOUBLE   **bopd,    /* B-operator for displacements    */
                DOUBLE    *functe,  /* Ansatz-funct. for equiv. strain */
                DOUBLE   **bope,    /* B-operator for equiv. strain    */
                DOUBLE    *strain,  /* stress vector                   */
                DOUBLE    *eps_vnl, /* nonlocal equivalent strain      */
                DOUBLE    *grad_eps_vnl) /* grad nonl. equiv. strain   */
{
#ifdef D_WALLGE

INT i,dnode, node_start,enode;
DOUBLE disd[4];


#ifdef DEBUG 
dstrc_enter("wge_strain");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<4; i++) disd[i]   = 0.0;
for (i=0; i<4; i++) strain[i] = 0.0;
*eps_vnl = 0.0;
grad_eps_vnl[0]=0.0;
grad_eps_vnl[1]=0.0;
/*----------------------------------------- displacement derivatives ---*/
for (dnode=0; dnode<ele->numnp; dnode++)
{
  node_start = dnode*2;
  disd[0] += bopd[0][node_start+0]*ele->node[dnode]->sol.a.da[0][0];
  disd[1] += bopd[1][node_start+1]*ele->node[dnode]->sol.a.da[0][1];
  disd[2] += bopd[2][node_start+0]*ele->node[dnode]->sol.a.da[0][0];
  disd[3] += bopd[2][node_start+1]*ele->node[dnode]->sol.a.da[0][1];
} /* end of loop over displacement-nodes */
/*------------------------------------------ strains exept strain_33 ---*/
strain[0]  = disd[0];
strain[1]  = disd[1];
strain[2]  = disd[2] + disd[3];
/*--------------------------------------- nonlocal equivalent strain ---*/
for (enode=0; enode<4; enode++)
{
  *eps_vnl += functe[enode]*ele->node[enode]->sol.a.da[0][2];
} /* end of loop over equivalent-strain-nodes */
/*--------------------------- gradient of nonlocal equivalent strain ---*/
for (enode=0; enode<4; enode++)
{
  grad_eps_vnl[0] += bope[0][enode]*ele->node[enode]->sol.a.da[0][2];
  grad_eps_vnl[1] += bope[1][enode]*ele->node[enode]->sol.a.da[0][2];
} /* end of loop over equivalent-strain-nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_strain */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
