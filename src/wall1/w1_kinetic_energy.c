/*!-----------------------------------------------------------------------------
\file
\brief contains the routine 'w1_kinetic_energy' which calculates the kinetic energy
linear and angular momentum of a wall element

*-----------------------------------------------------------------------------*/
#ifdef GEMM
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------------*
 | Calculation of kinetic energy and lin. and angular momentum for an element /                              |
 *---------------------------------------------------------------------------*/
void w1_kinetic_energy(ELEMENT *ele, DOUBLE **mass)

{
const INT           numdf  = 2;     /* number dof per node          */
const INT           numeps = 4;     /* number of strain components  */
INT                 iel;            /* numnp to this element        */
INT                 nd;             /* dof of this element          */
ARRAY      velocity_components_a;
DOUBLE    *velocity_components;
ARRAY      linear_momentum_a;
DOUBLE    *linear_momentum;
ARRAY      angular_momentum_a;
DOUBLE    *angular_momentum;
ARRAY      auxilliary_a;
DOUBLE   **auxilliary;
ARRAY      xcure_a;
DOUBLE   **xcure;
DOUBLE    ang_momentum, kinetic_energy;
INT       i,j,k;


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_kinetic_energy");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------Initialization phase*/
velocity_components = amdef("velocity", &velocity_components_a,(numdf*MAXNOD_WALL1),1,"DV");
linear_momentum     = amdef("linear_mom", &linear_momentum_a,(numdf*MAXNOD_WALL1),1,"DV");
angular_momentum    = amdef("ang_mom", &angular_momentum_a,MAXNOD_WALL1,1,"DV");
xcure               = amdef("xcure",  &xcure_a,2,MAXNOD_WALL1,"DA");
auxilliary          = amdef("auxil",  &auxilliary_a,MAXNOD_WALL1,(numdf*MAXNOD_WALL1),"DA");

ang_momentum = 0.0;
kinetic_energy = 0.0;
iel     = ele->numnp;
nd      = numdf * iel;

amzero(&velocity_components_a);
amzero(&linear_momentum_a);
amzero(&angular_momentum_a);
amzero(&xcure_a);
amzero(&auxilliary_a);
/*---------------------------Current position of the nodes of an element*/
for(k=0; k<iel; k++) {
 xcure[0][k]= ele->node[k]->x[0] + ele->node[k]->sol.a.da[0][0]; 
 xcure[1][k]= ele->node[k]->x[1] + ele->node[k]->sol.a.da[0][1];
 }


for(i=0; i<iel; i++){
  angular_momentum[i] = 0.0;
  for(j=0; j<nd; j++) auxilliary[i][j] = 0.0;
}

/*---Auxilliary matrix which is used in calculation of r cross product mv*/
k=0;
for(i=0; i<iel; i++){
  auxilliary[i][k]   = -xcure[1][i];
  auxilliary[i][k+1] =  xcure[0][i];
  k = k+2;
  }

for(i=0; i<nd; i++){
  velocity_components[i] = 0.0;
  linear_momentum[i] = 0.0;
  }

/*-------------------------Velocity components of the nodes of an element*/
i=0;  
for(k=0; k<iel; k++){
  for(j=0; j<2; j++){
  velocity_components[i] = ele->node[k]->sol.a.da[1][j];
  i++;
  }
}  

/*-----------------------------Calculation of linear momentum  : L = mv */
for(i=0; i<nd; i++)
  for(j=0;j<nd; j++)
  linear_momentum[i] += mass[i][j] * velocity_components[j]; 

 /*-------------------Calculation of angular momentum: J = Auxilliary*L */ 
for(i=0; i<iel; i++)
  for(j=0; j<nd; j++) 
  angular_momentum[i] += auxilliary[i][j]*linear_momentum[j];


/*-------------------------------------Linear momentum x-component : Lx */
for(i=0;i<nd; i=i+2)
 ele->e.w1->linmom[0] += linear_momentum[i]; 
 
/*-------------------------------------Linear momentum y-component : Ly */
for(i=1;i<nd; i=i+2)
 ele->e.w1->linmom[1] += linear_momentum[i]; 

/*-----------------------------------------------------Angular momentum */
for(i=0; i<iel; i++)
ang_momentum += angular_momentum[i];

ele->e.w1->angular_momentum = ang_momentum; 


/*-----------------------------------------Calculation of kinetic energy*/
for(i=0; i<nd; i++)
  for(j=0; j<nd; j++)
  kinetic_energy +=  0.5 * velocity_components[i] * mass[i][j] * velocity_components[j];

ele->e.w1->kinetic_energy = kinetic_energy;


/*----------------------------------------------------Cleaning-up phase */
amdel(&velocity_components_a);
amdel(&linear_momentum_a);
amdel(&angular_momentum_a);
amdel(&xcure_a);
amdel(&auxilliary_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
}/* end of w1_kinetic_energy */
#endif
/*! @} (documentation module close)*/
