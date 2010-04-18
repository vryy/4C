/*!----------------------------------------------------------------------
\file
\brief contains the routine

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
#include "s2ml_prototypes.h"
#include "s2ml.h"

/*!
\addtogroup MLSTRUCT
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine for calculation of equivalent strain in the macroelement
-> decision: only large scale or multiscale calculation in this macroelement

\param  *actmaele         ELEMENT     (I)   actual Macro-element
\param   wtype            WALL_TYPE   (I)   plane stress/strain...
\param **bopma            DOUBLE      (I)   Macro-B-operator
\param   nue              DOUBLE      (I)   poisson-ratio
\param  *eps_equiv        DOUBLE      (O)   equivalent strain

\return void

*----------------------------------------------------------------------*/
void s2ml_aequistrain(ELEMENT     *actmaele,  /* actual Macro-element   */
                      WALL_TYPE    wtype,    /* plane stress/strain... */
                      DOUBLE     **bopma,    /* Macro-B-operator       */
                      DOUBLE       nue,      /* poisson-ratio          */
                      DOUBLE      *eps_equiv)/* equivalent strain      */
{
INT i,j,inode, node_start;
DOUBLE I1,J2,scalar;
DOUBLE disd[5];
DOUBLE strain[4];
DOUBLE epsilon[3][3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_aequistrain");
#endif
/*------------------------------- calculate displacement derivatives ---*/
for (i=0; i<5; i++) disd[i] = 0.0;
for (inode=0; inode<actmaele->numnp; inode++)
{
  node_start = inode*2;
  disd[0] += bopma[0][node_start+0]*actmaele->node[inode]->sol.a.da[0][0];
  disd[1] += bopma[1][node_start+1]*actmaele->node[inode]->sol.a.da[0][1];
  disd[2] += bopma[2][node_start+0]*actmaele->node[inode]->sol.a.da[0][0];
  disd[3] += bopma[2][node_start+1]*actmaele->node[inode]->sol.a.da[0][1];
} /* end of loop over nodes */

/*--------------------------------------------------- strain vector ---*/
strain[0] = disd[0];
strain[1] = disd[1];
strain[2] = disd[2] + disd[3];

if(wtype==plane_strain)
{
  strain[3] = 0.;
}
else if (wtype == plane_stress)
{
 strain[3] = - (nue*(strain[0]+strain[1]))/(1.0 - nue);
}
/*--------------------------- calculate strainvector to straintensor ---*/
w1_4to9(strain,epsilon);
/*--------------------------- calculate equivalent strain meassure ---*/
I1 = epsilon[0][0]+epsilon[1][1]+epsilon[2][2];
scalar= 0.0;
for (i=0; i<3; i++)
{
 for (j=0; j<3; j++)
 {
  scalar += epsilon[i][j]*epsilon[i][j];
 }
}
J2 = (I1 * I1)/6.0 - scalar/2.0;
*eps_equiv = sqrt(- J2);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of s2ml_aequistrain */

/*!----------------------------------------------------------------------
\brief  routine for calculation of macro-B-operator-values at submesh-element
GP's and of macro-strains at submesh-element GP's

<pre>                                                             ah 06/04
This routine

</pre>
\param **bopma            DOUBLE      (O)   Macro-B-operator at submesh GP
\param  *strainma         DOUBLE      (O)   macro strain at submesh GP
\param  *functmi          DOUBLE      (I)   micro quad4-Ansatzfunc at submesh GP
\param  *actsmele         ELEMENT     (I)   actual submesh element

\return void

*----------------------------------------------------------------------*/

void s2ml_bopstrainma(DOUBLE     **bopma,     /* macro B-operator at GP       */
                      DOUBLE      *strainma,  /* macro strain at GP           */
                      DOUBLE      *functmi,   /* micro quad4-Ansatzfunc at GP */
                      ELEMENT     *actsmele)  /* actual submesh element       */
{
INT      nodema,nodemi,j,node_start;       /* loopers */
DOUBLE   derivma_xy[2][4];      /* derivatives of macro-Ansatzfunctions */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_bopstrainma");
#endif
/*--------------------------------------------------- Initialisation ---*/
/*for (nodemi=0; nodemi<actsmele->numnp; nodemi++)
  for (nodema=0; nodema<4; nodema++)
    for (j=0; j<2; j++)
     printf("Bop (smnode/manode/xy) %f \n",actsmele->node[nodemi]->sm_macroinfo->derivxy_ma[j][nodema]);
*/
for (j=0; j<2; j++)
{
  for (nodema=0; nodema<4; nodema++)
  {
      derivma_xy[j][nodema] = 0.0;
  }
}
for (j=0; j<4; j++)
{
  strainma[j] = 0.0;
}
/*-- macro B-operator in submesh = Nprime*(sm nodal value of Bmacro) ---*/
for (j=0; j<2; j++)
{
  for (nodema=0; nodema<4; nodema++)
  {
    for (nodemi=0; nodemi<actsmele->numnp; nodemi++)
    {
      derivma_xy[j][nodema] += functmi[nodemi] * actsmele->node[nodemi]->sm_macroinfo->derivxy_ma[j][nodema];
    }
  }
}
for (nodema=0; nodema<4; nodema++)
{
  node_start = nodema*2;
  bopma[0][node_start]   = derivma_xy[0][nodema];
  bopma[1][node_start+1] = derivma_xy[1][nodema];
  bopma[2][node_start]   = derivma_xy[1][nodema];
  bopma[2][node_start+1] = derivma_xy[0][nodema];
} /* end of loop over macronodes */

/*- macro strain in submesh = Nprime*(sm nodal value of macrostrain) ---*/

for (j=0; j<4; j++)
{
  for (nodemi=0; nodemi<actsmele->numnp; nodemi++)
  {
    strainma[j] += functmi[nodemi] * actsmele->node[nodemi]->sm_macroinfo->strain_ma[j];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of s2ml_bopstrainma */
/*----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief  routine for calculation of macro-B-operator-values at submesh-element
GP's and of macro-strains at submesh-element GP's

<pre>                                                             ah 06/04
This routine

</pre>
\param **bopmi            DOUBLE      (O)   micro-B-operator at submesh GP
\param  *strainmi         DOUBLE      (O)   micro strain at submesh GP
\param  *actsmele         ELEMENT     (I)   actual submesh element
\param  *actmaele         ELEMENT     (I)   actual macro element
\param   nue              DOUBLE      (I)   poisson ratio

\return void

*----------------------------------------------------------------------*/
void s2ml_strainmi(DOUBLE     **bopmi,     /* micro B-operator at GP    */
                   DOUBLE      *strainmi,  /* micro strain at GP        */
                   ELEMENT     *actsmele,  /* actual submesh element    */
                   ELEMENT     *actmaele,  /* actual macro element    */
                   DOUBLE       nue)       /* poisson ratio             */
{
INT      i,smnode,node_start; /* loopers */
INT      Id;                  /* Id of submesh nodes */
DOUBLE   disd_mi[4];          /* micro-displacement derivatives */
DOUBLE   d_mi[9][2];          /* micro-displacements at submesh-nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_strainmi");
#endif

/*--- get by use of submesh-element-node-Id its micro-displ.solution ---*/
for (smnode=0; smnode<actsmele->numnp; smnode++)
{
  Id = actsmele->node[smnode]->Id;
  d_mi[smnode][0] = actmaele->e.w1->sm_nodaldata[Id].displ_mi[0];
  d_mi[smnode][1] = actmaele->e.w1->sm_nodaldata[Id].displ_mi[1];
}
/*--------------------------------------------------- initialisation ---*/
for (i=0; i<4; i++) disd_mi[i] = 0.0;
/*------------------------ micro strain = sum (Bmicro * displ.micro) ---*/
for (smnode=0; smnode<actsmele->numnp; smnode++)
{
  node_start = smnode*2;
  disd_mi[0] += bopmi[0][node_start]   * d_mi[smnode][0];
  disd_mi[1] += bopmi[1][node_start+1] * d_mi[smnode][1];
  disd_mi[2] += bopmi[2][node_start]   * d_mi[smnode][0];
  disd_mi[3] += bopmi[2][node_start+1] * d_mi[smnode][1];
} /* end of loop over nodes */

strainmi[0] = disd_mi[0];
strainmi[1] = disd_mi[1];
strainmi[2] = disd_mi[2] + disd_mi[3];
/*------------------------------------------------------ strain z-z ---*/
if(actsmele->e.w1->wtype == plane_strain)  strainmi[3] = 0.0;
if(actsmele->e.w1->wtype == plane_stress)  strainmi[3] = -(nue*(strainmi[0]+strainmi[1]))/(1.0-nue);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of s2ml_strainmi */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief  routine for calculation of displacement jumps over interface-thickness
due to macro and micro displacements at a submesh-interface-elements GP

<pre>                                                             ah 06/04
This routine

</pre>
\param  *actsmele         ELEMENT     (I)   actual submesh element
\param  *actmaele         ELEMENT     (I)   actual macro element
\param **bopmi            DOUBLE      (I)   micro-B-operator at submesh GP
\param **bopma            DOUBLE      (I)   macro-B-operator at submesh GP
\param  *jumpuma          DOUBLE      (O)   displacement jump (macro) at submesh GP
\param  *DELTAjumpuma     DOUBLE      (O)   incre displacement jump (macro) at submesh GP
\param  *jumpumi          DOUBLE      (O)   displacement jump (micro) at submesh GP
\param  *DELTAjumpumi     DOUBLE      (O)   incre displacement jump (micro) at submesh GP

\return void

*----------------------------------------------------------------------*/
void s2ml_ifbopmajumpu(ELEMENT  *actsmele,    /* actual submesh element    */
                       ELEMENT  *actmaele,    /* actual macroe element     */
                       DOUBLE  **bopmi,       /* micro B-operator at GP    */
                       DOUBLE  **bopma,       /* macro B-operator at GP    */
                       DOUBLE   *jumpuma,     /* displacement jump (macro) */
                       DOUBLE   *DELTAjumpuma,/* incre disp jump (macro)   */
                       DOUBLE   *jumpumi,     /* displacement jump (micro) */
                       DOUBLE   *DELTAjumpumi)/* incre disp jump (micro)   */
{
INT     i,j,smnode,manode;     /* loopers */
INT     ID;                    /* ID of element-submeshnode */
DOUBLE  dis[16];               /* displacement, (used for macro and micro) */
DOUBLE  DELTAdis[16];          /* incremental displacement (used for macro and micro) */
/*DOUBLE  functma[4][8];*/   /* macroansatzfunctions of macronode i, evaluated at smnode j */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_ifbopmajumpu");
#endif
/*---------------------------------------------------- anschaulicher ---*/
# if 0
  for (manode=0; manode<4; manode++)
  {
    for (smnode=0; smnode<actsmele->numnp; smnode++)
    {
       functma[manode][smnode] = actsmele->node[smnode]->sm_macroinfo->funct_ma[manode];
    }
  }
  for(i=0; i<2; i++)
  {
    for(manode=0; manode<4; manode++)
    {
      for(smnode=0; smnode<actsmele->numnp; smnode++)
      {
        bopma[i][2*manode]   += bopmi[i][2*smnode]   * functma[manode][smnode];
        bopma[i][2*manode+1] += bopmi[i][2*smnode+1] * functma[manode][smnode];
      }
    }
  }
# endif
/*--------------------------------- calculate B-operator (schneller) ---*/
for(i=0; i<2; i++)
{
  for(manode=0; manode<4; manode++)
  {
    for(smnode=0; smnode<actsmele->numnp; smnode++)
    {
      bopma[i][2*manode]   += bopmi[i][2*smnode]   * actsmele->node[smnode]->sm_macroinfo->funct_ma[manode];
      bopma[i][2*manode+1] += bopmi[i][2*smnode+1] * actsmele->node[smnode]->sm_macroinfo->funct_ma[manode];
    }
  }
}
/*-- get nodal displ.+ increm.displ (=Delta d + sum delta d) [macro] ---*/

for (manode=0; manode<4; manode++)
{
  dis[2*manode]        = actmaele->node[manode]->sol.a.da[0][0];
  dis[2*manode+1]      = actmaele->node[manode]->sol.a.da[0][1];
  DELTAdis[2*manode]   = actmaele->node[manode]->sol_increment.a.da[0][0];
  DELTAdis[2*manode+1] = actmaele->node[manode]->sol_increment.a.da[0][1];
}
/*- calculate tangential and normal macro-displacementjump at interf ---*/
for (i=0; i<2; i++)
{
  jumpuma[i]      = 0.0;
  DELTAjumpuma[i] = 0.0;
  for (j=0; j<8; j++)
  {
    jumpuma[i]      += bopma[i][j] * dis[j];
    DELTAjumpuma[i] += bopma[i][j] * DELTAdis[j];
  }
}
/*-- get nodal displ.+ increm.displ (=Delta d + sum delta d) [micro] ---*/
for(smnode=0; smnode<actsmele->numnp; smnode++)
{
  ID = actsmele->node[smnode]->Id;
  dis[2*smnode]       = actmaele->e.w1->sm_nodaldata[ID].displ_mi[0];
  dis[2*smnode+1]     = actmaele->e.w1->sm_nodaldata[ID].displ_mi[1];
  DELTAdis[2*smnode]  = actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[0];
  DELTAdis[2*smnode+1]= actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[1];
}
/*- calculate tangential and normal micro-displacementjump at interf ---*/
for (i=0; i<2; i++)
{
  jumpumi[i]      = 0.0;
  DELTAjumpumi[i] = 0.0;
  for (j=0; j<2*actsmele->numnp; j++)
  {
    jumpumi[i]      += bopmi[i][j] * dis[j];
    DELTAjumpumi[i] += bopmi[i][j] * DELTAdis[j];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of s2ml_ifbopmajumpu */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
#endif
