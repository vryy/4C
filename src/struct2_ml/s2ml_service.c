/*!----------------------------------------------------------------------
\file
\brief contains the routine w1ml_bopstraintonode

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
#include "s2ml.h"
#include "s2ml_prototypes.h"

/*! 
\addtogroup MLSTRUCT 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine for calculation of macro-B-operator and macro-strains
at the nodes of the submesh

\param *actsmfield       FIELD      (I)    submesh field
\param *ele              ELEMENT    (I)    actual makroelement
\param  nue              DOUBLE     (I)    poisson ratio of macro material
\param  nue              DOUBLE     (I)    poisson ratio of macro material
\param  init             INT        (I)    init

\return void                                               

*----------------------------------------------------------------------*/
void s2ml_bopstraintonode(FIELD       *actsmfield, /*  submesh field           */
                          ELEMENT     *ele,        /* actual makroelement           */
                          DOUBLE       nue,       /* poisson ratio of macro material           */
                          INT          init) 
{
INT       cornernode;              /* Elementnodenumber of left-down Makroelement-node */
INT       i,node;                  /* some loopers */
INT       ielma;                   /* number of nodal points of this Macroelement */
DOUBLE    translation[2];          /* coordinates of left-down Makroelement-node */
DOUBLE    xi_macro=0.0;            /* coordinates of submesh-node in macroelement-local-COOS */
DOUBLE    eta_macro=0.0;           /* coordinates of submesh-node in macroelement-local-COOS */
DOUBLE    x_P,y_P;                 /* coordinates of submesh-node in X-Y System */
DOUBLE    functma[4];              /* Macroansatzfunctions */
DOUBLE    disdma[5];               /* Makro-displacement derivatives */
DOUBLE    strainma[4];             /* Macro-strain */
DOUBLE    detma;                   /* Macro-jacobi-determinante */
const DOUBLE    tol = 1E-05;

static ARRAY    bopma_a;           /*  Makro-B-operator */   
static DOUBLE **bopma;       
static ARRAY    derivma_a;         /* Derivatives of Macroansatzfunctions */  
static DOUBLE **derivma;       
static ARRAY    xjmma_a;           /* Jacobian Matrix (Macro) */  
static DOUBLE **xjmma;       

DOUBLE    A_1[2],A_2[2],A_3[2];    /* some interim values */

struct _NODE       *actsmnode;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_bopstraintonode");
#endif

/*----------------------------------------------------------------------*/
if (init==1)
{
  bopma   = amdef("bopma"    ,&bopma_a   ,3,(2*MAXNOD_WALL1),"DA");           
  derivma = amdef("derivma"  ,&derivma_a ,2,MAXNOD_WALL1 ,"DA");       
  xjmma   = amdef("xjmma"    ,&xjmma_a   ,2,2    ,"DA");           
  goto end;
}

/*--------------- find macro-node with smallest x and y coordinate ---*/
translation[0] = ele->node[0]->x[0];
translation[1] = ele->node[0]->x[1];
for (i=0; i<ele->numnp; i++)
{
  if(ele->node[i]->x[0] <= translation[0] && ele->node[i]->x[1] <= translation[1])
  {
    translation[0] = ele->node[i]->x[0];
    translation[1] = ele->node[i]->x[1];
    cornernode = i;
  }
}

/*------ Inversion of ansatzfunctions for RECTENGULAR macroelements ---*/
for (i=0; i<2; i++)
{
  A_1[i] = ele->node[0]->x[i] + ele->node[1]->x[i] + ele->node[2]->x[i] + ele->node[3]->x[i];
  A_2[i] = ele->node[0]->x[i] - ele->node[1]->x[i] - ele->node[2]->x[i] + ele->node[3]->x[i];
  A_3[i] = ele->node[0]->x[i] + ele->node[1]->x[i] - ele->node[2]->x[i] - ele->node[3]->x[i];
}

/*---------------------------------------------------------------------*/
ielma   = ele->numnp;

for (node=0; node<actsmfield->dis[0].numnp; node++)
{
  actsmnode = &(actsmfield->dis[0].node[node]);
  /*---------------------------- x and y coordinate of submesh nodes ---*/
  x_P       = actsmnode->x[0]+translation[0];
  y_P       = actsmnode->x[1]+translation[1];
  /*----- xi and eta coordinate of submesh nodes in local makro COOS ---*/
  if (fabs(A_3[0])<tol && fabs(A_2[1])<tol)  /* xi // x and eta // y ---*/
  {
    xi_macro  = (4 * x_P - A_1[0])/A_2[0];
    eta_macro = (4 * y_P - A_1[1])/A_3[1];
  }
  else if(fabs(A_2[0])<tol && fabs(A_3[1])<tol)/* xi // y and eta // x -*/
  {
    xi_macro  = (4 * y_P - A_1[1])/A_2[1];
    eta_macro = (4 * x_P - A_1[0])/A_3[0];
  }
  else
  {
    dserror("macroelements have to be rectengular");
  }
  /*-------------------  Makro shape functions and their derivatives ---*/
  w1_funct_deriv(functma,derivma,xi_macro,eta_macro,ele->distyp,1);
  /*---------------------------------- compute Makro-jacobian matrix ---*/       
  w1_jaco (derivma,xjmma,&detma,ele,ielma);                         
  /*------------------------------------- calculate Makro-operator B ---*/
  amzero(&bopma_a);
  w1_bop(bopma,derivma,xjmma,detma,ielma);
  /*----------------------------------------- calculate Makro-strain ---*/
  w1_disd(ele,bopma,NULL,NULL,ele->e.w1->wtype,disdma);                  
  strainma[0] = disdma[0];
  strainma[1] = disdma[1];
  strainma[2] = disdma[2] + disdma[3];
  strainma[3] = 0.0;
  if (ele->e.w1->wtype == plane_stress) 
  {
   strainma[3] = - (nue*(strainma[0]+strainma[1]))/(1.0 - nue);
  }
  /*------------------------------------------- save to submesh-node ---*/
  for (i=0; i<4; i++)
    actsmnode->sm_macroinfo->strain_ma[i]     = strainma[i]; 
  for (i=0; i<ielma; i++)
  {
    actsmnode->sm_macroinfo->funct_ma[i]      = functma[i]; 
    actsmnode->sm_macroinfo->derivxy_ma[0][i] = bopma[0][2*i]; 
    actsmnode->sm_macroinfo->derivxy_ma[1][i] = bopma[1][2*i+1]; 
  }
}

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of s2ml_bopstraintonode */


/*!----------------------------------------------------------------------
\brief  routine for calling the material law from a submeshelement 

<pre>                                                           ah 07/04
This routine 

</pre>
\param  *actmaele        ELEMENT     (I)    actual macro element
\param  *actsmele        ELEMENT     (I)    actual submesh element
\param  *actsmmat        MATERIAL    (I)    actual submesh material
\param  *strain          DOUBLE      (I)    total strain at GP
\param   ip              INT         (I)    IP of GP
\param  *stress          DOUBLE      (O)    stress at GP
\param **D               DOUBLE      (O)    material tangent at GP
\param   istore          INT         (I)    is it update step?

\return void                                               

*----------------------------------------------------------------------*/
void s2ml_callmat(ELEMENT     *actmaele,   /*  actual macro element   */
                  ELEMENT     *actsmele,   /*  actual submeshelement  */
                  MATERIAL    *actsmmat,   /*  actual sm_material     */
                  DOUBLE      *strain,     /* total strain at GP      */
                  INT          ip,         /* ip of GP                */
                  DOUBLE      *stress,     /*  stress                 */
                  DOUBLE     **D,          /*  tangent                */
                  INT          istore)     /*  is it update step?     */
{
INT smallscale = 1;
INT newval = 0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_callmat");
#endif
/*----------------------------------------------------------------------*/
switch(actsmmat->mattyp)
{
  case m_stvenant:/*------------------------------- linear elastic ---*/
    dserror("st.vernant not jet implemented for multiscale");
  break;
  case m_damage:   /*------------------------------------- CDM law ---*/
    w1_mat_damage(actsmmat->m.damage->youngs,
                  actsmmat->m.damage->possionratio,
                  actsmmat->m.damage->Equival,
                  actsmmat->m.damage->Damtyp,
                  actsmmat->m.damage->Kappa_0,
                  actsmmat->m.damage->Kappa_m,
                  actsmmat->m.damage->Alpha,
                  actsmmat->m.damage->Beta,
                  actsmmat->m.damage->k_fac,
                  actmaele,
                  plane_strain,
                  NULL,
                  NULL,
                  NULL,
                  ip,
                  stress,
                  D,
                  istore,
                  newval,
                  actsmele,
                  smallscale,
                  strain);
  break;
  default:
    dserror("unknown type of material law");
  break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of s2ml_callmat */

/*!----------------------------------------------------------------------
\brief  routine for integration of element stiffness(has not to be quadratic!)

<pre>                                                            ah 06/04
This routine 

</pre>
\param **stiffmatrix    DOUBLE       (O)   element stiffnes matrix 
\param **left_bop       DOUBLE       (I)   left B-operator (to be transp)
\param **right_bop      DOUBLE       (I)   right B-operator
\param **D              DOUBLE       (I)   material tangent
\param   fac            DOUBLE       (I)   integration factor
\param   left_nd        INT          (I)   number ele-dof left vec
\param   right_nd       INT          (I)   number ele-dof right vec
\param   numeps         INT          (I)   number of strains-compon

\return void                                               

*----------------------------------------------------------------------*/
void s2ml_cal_stiff(DOUBLE   **stiffmatrix, /* element stiffnes matrix */
                    DOUBLE   **left_bop,    /* left B-operator (to be transp) */
                    DOUBLE   **right_bop,   /* right B-operator        */
                    DOUBLE   **D,           /* material tangent        */
                    DOUBLE   fac,           /* integration factor      */
                    INT      left_nd,       /* number ele-dof left vec */
                    INT      right_nd,      /* number ele-dof right vec*/
                    INT      numeps)        /* number of strains-compon*/
{
INT            i, j, k, l, m;
DOUBLE         dum;
DOUBLE         DB_right[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_cal_stiff");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<right_nd; j++)
{
  for (k=0; k<numeps; k++)
  {
    DB_right[k]= 0.0;
    for(l=0; l<numeps; l++)
    {
      DB_right[k] = DB_right[k] + D[k][l] * right_bop[l][j] * fac;
    }  
  }
  for(i=0; i<left_nd; i++)
  {
    dum = 0.0;
    for(m=0; m<numeps; m++)
    {
      dum = dum + left_bop[m][i] * DB_right[m];
    }
    stiffmatrix[i][j] = stiffmatrix[i][j] + dum;
  }/* end loop over dof left (i) */
}/* end loop over dof right (j) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of s2ml_cal_stiff */
/*----------------------------------------------------------------------*/


/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
