/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains a lot of functions for the axishell element


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/


#ifdef D_AXISHELL

#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"


/*!
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/




/*-----------------------------------------------------------------------*/
/*!
  \brief initialization routine for the axisymmetric shell element

  This routine acts according to the action and either initializes the element
  or computes the linear stiffness matrix, the stresses or the right hand side
  vector.

  \param *actpart      PARTITION   (i)   my partition

  \return void
  \sa calling:   saxiintg(); called by: axishell();

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxiinit(
    PARTITION *actpart
    )
{

  INT          i;
  ELEMENT     *actele;
  SAXI_DATA      data;

  for (i=0; i<actpart->pdis[0].numele; i++)
  {
    actele = actpart->pdis[0].element[i];
    if (actele->eltyp != el_axishell) continue;

    /* init integration points */
    saxiintg(&data);

    /* allocate the space for stresses */
    am4def("stress_GP",&(actele->e.saxi->stress_GP),1,5,1,0,"D3");
    am4def("stress_ND",&(actele->e.saxi->stress_ND),1,5,1,0,"D3");
  }

#ifdef DEBUG
  dstrc_enter("saxiinit");
#endif

  return;
} /* end of saxiinit */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculate operator matrix at point xsi=s/l

  This routine calcuates the operator matrix B at the given point xsi=s/l
  for an axisymmetric shell element.

  \param   B     DOUBLE[][] (o)  the calculated operator matrix
  \param   xsi   DOUBLE     (i)  the integration point where B is computed
  \param   r     DOUBLE     (i)  the radius of the integration point
  \param   dl    DOUBLE     (i)  the length of the current element
  \param   cosa  DOUBLE     (i)  the cosine of the angle of the current element
  \param   sina  DOUBLE     (i)  the sine of the angle of the current element

  \return void
  \sa calling: ---; called by: saxi_static_ke(), saxi_cal_stress()

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxi_B(
    DOUBLE      B[4][7],
    DOUBLE      xsi,
    DOUBLE      r,
    DOUBLE      dl,
    DOUBLE      cosa,
    DOUBLE      sina
    )
{

  INT i,j;

#ifdef DEBUG
  dstrc_enter("saxi_B");
#endif

  for (i=0; i<4; i++)
  {
    for (j=0; j<7; j++)
    {
      B[i][j] = 0.0;
    }
  }

  /* compute operator B  */
  /* compare with Schalen-Skriptum S.11.8  */
  B[0][0] = (-3.0+4.0*xsi)/dl;
  B[0][3] = (-1.0+4.0*xsi)/dl;
  B[1][0] = (1.0-3.0*xsi+2.0*xsi*xsi)*cosa/r;
  B[1][1] = (1.0-3.0*xsi*xsi+2.0*xsi*xsi*xsi)*sina/r;
  B[1][2] = (1.0-2.0*xsi+xsi*xsi)*sina*dl*xsi/r;
  B[1][3] = (-xsi+2.0*xsi*xsi)*cosa/r;
  B[1][4] = (3.0*xsi*xsi-2.0*xsi*xsi*xsi)*sina/r;
  B[1][5] = (-1.0+xsi)*xsi*xsi*sina*dl/r;
  B[2][1] = (-6.0+12.0*xsi)/dl/dl;
  B[2][2] = (-4.0+6.0*xsi)/dl;
  B[2][4] = -B[2][1];
  B[2][5] = (-2.0+6.0*xsi)/dl;
  B[3][1] = 6.0*xsi*(-1.0+xsi)*cosa/(dl*r);
  B[3][2] = (1.0-4.0*xsi+3.0*xsi*xsi)*cosa/r;
  B[3][4] = -B[3][1];
  B[3][5] = xsi*(-2.0+3.0*xsi)*cosa/r;
  B[0][6] = (4.0-8.0*xsi)/dl;
  B[1][6] = 4.0*xsi*(1.0-xsi)*cosa/r;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_B */




/*-----------------------------------------------------------------------*/
/*!
  \brief coordinates and weight factors for numerical integration

  This routine  gives the coordinates and weight factors for numerical
  integration of a axisymmetric shell element. Gauss-Integration with 5 GP
  from 0 to 1 is used here.

  \param *data   SAXI_DATA  (o)   structure containing the coordinates and
                                  weighting factors

  \return void
  \sa calling:   ---;
      called by: saxiinit, saxi_static_ke, saxi_cal_fext, saxi_cal_stress;

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxiintg(
    SAXI_DATA   *data
    )
{

#ifdef DEBUG
  dstrc_enter("saxiintg");
#endif

  /* gauss sampling points in xsi-direction */
  data->xgr[0] = 0.046910077030668;
  data->xgr[1] = 0.2307653449471585;
  data->xgr[2] = 0.5;
  data->xgr[3] = 1.0-data->xgr[1];
  data->xgr[4] = 1.0-data->xgr[0];

  /* weighting factors */
  data->wgt[0] = 0.236926885056189 / 2.0;
  data->wgt[1] = 0.478628670499366 / 2.0;
  data->wgt[2] = 0.568888888888889 / 2.0;
  data->wgt[3] = data->wgt[1];
  data->wgt[4] = data->wgt[0];

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxiintg */




/*-----------------------------------------------------------------------*/
/*!
  \brief Computation of the rhs of the axialsymmetric shell element

  This routine computes the right hand side of the axialsymmetric shell element
  for loads acting on the actual element. The element loads used here are:

  <pre>
  /-----------------------------------------------------------------------*
  |  From Input-file:                                                     |
  |    px...load parallel to xsi direction (meridian direction)           |
  |    pw...load normal to xsi direction (normal direction)               |
  |    pv...load in vertical direction                                    |
  |    ph...load in horizontal direction                                  |
  |  Computed:                                                            |
  |    qx...resultant load in xsi direction (meridian direction)          |
  |    pw...resultant load normal to xsi direction (normal direction)     |
  *-----------------------------------------------------------------------/

  /----------------- SIGN CONVENTION FOR ELEMENT FORCES ------------------*
  |                                                                       |
  |                                   node i                              |
  |                +pv    +pw           o          |                      |
  |                 ^        \         /           |                      |
  |                 |          \      /            |                      |
  |                 |            \   /             |                      |
  |                 |              \/              |                      |
  |                 |------>+ph   //               |                      |
  |                              //                |                      |
  |                             //                 |                      |
  |                            //                  |                      |
  |                         +px/                   |                      |
  |                           /                    |                      |
  |                          o                     |                      |
  |                       node j                   rotational axis        |
  |    +y                                                                 |
  |    | global                                                           |
  |    | coordinate                                                       |
  |    | system                                                           |
  |    |----->+x                                                          |
  |                                                                       |
  *-----------------------------------------------------------------------/
  </pre>

  \param *ele       ELEMENT   (i)  my element
  \param *data      SAXI_DATA (i)  structure containing gaussian point and weight
  \param *loadvec   DOUBLE    (i)  global element load vector
  \param  init      INT       (i)  flag init == 1 : initialization
                                        init != 1 : integration

  \return void
  \sa calling:   saxiintg(); called by: axishell();

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxi_eleload(
    ELEMENT    *ele,
    SAXI_DATA  *data,
    DOUBLE     *loadvec,
    INT	      init
    )
{

  INT     i,lxsi,nixsi;
  DOUBLE  xsi,fac;
  DOUBLE  ri,rj,dr,dz,dl,cosa,sina,r,px,pw,pv,ph,qx,qw;
  DOUBLE  px_ij[2],pw_ij[2],pv_ij[2],ph_ij[2];    /* 0..node i, 1..node j */
  DOUBLE  loadvec_7[7];     /*local  loadvector before static condensation*/
  DOUBLE  loadvec_local[6];   /*local loadvector after static condensation*/

  static ARRAY eload_a;
  static DOUBLE **eload;                    /* static element load vector */
  DOUBLE *statcond;

#ifdef DEBUG
  dstrc_enter("saxi_eleload");
#endif

  /* init phase (init==1) */
  if (init==1)
  {
    eload   = amdef("eload"  ,&eload_a,MAXDOFPERNODE,MAXNOD_AXISHELL,"DA");
    goto end;
  }

  /* uninit phase (init==-1) */
  else if (init==-1)
  {
    amdel(&eload_a);
    goto end;
  }


  /* initialize eload */
  amzero(&eload_a);

  /* integration parameters */
  saxiintg(data);

  /* some element infos*/
  ri = ele->node[0]->x[0];
  rj = ele->node[1]->x[0];
  dr = rj - ri;                                 /* delta r */
  dz = ele->node[1]->x[1] - ele->node[0]->x[1]; /* delta z */
  dl = sqrt(dr*dr+dz*dz);                       /* delta l */
  cosa = dr/dl;
  sina = dz/dl;

  /* general case for a linear distribution of the loads in the element */
  /* Value of the loads at node i and j */
  pv_ij[0] = ele->e.saxi->pv[0];
  pv_ij[1] = ele->e.saxi->pv[1];
  ph_ij[0] = ele->e.saxi->ph[0];
  ph_ij[1] = ele->e.saxi->ph[1];
  px_ij[0] = ele->e.saxi->px[0];
  px_ij[1] = ele->e.saxi->px[1];
  pw_ij[0] = ele->e.saxi->pw[0];
  pw_ij[1] = ele->e.saxi->pw[1];

  for (i=0; i<7; i++)
  {
    loadvec_7[i] = 0.0;
  }
  nixsi   = 5; /* number of Gauss-points */

  /* integration loop */
  for (lxsi=0; lxsi<nixsi; lxsi++)
  {
    /* gaussian point, weight and thickness at it */
    xsi = data->xgr[lxsi];
    fac = data->wgt[lxsi];
    r = ele->node[1]->x[0] + cosa*xsi*dl;

    pv = pv_ij[0]+xsi/1.0*(pv_ij[1]-pv_ij[0]);
    ph = ph_ij[0]+xsi/1.0*(ph_ij[1]-ph_ij[0]);
    px = px_ij[0]+xsi/1.0*(px_ij[1]-px_ij[0]);
    pw = pw_ij[0]+xsi/1.0*(pw_ij[1]-pw_ij[0]);

    qx = px - ph*cosa + pv*sina;
    qw = pw + ph*sina - pv*cosa;

    loadvec_7[0] += dl*r*qx*(-1.0+3.0*xsi-2.0*xsi*xsi)        *fac;
    /*zusaetzlicher term bei loadvec_7[0] ?*/
    loadvec_7[1] += dl*r*qw*(-1.0+3.0*xsi*xsi-2.0*xsi*xsi*xsi)*fac;
    loadvec_7[2] += dl*r*qw*dl*xsi*(-1.0+2.0*xsi-xsi*xsi)     *fac;
    loadvec_7[3] += dl*r*qx*xsi*(1.0-2.0*xsi)                 *fac;
    /*zusaetzlicher term bei loadvec_7[3] ?*/
    loadvec_7[4] += dl*r*qw*xsi*xsi*(-3.0+2.0*xsi)            *fac;
    loadvec_7[5] += dl*r*qw*dl*xsi*xsi*(1.0-xsi)              *fac;
    loadvec_7[6] += dl*r*qx*4.0*xsi*(-1.0+xsi)                *fac;
  }   /* end of integration loop */

  ele->e.saxi->Sk *= loadvec_7[6]; /*Sk = S^0_7/k77 (siehe Skript S.11.11)*/

  /* reduction of the rhs-vector from 7 to 6 entries */
  statcond = ele->e.saxi->statcond;
  for (i=0; i<6; i++)
  {
    /*auf statcond soll schon -k_17/k77 stehen!*/
    loadvec_local[i] = loadvec_7[i] + statcond[i] * loadvec_7[6];/*S_ausFlaechenlasten)*/
    loadvec_local[i] = -loadvec_local[i];/* K*u=(S_ausKnotenkraeften-S_ausFlaechenlasten)*/
  }

  /* equivalent loads in global coordinates */
  if (ele->node[0]->gnode->ondesigntyp==ondnode && ele->node[0]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    loadvec[0] = loadvec_local[0];
    loadvec[1] = loadvec_local[1];
    loadvec[2] = loadvec_local[2];
  }
  else
  {
    /* global cos */
    loadvec[0] =  loadvec_local[0]*cosa + loadvec_local[1]*sina;/*load in horizontal direction*/
    loadvec[1] =  loadvec_local[0]*sina - loadvec_local[1]*cosa;/*load in vertical direction*/
    loadvec[2] = -loadvec_local[2];                             /*moment */
  }


  if (ele->node[1]->gnode->ondesigntyp==ondnode && ele->node[1]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    loadvec[3] = loadvec_local[3];
    loadvec[4] = loadvec_local[4];
    loadvec[5] = loadvec_local[5];
  }
  else
  {
    /* global cos */
    loadvec[3] =  loadvec_local[3]*cosa + loadvec_local[4]*sina;/*load in horizontal direction*/
    loadvec[4] =  loadvec_local[3]*sina - loadvec_local[4]*cosa;/*load in vertical direction*/
    loadvec[5] = -loadvec_local[5];                             /*moment */
  }


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_eleload */




/*-----------------------------------------------------------------------*/
/*!
  \brief contains the material law for the axisymmetric shell element

  This routine computes the constitutive  matrix of the axisymmetric shell
  element

  \param *mat       STVENANT   (i)   my material model
  \param  D         DOUBLE[][] (i/o) the constitutive matrix (4/4)
  \param  thick     DOUBLE     (i)   thickness of the element

  \return void
  \sa calling:   ---; called by: saxi_static_ke(), saxi_cal_stress();

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxi_mat(
    STVENANT  *mat,
    DOUBLE     D[4][4],
    DOUBLE     thick
    )
{

  INT i,j;
  DOUBLE E, mu, fac;

#ifdef DEBUG
  dstrc_enter("saxi_mat");
#endif

  E  = mat->youngs;
  mu = mat->possionratio;
  fac= E*thick/(1.0-mu*mu);

  /* Zero out all entries of the 4/4 constitutive matrix */
  for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      D[i][j] = 0.0;
    }
  }

  D[0][0] = fac;
  D[0][1] = fac*mu;
  D[1][0] = D[0][1];
  D[1][1] = D[0][0];
  D[2][2] = fac*thick*thick/12.0;
  D[2][3] = fac*mu*thick*thick/12.0;
  D[3][2] = D[2][3];
  D[3][3] = D[2][2];

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_mat */




/*-----------------------------------------------------------------------*/
/*!
  \brief Computation of the stresses for the axialsymmetric shell element

  This routine computes the stresses for the axialsymmetric shell element
  after computation of the displacements.

  \param *ele       ELEMENT   (i)   my element
  \param *data      SAXI_DATA (i)   structure containing gaussian point and weight
  \param *mat       MATERIAL  (i)   my material
  \param  init      INT       (i)   flag init == 1 : initialization
                                         init != 1 : integration

  \return void
  \sa calling:   saxiintg, saxi_B, saxi_mat; called by: axishell();

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxi_cal_stress(
    ELEMENT   *ele,
    SAXI_DATA *data,
    MATERIAL  *mat,
    INT        init
    )
{

  INT                 i,j,k;                    /* some loopers */
  INT                 nixsi;                    /* num GP in xsi direction*/
  INT                 lxsi;                     /* loopers over GP */
  INT                 iel;                      /* numnp to this element */
  INT                 nd;
  INT                 diff,max;

  DOUBLE              thick;                    /* thickness */
  DOUBLE              xsi;                      /* GP-coords */
  DOUBLE              r,dr,dz,dl,cosa,sina;     /* infos about the element*/
  DOUBLE              fac;                      /* weights at GP */

  DOUBLE              disp_global[7];/* vector of global nodal displ. */
  DOUBLE              disp_local[7]; /* vector of local nodal displa. */
  DOUBLE              S[4][7];       /* stress matrix S=D*B(xsi=0.5) */
  DOUBLE              n[5];          /* Schnittkraftvektor */
  DOUBLE              T[6][6];       /* Transformation matrix */
  DOUBLE              Bs[4][7];      /* Ableitung der Matrix B nach xsi */
  DOUBLE              DBs[7];        /* Multiplication of D and Bs */

  static DOUBLE       D[4][4];       /* material tensor */
  static DOUBLE       B[4][7];       /* B-operator */

  DOUBLE ***gp_stress;   /* pointer to array of stresses at the middle GP */
  DOUBLE * statcond;

#ifdef DEBUG
  dstrc_enter("saxi_cal_stress");
#endif

  /* init phase */
  if (init==1)
    goto end;

  /* integration parameters */
  saxiintg(data);

  /* integration parameters */
  nixsi   = 5; /* number of Gauss-points */
  iel     = ele->numnp;
  nd      = 2 * iel;

  /* datas about the element under consideration */
  dr = ele->node[1]->x[0] - ele->node[0]->x[0]; /* delta r */
  dz = ele->node[1]->x[1] - ele->node[0]->x[1]; /* delta z */
  dl = sqrt(dr*dr+dz*dz);                       /* delta l */
  cosa = dr/dl;
  sina = dz/dl;

  /* store node displacements (=local) on local field --*/
  /* compute  transfomation matrix T */
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      T[i][j] = 0.0;
    }
  }
  T[0][0] =  cosa;
  T[0][1] =  sina;
  T[1][0] =  sina;
  T[1][1] = -cosa;
  T[2][2] =  -1.0;
  T[3][3] =  cosa;
  T[3][4] =  sina;
  T[4][3] =  sina;
  T[4][4] = -cosa;
  T[5][5] =  -1.0;


  /* node 0 */
  if (ele->node[0]->gnode->ondesigntyp==ondnode && ele->node[0]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    disp_local[0] = ele->node[0]->sol.a.da[0][0];
    disp_local[1] = ele->node[0]->sol.a.da[0][1];
    disp_local[2] = ele->node[0]->sol.a.da[0][2];
    /* calculate global disps */
    for (i=0; i<3; i++)
    {
      disp_global[i] = 0.0;
      for (j=0; j<3; j++)
      {
        disp_global[i] += T[j][i] * disp_local[j];
      }
    }
  }
  else
  {
    /* global cos */
    disp_global[0] = ele->node[0]->sol.a.da[0][0];
    disp_global[1] = ele->node[0]->sol.a.da[0][1];
    disp_global[2] = ele->node[0]->sol.a.da[0][2];
    /* calculate local disps */
    for (i=0; i<3; i++)
    {
      disp_local[i] = 0.0;
      for (j=0; j<3; j++)
      {
        disp_local[i] += T[i][j] * disp_global[j];
      }
    }
  }

  /* node 1 */
  if (ele->node[1]->gnode->ondesigntyp==ondnode && ele->node[1]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    disp_local[3] = ele->node[1]->sol.a.da[0][0];
    disp_local[4] = ele->node[1]->sol.a.da[0][1];
    disp_local[5] = ele->node[1]->sol.a.da[0][2];
    /* calculate global disps */
    for (i=3; i<6; i++)
    {
      disp_global[i] = 0.0;
      for (j=3; j<6; j++)
      {
        disp_global[i] += T[j][i] * disp_local[j];
      }
    }
  }
  else
  {
    /* global cos */
    disp_global[3] = ele->node[1]->sol.a.da[0][0];
    disp_global[4] = ele->node[1]->sol.a.da[0][1];
    disp_global[5] = ele->node[1]->sol.a.da[0][2];
    /* calculate local disps */
    for (i=3; i<6; i++)
    {
      disp_local[i] = 0.0;
      for (j=3; j<6; j++)
      {
        disp_local[i] += T[i][j] * disp_global[j];
      }
    }
  }

  /* write global dis. to the nodes */
  ele->node[0]->sol.a.da[0][0] = disp_global[0];
  ele->node[0]->sol.a.da[0][1] = disp_global[1];
  ele->node[0]->sol.a.da[0][2] = disp_global[2];
  ele->node[1]->sol.a.da[0][0] = disp_global[3];
  ele->node[1]->sol.a.da[0][1] = disp_global[4];
  ele->node[1]->sol.a.da[0][2] = disp_global[5];

  /* write local dis. to the nodes */
  /* enlarge sol, if necessary */
  if (1 >= ele->node[0]->sol.fdim)
  {
    diff = 1 - ele->node[0]->sol.fdim;
    max  = IMAX(diff,3);
    amredef(&(ele->node[0]->sol),ele->node[0]->sol.fdim+max+1,ele->node[0]->sol.sdim,"DA");
  }
  if (1 >= ele->node[1]->sol.fdim)
  {
    diff = 1 - ele->node[1]->sol.fdim;
    max  = IMAX(diff,3);
    amredef(&(ele->node[1]->sol),ele->node[1]->sol.fdim+max+1,ele->node[1]->sol.sdim,"DA");
  }

  ele->node[0]->sol.a.da[1][0] = disp_local[0];
  ele->node[0]->sol.a.da[1][1] = disp_local[1];
  ele->node[0]->sol.a.da[1][2] = disp_local[2];
  ele->node[1]->sol.a.da[1][0] = disp_local[3];
  ele->node[1]->sol.a.da[1][1] = disp_local[4];
  ele->node[1]->sol.a.da[1][2] = disp_local[5];




  /* compute mid-side displacement */
  statcond = ele->e.saxi->statcond;
  disp_local[6]= 0.0;
  for (i=0; i<6; i++)
  {
    disp_local[6] += statcond[i]*disp_local[i];
  }
  disp_local[6] -= ele->e.saxi->Sk;


  /* loop over GP (=Auswertungsstellen) */
  for (lxsi=2; lxsi<3; lxsi++) /* Auswertung nur bei xsi=0.5 */
  {
    /* gaussian point, weight and thickness at it */
    xsi   = data->xgr[lxsi];                  /* xsi */
    fac   = data->wgt[lxsi];                  /* weight */
    r     = ele->node[0]->x[0] + cosa*xsi*dl; /* radius */
    thick = ele->e.saxi->thick[0] +           /* thickness */
      xsi/1.0*(ele->e.saxi->thick[1]-ele->e.saxi->thick[0]);

    /* calculate operator B */
    saxi_B(B,xsi,r,dl,cosa,sina);

    /* compute constitutive matrix D */
    saxi_mat(mat->m.stvenant,D,thick);

    /* compute stress matrix D*B */
    for (k=0; k<7; k++)
    {
      for (i=0; i<4; i++)
      {
        S[i][k] = 0.0;
        for (j=0; j<4; j++)
        {
          S[i][k] += D[i][j] * B[j][k];
        }
      }
    }

    /* compute schnittkraefte n */
    /* normal forces n_s, n_deta and moments m_s, m_deta */
    for (i=0; i<4; i++)
    {
      n[i] = 0.0;
      for (j=0; j<7; j++)
      {
        n[i] += S[i][j] * disp_local[j];
      }
    }

    /* Rueckrechnung der Querkraft q_s */
    /* Ermittlung der Ableitung von ms nach ds                           */
    /* Bs[i][j]...Ableitung der Matrix B nach xsi (nur 3. und 4. Zeile)  */
    for (i=0; i<4; i++)
    {
      for (j=0; j<7; j++)
      {
        Bs[i][j] = 0.0;
      }
    }
    Bs[2][1] = 12.0/dl/dl;
    Bs[2][2] = 6.0/dl;
    Bs[2][4] = -Bs[2][1];
    Bs[2][5] = 6.0/dl;
    Bs[3][1] = 6.0*(-1.0+2.0*xsi)*cosa/(dl*r) -
      6.0*cosa*xsi*(-1.0+xsi)*dr/dl/r/r;
    Bs[3][2] = (-4.0+6.0*xsi)*cosa/r -
      cosa/r/r*dr*dl*(1.0-4.0*xsi+3.0*xsi*xsi);
    Bs[3][4] = -B[3][1];
    Bs[3][5] = (-2.0+6.0*xsi)*cosa/r -
      cosa/r/r*dr*dl*(-2.0*xsi+3.0*xsi*xsi);

    /* since we have to derive B with respect to s we have to divide the
       derivative of B with respect to xsi additionally bz dl !!!        */
    for (i=0; i<4; i++)
    {
      for (j=0; j<7; j++)
      {
        Bs[i][j] /= dl;
      }
    }

    /* Multiplication of D with Bs */
    for (i=0; i<7; i++)
    {
      DBs[i] = 0.0;
      for (j=2; j<4; j++) /* since D[2][0] and D[2][1] are 0 */
      {
        DBs[i] += D[2][j] * Bs[j][i];
      }
    }

    /* Computation of the first part of q_s = DBs * disp_local */
    n[4] = 0.0;
    for (i=0; i<7; i++)
    {
      n[4] += DBs[i] * disp_local[i];
    }
    /* final solution of q_s */
    n[4] += (n[2]-n[3])*cosa/r;
  }/* end of loop over GP (=Auswertungsstellen) */


  /* store the local n-vector onto stress_GP for postprocessing */
  gp_stress = ele->e.saxi->stress_GP.a.d3;
  for (i=0; i<5; i++)
  {
    gp_stress[0][i][0] = n[i];
  }

end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_cal_stress */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculates usual stiffness matrix in total lagrangian formulation

  This routine calculates usual stiffness matrix in total lagrangian
  formulation.

  \param   s       DOUBLE[][]  (o)  element stiffness matrix
  \param   b       DOUBLE[][]  (i)  derivative operator
  \param   d       DOUBLE[][]  (i)  constitutive matrix
  \param   work1   DOUBLE[][]  (i)  working array (4,7)
  \param   work2   DOUBLE[][]  (i)  working array (7,7)
  \param   fac     DOUBLE      (i)  integration factor
  \param   r       DOUBLE      (i)  radius of current integration point
  \param   dl      DOUBLE      (i)  length of current element

  \return void
  \sa calling:  ---; caled by: saxi_static_ke();

  \author mn
  \date   05/03

 */
/*-----------------------------------------------------------------------*/
void saxi_keku(
    DOUBLE    s[7][7],
    DOUBLE    b[4][7],
    DOUBLE    d[4][4],
    DOUBLE    work1[4][7],
    DOUBLE    work2[7][7],
    DOUBLE    fac,
    DOUBLE    r,
    DOUBLE    dl
    )
{

  INT            i, j, k;

#ifdef DEBUG
  dstrc_enter("saxi_keku");
#endif

  /*  make multiplication: work1 = D * B */
  for (k=0; k<7; k++)
  {
    for (i=0; i<4; i++)
    {
      work1[i][k] = 0.0;
      for (j=0; j<4; j++)
      {
        work1[i][k] += d[i][j] * b[j][k];
      }
    }
  }

  /*  make multiplication: work2 = B^T * work1 */
  for (k=0; k<7; k++)
  {
    for (i=0; i<7; i++)
    {
      work2[i][k] = 0.0;
      for (j=0; j<4; j++)
      {
        work2[i][k] += b[j][i] * work1[j][k];
      }
    }
  }

  /* compute: estiff[i][j] += fac*r*dl/2*work2[i][j] */
  for (i=0; i<7; i++)
  {
    for (j=0; j<7; j++)
    {
      s[i][j] += fac*r*dl*work2[i][j]/1.0;
      /* im fortran code wird durch 2 dividiert, dieser faktor ist bei
         mir schon in den Gewichten der numerischen integration */
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxi_keku */



/*! @} (documentation module close)*/


#endif /*D_AXISHELL*/


