/*!----------------------------------------------------------------------
\file
\brief contains the routine 'saxistatic_ke' which forms the linear stiffness
ke for axisymmetric shell element and the routine 'saxi_statcond' which
condenses the stiffness matrix

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_AXISHELL

#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*!
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  Computation of linear stiffness ke for axialsymmetric shell element

<pre>                                                              mn 05/03
This routine computes the linear stiffness for the axialsymmetric shell
element by computing first a (7,7)-stiffness matrix. Then this matrix is
condensed from (7,7) to (6,6) by static condensation.


</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          SAXI_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stiffness matrix (6,6)
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void
\sa calling:   saxiintg, saxi_B, saxi_mat, saxi_keku, saxi_statcond;
    called by: axishell();

*----------------------------------------------------------------------*/
void saxistatic_ke(
    ELEMENT   *ele,
    SAXI_DATA *data,
    MATERIAL  *mat,
    ARRAY     *estif_global,
    INT        init
    )
{

  INT                 i,j;                /* some loopers                 */
  INT                 nixsi;              /* num GP in xsi direction      */
  INT                 lxsi;               /* looper over GP               */
  INT                 iel;                /* numnp to this element        */
  INT                 nd;
  INT                 istore = 0;/* controls storing of new stresses to wa*/

  DOUBLE              thick;              /* thickness                    */
  DOUBLE              xsi;                /* GP-coord                     */
  DOUBLE              r,dr,dz,dl,cosa,sina; /* infos about the element    */
  DOUBLE              fac;                /* weights at GP                */

  static DOUBLE   D[4][4];
  static DOUBLE   work1[4][7];
  static DOUBLE   work2[7][7];
  static DOUBLE   B[4][7];
  static DOUBLE   estif_7[7][7];          /* ke(7,7) before static condensation */

  static DOUBLE **estif;            /* ke(6,6) after  static condensation */

#ifdef DEBUG
  dstrc_enter("saxistatic_ke");
#endif

  istore = 0;

  /* init phase (init==1) */
  if (init==1)
    goto end;

  /* uninit phase (init==-1) */
  else if (init==-1)
    goto end;

  else if(init==2)
  {
    istore = 1;
  }

  /* integration parameters */
  saxiintg(data);

  /* some of the fields have to be reinitialized to zero */
  amzero(estif_global);
  estif     = estif_global->a.da;

  for (i=0; i<7; i++)
  {
    for (j=0; j<7; j++)
    {
      estif_7[i][j] = 0.0;
    }
  }
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

  /* integration loops */
  for (lxsi=0; lxsi<nixsi; lxsi++)
  {
    /* gaussian point, weight and thickness at it */
    xsi   = data->xgr[lxsi];                 /* xsi       */
    fac   = data->wgt[lxsi];                 /* weight    */
    r     = ele->node[0]->x[0] + cosa*xsi*dl;/* radius    */
    thick = ele->e.saxi->thick[0] +          /* thickness */
      xsi/1.0*(ele->e.saxi->thick[1]-ele->e.saxi->thick[0]);

    /* calculate operator B */
    saxi_B(B,xsi,r,dl,cosa,sina);

    /* compute constitutive matrix D */
    saxi_mat(mat->m.stvenant,D,thick);

    if(istore==0)
    {
      /* element stiffness matrix ke */
      saxi_keku(estif_7,B,D,work1,work2,fac,r,dl);
    }

  }/* end of loop over lxsi */

  /* static condensation to reduce the stiffness matrix
     from (7,7) to (6,6) by elemination of the middle node */
  saxi_statcond(ele,estif,estif_7,cosa,sina);

/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of saxistatic_ke */



/*!----------------------------------------------------------------------
  \brief  Static condensation to reduce the stiffnes matrix

  <pre>                                                              mn 05/03
  This routine reduces the linear stiffness for the axialsymmetric shell from
  a (7,7)-stiffness matrix to a (6,6)-stiffness matrix by static condensation.


  </pre>
  \param  *ele      ELEMENT (i)   my element
  \param **estif    DOUBLE  (o)   the stiffness matrix (6,6)
  \param   estif_7  DOUBLE[][]  (i)   the stiffness matrix (7,7)
  \param   cosa     DOUBLE  (i)   the cosine of the angle of the current element
  \param   sina     DOUBLE  (i)   the sine of the angle of the current element

  \warning There is nothing special to this routine
  \return void
  \sa calling:   ---;
  called by: saxistatic_ke();

 *----------------------------------------------------------------------*/
void saxi_statcond(
    ELEMENT *ele,
    DOUBLE **estif,
    DOUBLE   estif_7[7][7],
    DOUBLE    cosa,
    DOUBLE    sina
    )
{

  INT i,j,k;
  DOUBLE k77;

  DOUBLE estif_local[6][6];   /* local stiffness matrix (6,6) */
  DOUBLE T[6][6];             /* Transformation matrix        */
  DOUBLE k71[6];
  DOUBLE k17[6];
  DOUBLE work[6][6];
  DOUBLE *statcond;          /* pointer to statcond vector of the element */

#ifdef DEBUG
  dstrc_enter("saxi_statcond");
#endif


  /* compute the transformation matrix T */
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      T[i][j] = 0.0;
    }
  }

  if (ele->node[0]->gnode->ondesigntyp==ondnode && ele->node[0]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    T[0][0] =  1.0;
    T[1][1] =  1.0;
    T[2][2] =  1.0;
  }
  else
  {
    /* global cos */
    T[0][0] =  cosa;
    T[0][1] =  sina;
    T[1][0] =  sina;
    T[1][1] = -cosa;
    T[2][2] =  -1.0;
  }

  if (ele->node[1]->gnode->ondesigntyp==ondnode && ele->node[1]->gnode->d.dnode->cos_type == 1)
  {
    /* local cos */
    T[3][3] =  1.0;
    T[4][4] =  1.0;
    T[5][5] =  1.0;
  }
  else
  {
    /* global cos */
    T[3][3] =  cosa;
    T[3][4] =  sina;
    T[4][3] =  sina;
    T[4][4] = -cosa;
    T[5][5] =  -1.0;
  }

  /*
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
     */

  /* compute the local 6/6 stiffnes matrix */
  k77 = estif_7[6][6];
  ele->e.saxi->Sk = 1.0/k77;
  for (i=0; i<6; i++)
  {
    k71[i] = -(1.0/k77)*estif_7[6][i];
    k17[i] = estif_7[i][6];
  }

  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      estif_local[i][j] = estif_7[i][j] + k17[i] * k71[j];
    }
  }

  /* compute the final global 6/6 stiffnes matrix */
  /* make multiplication: work = esfif_local * T */
  for (k=0; k<6; k++)
  {
    for (i=0; i<6; i++)
    {
      work[i][k] = 0.0;
      for (j=0; j<6; j++)
      {
        work[i][k] += estif_local[i][j] * T[j][k];
      }
    }
  }

  /* make multiplication: estif = T^T * work */
  for (k=0; k<6; k++)
  {
    for (i=0; i<6; i++)
    {
      estif[i][k] = 0.0;
      for (j=0; j<6; j++)
      {
        estif[i][k] += T[j][i] * work[j][k];
      }
    }
  }

  /* compute statcond which is substracted from the rhs */
  statcond = ele->e.saxi->statcond;
  for (i=0; i<6; i++)
  {
    statcond[i] = -(1.0/k77)*k17[i]; /* local */
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}/* end of saxi_statcond */

/*! @} (documentation module close)*/
#endif /*D_AXISHELL*/
#endif
