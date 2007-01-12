/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_static_ke' which integrates the
linear stiffness for the 3d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element

<pre>                                                              mn 06/02
This routine integrates the linear stiffness for the 3d ale element


</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE3_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void
\sa calling: ale3_intg(), ale3_funct_deriv(), ale3_jaco(), ale3_bop(),
             ale3_mat_linel(), ale3_keku(), ale3_hourglass();
             called by: ale3()

*----------------------------------------------------------------------*/
void ale3_static_ke(
    ELEMENT    *ele,
    ALE3_DATA  *data,
    MATERIAL   *mat,
    ARRAY      *estif_global,
    INT         init
    )
{
  INT                 nir=0,nis=0,nit=0;    /* num GP in r/s/t direction */
  INT                 lr, ls, lt;           /* loopers over GP */
  INT                 iel;                  /* numnp to this element */
  INT                 nd;

  DOUBLE              fac;
  DOUBLE              e1=0.0,e2=0.0,e3=0.0; /*GP-coords*/
  DOUBLE              facr=0.0,facs=0.0,fact=0.0;   /* weights at GP */
  DOUBLE              vol;                  /* element volume */

  static DOUBLE   D[6][6];
  static DOUBLE   funct[MAXNOD_ALE3];
  static DOUBLE   deriv[3][MAXNOD_ALE3];

  static DOUBLE   xjm[3][3];
  static DOUBLE   bop[6][3*MAXNOD_ALE3];

  static DOUBLE **estif;    /* element stiffness matrix ke */

  DOUBLE det;


#ifdef DEBUG
  dstrc_enter("ale3_static_ke");
#endif

  if (init == 1)
    goto end;

  /* integration parameters */
  ale3_intg(ele,data);

  /* some of the fields have to be reinitialized to zero */
  amzero(estif_global);
  estif     = estif_global->a.da;

  /* integration parameters */
  switch(ele->distyp)
  {
    case hex8:
    case hex20:
    case hex27:
      nir     = ele->e.ale3->nGP[0];
      nis     = ele->e.ale3->nGP[1];
      nit     = ele->e.ale3->nGP[2];
      break;
    case tet4:
    case tet10:
      nir = ele->e.ale3->nGP[0];
      nis = 1;
      nit = 1;
      break;
    default:
      dserror("unknown number of gaussian points in ale2_intg");
      break;
  }
  iel     = ele->numnp;
  nd      = 3 * iel;
  vol     = 0.0;


  /* integration loops */
  for (lr=0; lr<nir; lr++)
  {
    /* gaussian point and weight at it */
    e1   = data->xgpr[lr];
    facr = data->wgtr[lr];
    for (ls=0; ls<nis; ls++)
    {
      /* gaussian point and weight at it */
      for (lt=0; lt<nit; lt++)
      {
        /* gaussian point and weight at it */
        switch(ele->distyp)
        {
          case hex8:
          case hex20:
          case hex27:
            e2   = data->xgps[ls];
            facs = data->wgts[ls];
            e3   = data->xgpt[lt];
            fact = data->wgtt[lt];
            break;
          case tet4:
          case tet10:
            e2   = data->xgps[lr];
            facs = ONE;
            e3   = data->xgpt[lr];
            fact = ONE;
            break;
          default:
            dserror("unknown number of gaussian points in ale2_intg");
            break;
        }

        /* shape functions and their derivatives */
        ale3_funct_deriv(funct,deriv,e1,e2,e3,ele->distyp,1);

        /* compute jacobian matrix */
        ale3_jaco(deriv,xjm,&det,ele,iel);

        /* calculate element volume */
        vol = vol + facr*facs*fact*det;

        /* use jacobian determinant or not */
        if(ele->e.ale3->jacobi==1)
          fac = facr * facs *  fact * det;
        else
          fac = facr * facs *  fact;

        /* calculate operator B */
        ale3_bop(bop,deriv,xjm,det,iel);

        /* call material law */
        ale3_mat_linel(mat->m.stvenant,D);

        /* elastic stiffness matrix ke */
        ale3_keku(estif,bop,D,fac,nd);

        /* hourglass stabalization  stiffness matrix ke */
        if(ele->distyp==hex8 && nir == 1 && nis == 1 && nit ==1)
          ale3_hourglass(ele,estif,vol);

      }  /* end of loop over lt */
    }  /* end of loop over ls */
  }  /* end of loop over lr */


  /* local co-system */
  dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_static_ke */




/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element

<pre>                                                              mn 06/02
This routine integrates the linear stiffness for the 3d ale element


</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE3_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void
\sa calling: ale3_intg(), ale3_funct_deriv(), ale3_jaco(), ale3_bop(),
             ale3_mat_linel(), ale3_keku(), ale3_hourglass();
             called by: ale3()

*----------------------------------------------------------------------*/
void ale3_static_ke_red(
    ELEMENT    *ele,
    ALE3_DATA  *data,
    MATERIAL   *mat,
    ARRAY      *estif_global,
    INT         init
    )
{

  INT             nir=0,nis=0,nit=0;    /* num GP in r/s/t direction */
  INT             lr, ls, lt;           /* loopers over GP */
  INT             iel=0;                /* numnp to this element */
  INT             nd;

  DOUBLE          fac;
  DOUBLE          e1=0.0,e2=0.0,e3=0.0; /*GP-coords*/
  DOUBLE          facr=0.0,facs=0.0,fact=0.0;   /* weights at GP */
  DOUBLE          vol;                  /* element volume */
  DOUBLE          det;

  static DOUBLE   D[6][6];
  static DOUBLE   funct[MAXNOD_ALE3];
  static DOUBLE   deriv[3][MAXNOD_ALE3];

  static DOUBLE   xjm[3][3];
  static DOUBLE   bop[6][3*MAXNOD_ALE3];

  static DOUBLE **estif;    /* element stiffness matrix ke */



#ifdef DEBUG
  dstrc_enter("ale3_static_ke_red");
#endif

  if (init == 1)
    goto end;


  /* integration parameters */
  /* identical for hex8 and hex20 */
  ale3_intg(ele,data);


  /* some of the fields have to be reinitialized to zero */
  amzero(estif_global);
  estif     = estif_global->a.da;


  /* integration parameters */
  switch(ele->distyp)
  {
    case hex8:
    case hex20:
    case hex27:
      nir     = ele->e.ale3->nGP[0];
      nis     = ele->e.ale3->nGP[1];
      nit     = ele->e.ale3->nGP[2];
      /* use only 8 nodes */
      iel     = 8;
      break;

    case tet4:
    case tet10:
      nir = ele->e.ale3->nGP[0];
      nis = 1;
      nit = 1;
      /* use only 4 nodes */
      iel     = 4;
      break;

    default:
      dserror("unknown number of gaussian points in ale2_intg");
      break;
  }

  nd      = 3 * iel;
  vol     = 0.0;


  /* integration loops */
  for (lr=0; lr<nir; lr++)
  {
    /* gaussian point and weight at it */
    e1   = data->xgpr[lr];
    facr = data->wgtr[lr];
    for (ls=0; ls<nis; ls++)
    {
      /* gaussian point and weight at it */
      for (lt=0; lt<nit; lt++)
      {
        /* gaussian point and weight at it */
        switch(ele->distyp)
        {
          case hex8:
          case hex20:
          case hex27:
            e2   = data->xgps[ls];
            facs = data->wgts[ls];
            e3   = data->xgpt[lt];
            fact = data->wgtt[lt];
            break;
          case tet4:
          case tet10:
            e2   = data->xgps[lr];
            facs = ONE;
            e3   = data->xgpt[lr];
            fact = ONE;
            break;
          default:
            dserror("unknown number of gaussian points in ale2_intg");
            break;
        }

        /* shape functions and their derivatives */
        ale3_funct_deriv(funct,deriv,e1,e2,e3,hex8,1);

        /* compute jacobian matrix */
        ale3_jaco (deriv,xjm,&det,ele,iel);

        /* calculate element volume */
        vol = vol + facr*facs*fact*det;

        /* use jacobian determinant or not */
        if(ele->e.ale3->jacobi==1)
          fac = facr * facs *  fact * det;
        else
          fac = facr * facs *  fact;

        /* calculate operator B */
        ale3_bop(bop,deriv,xjm,det,iel);

        /* call material law */
        ale3_mat_linel(mat->m.stvenant,D);

        /* elastic stiffness matrix ke */
        ale3_keku(estif,bop,D,fac,nd);

        /* hourglass stabalization  stiffness matrix ke */
        if(ele->distyp==hex8 && nir == 1 && nis == 1 && nit ==1)
          ale3_hourglass(ele,estif,vol);

      }  /* end of loop over lt */
    }  /* end of loop over ls */
  }  /* end of loop over lr */




  /* node 8                  node 0                node 1      */
  estif[24][24] = 1.0; estif[24][ 0] = -0.5; estif[24][ 3] = -0.5;
  estif[25][25] = 1.0; estif[25][ 1] = -0.5; estif[25][ 4] = -0.5;
  estif[26][26] = 1.0; estif[26][ 2] = -0.5; estif[26][ 5] = -0.5;

  /* node 9                  node 1                node 2       */
  estif[27][27] = 1.0; estif[27][ 3] = -0.5; estif[27][ 6] = -0.5;
  estif[28][28] = 1.0; estif[28][ 4] = -0.5; estif[28][ 7] = -0.5;
  estif[29][29] = 1.0; estif[29][ 5] = -0.5; estif[29][ 8] = -0.5;

  /* node 10                 node 2                node 3       */
  estif[30][30] = 1.0; estif[30][ 6] = -0.5; estif[30][ 9] = -0.5;
  estif[31][31] = 1.0; estif[31][ 7] = -0.5; estif[31][10] = -0.5;
  estif[32][32] = 1.0; estif[32][ 8] = -0.5; estif[32][11] = -0.5;

  /* node 11                 node 3                node 0       */
  estif[33][33] = 1.0; estif[33][ 9] = -0.5; estif[33][ 0] = -0.5;
  estif[34][34] = 1.0; estif[34][10] = -0.5; estif[34][ 1] = -0.5;
  estif[35][35] = 1.0; estif[35][11] = -0.5; estif[35][ 2] = -0.5;

  /* node 12                 node 0                node 4       */
  estif[36][36] = 1.0; estif[36][ 0] = -0.5; estif[36][12] = -0.5;
  estif[37][37] = 1.0; estif[37][ 1] = -0.5; estif[37][13] = -0.5;
  estif[38][38] = 1.0; estif[38][ 2] = -0.5; estif[38][14] = -0.5;

  /* node 13                 node 1                node 5       */
  estif[39][39] = 1.0; estif[39][ 3] = -0.5; estif[39][15] = -0.5;
  estif[40][40] = 1.0; estif[40][ 4] = -0.5; estif[40][16] = -0.5;
  estif[41][41] = 1.0; estif[41][ 5] = -0.5; estif[41][17] = -0.5;

  /* node 14                 node 2                node 6       */
  estif[42][42] = 1.0; estif[42][ 6] = -0.5; estif[42][18] = -0.5;
  estif[43][43] = 1.0; estif[43][ 7] = -0.5; estif[43][19] = -0.5;
  estif[44][44] = 1.0; estif[44][ 8] = -0.5; estif[44][20] = -0.5;

  /* node 15                 node 3                node 7       */
  estif[45][45] = 1.0; estif[45][ 9] = -0.5; estif[45][21] = -0.5;
  estif[46][46] = 1.0; estif[46][10] = -0.5; estif[46][22] = -0.5;
  estif[47][47] = 1.0; estif[47][11] = -0.5; estif[47][23] = -0.5;

  /* node 16                 node 4                node 5       */
  estif[48][48] = 1.0; estif[48][12] = -0.5; estif[48][15] = -0.5;
  estif[49][49] = 1.0; estif[49][13] = -0.5; estif[49][16] = -0.5;
  estif[50][50] = 1.0; estif[50][14] = -0.5; estif[50][17] = -0.5;

  /* node 17                 node 5                node 6       */
  estif[51][51] = 1.0; estif[51][15] = -0.5; estif[51][18] = -0.5;
  estif[52][52] = 1.0; estif[52][16] = -0.5; estif[52][19] = -0.5;
  estif[53][53] = 1.0; estif[53][17] = -0.5; estif[53][20] = -0.5;

  /* node 18                 node 6                node 7       */
  estif[54][54] = 1.0; estif[54][18] = -0.5; estif[54][21] = -0.5;
  estif[55][55] = 1.0; estif[55][19] = -0.5; estif[55][22] = -0.5;
  estif[56][56] = 1.0; estif[56][20] = -0.5; estif[56][23] = -0.5;

  /* node 19                 node 7                node 4       */
  estif[57][57] = 1.0; estif[57][21] = -0.5; estif[57][12] = -0.5;
  estif[58][58] = 1.0; estif[58][22] = -0.5; estif[58][13] = -0.5;
  estif[59][59] = 1.0; estif[59][23] = -0.5; estif[59][14] = -0.5;




  /* local co-system */
  dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");


end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_static_ke_red */


#endif
/*! @} (documentation module close)*/
