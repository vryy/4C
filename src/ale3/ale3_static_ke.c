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
  static DOUBLE   funct[MAXNOD_BRICK1];
  static DOUBLE   deriv[3][MAXNOD_BRICK1];

  static DOUBLE   xjm[3][3];
  static DOUBLE   bop[6][6*MAXNOD_BRICK1];

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
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
end:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_static_ke */

#endif
/*! @} (documentation module close)*/
