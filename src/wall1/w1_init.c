/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1init' which initializes the element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | initialize the element                                    al 6/01    |
 *----------------------------------------------------------------------*/
void w1init(PARTITION *actpart,MATERIAL *mat)
{
INT          i,j,k,l,m,ncm;
INT          size_i, size_j;
ELEMENT     *actele;
W1_DATA      data;

ARRAY    funct_a_h;  /* shape functions */
DOUBLE  *funct_h;
ARRAY    deriv_a_h;  /* derivatives of shape functions */
DOUBLE **deriv_h;
ARRAY    xjm_a_h;    /* jacobian matrix */
DOUBLE **xjm_h;

#ifdef GEMM
DOUBLE ***b_bar_history; /* previous B_hat operator*/
DOUBLE ***PK_history; /* previous 2nd PK stresses*/
#endif

#ifdef DEBUG
dstrc_enter("w1init");
#endif
/*----------------------------------------------------------------------*/
funct_h     = amdef("funct_h"  ,&funct_a_h,MAXNOD_WALL1,1 ,"DV");
deriv_h     = amdef("deriv_h"  ,&deriv_a_h,2,MAXNOD_WALL1 ,"DA");
xjm_h       = amdef("xjm_h"    ,&xjm_a_h  ,2,2            ,"DA");
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numele; i++)
{/*matplast00*/
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_wall1) continue;
  /*---------------------------------------- init integration points ---*/
  w1intg(actele,&data,0);
/*---------------------------------------------------------fh 06/02-----*/
  /*-------------------------------- allocate the space for stresses ---*/
  am4def("stress_GP",&(actele->e.w1->stress_GP),1,7,MAXGAUSS,0,"D3");
  am4def("stress_ND",&(actele->e.w1->stress_ND),1,7,MAXNOD,0,"D3");
/*------------------------------------------init history for GEMM scheme*/
#ifdef GEMM
  b_bar_history = am4def("b_bar",&(actele->e.w1->b_bar_history),MAXGAUSS,4,MAXNOD_WALL1*2,0,"D3");
  PK_history    = am4def("PK",&(actele->e.w1->PK_history),MAXGAUSS,4,4,0,"D3");
  am4zero(&(actele->e.w1->b_bar_history));
  am4zero(&(actele->e.w1->PK_history));
#endif
/*---------------------------------------- init info for multiscale ---*/  
#ifdef D_MLSTRUCT
    actele->e.w1->isinomegaprime = 0;
    actele->e.w1->firstinomegaprime = 1;
    actele->e.w1->translation[0] = 0.0;
    actele->e.w1->translation[1] = 0.0;
#endif /* D_MLSTRUCT */
  /*--------------------------------------------- init working array ---*/
  if(mat[actele->mat-1].mattyp == m_stvenpor)
  {
    size_i = 1;
    actele->e.w1->elewa = (W1_ELE_WA*)CCACALLOC(size_i,sizeof(W1_ELE_WA));
    if (actele->e.w1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    }
    /*----------------------------------------------------------*
     | actele->e.w1->elewa->matdata[0] = current density value  |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.w1->elewa->matdata = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    if (actele->e.w1->elewa->matdata==NULL)
    {
      dserror("Allocation of matdata in ELEMENT failed");
      break;
    }
    actele->e.w1->elewa->matdata[0] = mat[actele->mat-1].m.stvenpor->density;
    /*----------------------------------------------------------*
     | actele->e.w1->elewa->optdata[0] = current opt.var.num.   |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.w1->elewa->optdata = (INT*)CCACALLOC(size_j,sizeof(INT));
    if (actele->e.w1->elewa->optdata==NULL)
    {
      dserror("Allocation of optdata in ELEMENT failed");
      break;
    }
    actele->e.w1->elewa->optdata[0] = 0;
  }
  /*--------------------------------------------- damage ---*/
  if(mat[actele->mat-1].mattyp == m_dam_mp )
  {/*matdam01*/
    size_i = 1;
    actele->e.w1->elewa = (W1_ELE_WA*)CCACALLOC(size_i,sizeof(W1_ELE_WA));
    if (actele->e.w1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    }
    size_j = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
    actele->e.w1->elewa[0].ipwa =
                               (W1_IP_WA*)CCACALLOC(size_j,sizeof(W1_IP_WA));
    if (actele->e.w1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {
      actele->e.w1->elewa[0].ipwa[k].yip    = -1;
      actele->e.w1->elewa[0].ipwa[k].kappa  = 0.0;
      actele->e.w1->elewa[0].ipwa[k].damage = 0.0;
      actele->e.w1->elewa[0].ipwa[k].aequistrain = 0.0;
      for (j=0; j<4; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].sig[j] = 0.0;
      }
    }
  }/*matdam01*/
  /*-----------------------------------------------------------*/
  /* for plasticity and 3D-damage */
  if(mat[actele->mat-1].mattyp == m_pl_mises ||
     mat[actele->mat-1].mattyp == m_pl_mises_3D ||  /*Stefan's mises 3D*/
     mat[actele->mat-1].mattyp == m_pl_dp ||
     mat[actele->mat-1].mattyp == m_pl_epc ||
     mat[actele->mat-1].mattyp == m_pl_epc3D ||
     mat[actele->mat-1].mattyp == m_damage )
  {/*matplast01*/
    size_i = 1;
    actele->e.w1->elewa = (W1_ELE_WA*)CCACALLOC(size_i,sizeof(W1_ELE_WA));
    if (actele->e.w1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    }

    size_j = actele->e.w1->nGP[0] * actele->e.w1->nGP[1];
    actele->e.w1->elewa[0].ipwa =
                               (W1_IP_WA*)CCACALLOC(size_j,sizeof(W1_IP_WA));
    if (actele->e.w1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {/*matplast02*/
      actele->e.w1->elewa[0].ipwa[k].epstn = 0.;
      actele->e.w1->elewa[0].ipwa[k].yip   = -1;
      actele->e.w1->elewa[0].ipwa[k].kap   = 0.;
      actele->e.w1->elewa[0].ipwa[k].dam   = 0.;
      actele->e.w1->elewa[0].ipwa[k].damage= 0.;
      actele->e.w1->elewa[0].ipwa[k].aequistrain = 0.;
      actele->e.w1->elewa[0].ipwa[k].qn = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));

      /*additional values needed for condensation       sh 08/02*/
      if(mat[actele->mat-1].mattyp == m_pl_mises_3D)
      {
      actele->e.w1->elewa[0].ipwa[k].sigc = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].sigi = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].epsi = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].di   = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      }

      if(mat[actele->mat-1].mattyp == m_pl_epc ||
         mat[actele->mat-1].mattyp == m_pl_epc3D)
      {
      actele->e.w1->elewa[0].ipwa[k].sigc = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].grad = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].dlam = (DOUBLE*)CCACALLOC(2,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].dlam[0] = 0.;
      actele->e.w1->elewa[0].ipwa[k].dlam[1] = 0.;
      actele->e.w1->elewa[0].ipwa[k].kappa_t = 0.;
      actele->e.w1->elewa[0].ipwa[k].kappa_c = 0.;
      actele->e.w1->elewa[0].ipwa[k].sigi = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].epsi = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].di   = (DOUBLE*)CCACALLOC(4,sizeof(DOUBLE));

      ncm = mat[actele->mat-1].m.pl_epc->maxreb;

      actele->e.w1->elewa[0].ipwa[k].rsig   = (DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].reps   = (DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].repstn = (DOUBLE*)CCACALLOC(ncm,sizeof(DOUBLE));
      actele->e.w1->elewa[0].ipwa[k].ryip   = (INT*)CCACALLOC(ncm,sizeof(INT));
      for (j=0; j<ncm; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].rsig[j]   = 0.;
        actele->e.w1->elewa[0].ipwa[k].reps[j]   = 0.;
        actele->e.w1->elewa[0].ipwa[k].repstn[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].ryip[j]   = -1;
      }
      }
      for (j=0; j<2; j++)
      {
       actele->e.w1->elewa[0].ipwa[k].t_d[j] = 0.;
      }
      for (j=0; j<4; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].sig[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].sig_esz[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].eps[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].eps_esz[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].d4_esz[j]  = 0.;
        actele->e.w1->elewa[0].ipwa[k].qn[ j] = 0.;
        if(mat[actele->mat-1].mattyp == m_pl_mises_3D  )
        {
        actele->e.w1->elewa[0].ipwa[k].sigc[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].sigi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].epsi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].di[  j] = 0.;
        }
        if(mat[actele->mat-1].mattyp == m_pl_epc ||
           mat[actele->mat-1].mattyp == m_pl_epc3D )
        {
        actele->e.w1->elewa[0].ipwa[k].sigc[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].grad[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].sigi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].epsi[j] = 0.;
        actele->e.w1->elewa[0].ipwa[k].di[  j] = 0.;
        }
      }
    }/*matplast02*/
  /*------------------------------------- calculate element diameter ---*/
    if( (mat[actele->mat-1].mattyp == m_pl_mises &&
        (fabs(0.0001 - mat[actele->mat-1].m.pl_mises->GF) > 0.0001) )
      || (mat[actele->mat-1].mattyp == m_pl_mises_3D &&
         (fabs(0.0001 - mat[actele->mat-1].m.pl_mises->GF) > 0.0001) ) )
    {
       w1cdia(actele, &data, funct_h, deriv_h, xjm_h);
    }
    else if(mat[actele->mat-1].mattyp == m_pl_epc ||
            mat[actele->mat-1].mattyp == m_pl_epc3D )
    {
       w1cdia(actele, &data, funct_h, deriv_h, xjm_h);
    }
  /*--------------------------------------------------------------------*/
  }/*matplast01*/
  /*-------------------------------------------------------------------*/
    if(actele->e.w1->modeltype == incomp_mode)
    {
     if(mat->mattyp == m_stvenant)
     {
      size_i = 1;
      actele->e.w1->elewa = (W1_ELE_WA*)CCACALLOC(size_i,sizeof(W1_ELE_WA));
      if (actele->e.w1->elewa==NULL)
      {
       dserror("Allocation of elewa in ELEMENT failed");
       break;
      }
     }
     size_i = 1;
     actele->e.w1->elewa[0].imodewa = (W1_IMODE_WA*)CCACALLOC(size_i,sizeof(W1_IMODE_WA));
     if (actele->e.w1->elewa[0].imodewa==NULL)
     {
      dserror("Allocation of imodewa in ELEMENT failed");
      break;
     }
     for(l=0;l<4;l++)
     {
       if(mat->mattyp != m_stvenant)
       {
         actele->e.w1->elewa[0].imodewa[0].fintn[l] = 0.0;
         actele->e.w1->elewa[0].imodewa[0].alpha[l] = 0.0;
       }
       for(m=0;m<4;m++)
       {
        actele->e.w1->elewa[0].imodewa[0].knninv[l][m]  = 0.0;
        actele->e.w1->elewa[0].imodewa[0].knc[2*l][m]   = 0.0;
        actele->e.w1->elewa[0].imodewa[0].knc[2*l+1][m] = 0.0;
       }
     }
    }
   /*-------------------------------------------------------------------*/
}/*matplast00*/
/*----------------------------------------------------------------------*/
amdel(&funct_a_h);
amdel(&deriv_a_h);
amdel(&xjm_a_h  );
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1init */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
