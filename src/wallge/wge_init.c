/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wgeinit' which initializes the element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALLGE
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  initialization routine for the gradient enhanced wall element

<pre>                                                              ah 05/03
This routine inits the integration points and allocates stresses.
</pre>

\param *actpart      PARTITION   (I)   my partition
\param *mat          MATERIAL    (I)   the element's material

\warning There is nothing special to this routine.
\return void

*----------------------------------------------------------------------*/
void wgeinit(PARTITION *actpart,MATERIAL *mat)
{
INT          i,j,k;
INT          size_i, size_j;
ELEMENT     *actele;
WALLGE_DATA    data;

ARRAY    funct_a_h;  /* shape functions */
DOUBLE  *funct_h;
ARRAY    deriv_a_h;  /* derivatives of shape functions */
DOUBLE **deriv_h;
ARRAY    xjm_a_h;    /* jacobian matrix */
DOUBLE **xjm_h;


/*----------------------------------------------------------------------*/
funct_h     = amdef("funct_h"  ,&funct_a_h,MAXNOD_WALL1,1 ,"DV");
deriv_h     = amdef("deriv_h"  ,&deriv_a_h,2,MAXNOD_WALL1 ,"DA");
xjm_h       = amdef("xjm_h"    ,&xjm_a_h  ,2,2            ,"DA");
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numele; i++)
{/*matdam00*/
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_wallge) continue;

  /*---------------------------------------- init integration points ---*/
  wgeintg(actele,&data,0);
  /*-------------------------------- allocate the space for stresses ---*/
  am4def("stress_GP",&(actele->e.wallge->stress_GP),1,7,25,0,"D3");
  am4def("stress_ND",&(actele->e.wallge->stress_ND),1,7,MAXNOD,0,"D3");
/*----------------------------------------------------------------------*/

  /*--------------------------------------------- init working array ---*/
  /*--------------------------------------------- damage ---*/
  if(mat[actele->mat-1].mattyp == m_damage_ge )
  {/*matdam01*/
    size_i = 1;
    actele->e.wallge->elwa = (WGE_ELE_WA*)CCACALLOC(size_i,sizeof(WGE_ELE_WA));
    if (actele->e.wallge->elwa==NULL)
    {
      dserror("Allocation of elwa in ELEMENT failed");
      break;
    }
    size_j = actele->e.wallge->nGP[0] * actele->e.wallge->nGP[1];
    actele->e.wallge->elwa[0].iptwa =
                               (WGE_IP_WA*)CCACALLOC(size_j,sizeof(WGE_IP_WA));
    if (actele->e.wallge->elwa[0].iptwa==NULL)
    {
      dserror("Allocation of iptwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {
      actele->e.wallge->elwa[0].iptwa[k].yip    = -1;
      actele->e.wallge->elwa[0].iptwa[k].kappa  = 0.0;
      actele->e.wallge->elwa[0].iptwa[k].damage = 0.0;
      actele->e.wallge->elwa[0].iptwa[k].aequistrain = 0.0;
      actele->e.wallge->elwa[0].iptwa[k].aequistrain_nl = 0.0;
      for (j=0; j<4; j++)
      {
        actele->e.wallge->elwa[0].iptwa[k].sig[j] = 0.0;
      }
    }
  }/*matdam01*/
  /*-----------------------------------------------------------*/
   /*-------------------------------------------------------------------*/
}/*matdam00*/
/*----------------------------------------------------------------------*/
amdel(&funct_a_h);
amdel(&deriv_a_h);
amdel(&xjm_a_h  );
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_enter("wgeinit");
#endif


return;
} /* end of wgeinit */
/*----------------------------------------------------------------------*/
#endif /*D_WALLGE*/
/*! @} (documentation module close)*/
#endif
