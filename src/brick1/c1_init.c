/*!----------------------------------------------------------------------
\file
\brief initialize the 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief initialize the 3D hex element

<pre>                                                              al 06/02
This routine initializes the 3D-hex-element.

</pre>
\param *actpart      PARTITION   (i)   my partition
\param *mat            MATERIAL  (i)   my material

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: brick1()

*----------------------------------------------------------------------*/
void c1init(PARTITION *actpart,MATERIAL    *mat )
{
INT          i,j,k,ngp,nnp;
INT          size_i, size_j;
ELEMENT     *actele;

#ifdef DEBUG
dstrc_enter("c1init");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numele; i++)
{
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_brick1) continue;
  /*--------------------------------------------------------------------*/
  /* number of nodal points */
  nnp = actele->numnp;
  /* number of gaussian points */
  ngp = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
  /*-------------------------------- allocate the space for stresses ---*/
  if(actele->e.c1->stresstyp!=c1_nostr)
  {
    am4def("stress_GP",&(actele->e.c1->stress_GP),1,27,ngp,0,"D3");
    am4def("stress_ND",&(actele->e.c1->stress_ND),1,27,nnp,0,"D3");
  }
  /*--------------------------------------------- init working array ---*/
  if(mat[actele->mat-1].mattyp == m_pl_mises    ||
     mat[actele->mat-1].mattyp == m_pl_mises_ls ||
     mat[actele->mat-1].mattyp == m_pl_foam     ||
     mat[actele->mat-1].mattyp == m_nhmfcc      ||
     mat[actele->mat-1].mattyp == m_stvenpor    ||
     mat[actele->mat-1].mattyp == m_mfoc        ||
     mat[actele->mat-1].mattyp == m_mfcc        ||
     mat[actele->mat-1].mattyp == m_pl_hash     ||
     actele->e.c1->nhyb>0)
  {
    size_i = 1;
    actele->e.c1->elewa = (C1_ELE_WA*)CCACALLOC(size_i,sizeof(C1_ELE_WA));
    if (actele->e.c1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    }
  }
  /*------------------------- init working array for porous material ---*/
  if(mat[actele->mat-1].mattyp == m_stvenpor)
  {
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->matdata[0] = current density value  |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->matdata = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    if (actele->e.c1->elewa->matdata==NULL)
    {
      dserror("Allocation of matdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->matdata[0] = mat[actele->mat-1].m.stvenpor->density;
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->optdata[0] = current opt.var.num.   |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->optdata = (INT*)CCACALLOC(size_j,sizeof(INT));
    if (actele->e.c1->elewa->optdata==NULL)
    {
      dserror("Allocation of optdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->optdata[0] = 0;
  }
  /*------------------------- init working array for foam materials ---*/
  if(mat[actele->mat-1].mattyp == m_mfoc)
  {
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->matdata[0] = current density value  |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->matdata = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    if (actele->e.c1->elewa->matdata==NULL)
    {
      dserror("Allocation of matdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->matdata[0] = mat[actele->mat-1].m.mfoc->dens;
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->optdata[0] = current opt.var.num.   |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->optdata = (INT*)CCACALLOC(size_j,sizeof(INT));
    if (actele->e.c1->elewa->optdata==NULL)
    {
      dserror("Allocation of optdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->optdata[0] = 0;
  }
  if(mat[actele->mat-1].mattyp == m_mfcc)
  {
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->matdata[0] = current density value  |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->matdata = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    if (actele->e.c1->elewa->matdata==NULL)
    {
      dserror("Allocation of matdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->matdata[0] = mat[actele->mat-1].m.mfcc->dens;
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->optdata[0] = current opt.var.num.   |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->optdata = (INT*)CCACALLOC(size_j,sizeof(INT));
    if (actele->e.c1->elewa->optdata==NULL)
    {
      dserror("Allocation of optdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->optdata[0] = 0;
  }
  if(mat[actele->mat-1].mattyp == m_nhmfcc)
  {
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->matdata[0] = current density value  |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->matdata = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    if (actele->e.c1->elewa->matdata==NULL)
    {
      dserror("Allocation of matdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->matdata[0] = mat[actele->mat-1].m.nhmfcc->dens;
    /*----------------------------------------------------------*
     | actele->e.c1->elewa->optdata[0] = current opt.var.num.   |
     *----------------------------------------------------------*/
    size_j = 1;
    actele->e.c1->elewa->optdata = (INT*)CCACALLOC(size_j,sizeof(INT));
    if (actele->e.c1->elewa->optdata==NULL)
    {
      dserror("Allocation of optdata in ELEMENT failed");
      break;
    }
    actele->e.c1->elewa->optdata[0] = 0;
  }
  /*------------------------------------- init working array for eas ---*/
  if(actele->e.c1->nhyb>0)
  {
    /* */
    size_i = 1;
    actele->e.c1->elewa[0].eas =
                             (C1_EASDAT*)CCACALLOC(size_i,sizeof(C1_EASDAT));
    /* */
    size_j = actele->e.c1->nhyb*3;
    actele->e.c1->elewa[0].eas[0].ehdis
                                  = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    for (k=0; k<size_j; k++) actele->e.c1->elewa[0].eas[0].ehdis[k]=0.;
    /* */
    size_j = 24; /* 8noded element only!*/
    actele->e.c1->elewa[0].eas[0].disl
                                  = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    for (k=0; k<size_j; k++) actele->e.c1->elewa[0].eas[0].disl[k]=0.;
    /* */
    size_j = 24*actele->e.c1->nhyb*3;
    actele->e.c1->elewa[0].eas[0].hil
                                  = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    for (k=0; k<size_j; k++) actele->e.c1->elewa[0].eas[0].hil[k]=0.;
    /* */
    size_j = actele->e.c1->nhyb*3;
    actele->e.c1->elewa[0].eas[0].hih
                                  = (DOUBLE*)CCACALLOC(size_j,sizeof(DOUBLE));
    for (k=0; k<size_j; k++) actele->e.c1->elewa[0].eas[0].hih[k]=0.;
    /* */
  }
  /*------------------------------ init working array for plasticity ---*/
  if(mat[actele->mat-1].mattyp == m_pl_mises ||
     mat[actele->mat-1].mattyp == m_pl_foam  ||
     mat[actele->mat-1].mattyp == m_pl_hash)
  {
    size_j = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
    actele->e.c1->elewa[0].ipwa =
                               (C1_IP_WA*)CCACALLOC(size_j,sizeof(C1_IP_WA));
    if (actele->e.c1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {
      actele->e.c1->elewa[0].ipwa[k].kappa = 0.;
      actele->e.c1->elewa[0].ipwa[k].epstn = 0.;
      actele->e.c1->elewa[0].ipwa[k].yip   = -1;
      actele->e.c1->elewa[0].ipwa[k].imod  = 0;

      for (j=0; j<6; j++)
      {
        actele->e.c1->elewa[0].ipwa[k].sig[j] = 0.;
        actele->e.c1->elewa[0].ipwa[k].eps[j] = 0.;
      }
    }
   /*-------------------------------------------------------------------*/
  }
  /*------- init working array for mises plasticity - finite strains ---*/
  if(mat[actele->mat-1].mattyp == m_pl_mises_ls)
  {
    size_j = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
    actele->e.c1->elewa[0].ipwa =
                               (C1_IP_WA*)CCACALLOC(size_j,sizeof(C1_IP_WA));
    if (actele->e.c1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {
      actele->e.c1->elewa[0].ipwa[k].epstn = 0.;
      actele->e.c1->elewa[0].ipwa[k].yip   = -1;

      for (j=0; j<6; j++) actele->e.c1->elewa[0].ipwa[k].sig[j] = 0.;
      for (j=0; j<9; j++) actele->e.c1->elewa[0].ipwa[k].eps[j] = 0.;
      for (j=0; j<3; j++) actele->e.c1->elewa[0].ipwa[k].sig[j] = 1.;
      for (j=0; j<3; j++) actele->e.c1->elewa[0].ipwa[k].eps[j] = 1.;
    }
  }
  /*------- init working array for foam plasticity - finite strains ---*/
  if(mat[actele->mat-1].mattyp == m_pl_foam)
  {
    size_j = actele->e.c1->nGP[0]*actele->e.c1->nGP[1]*actele->e.c1->nGP[2];
    actele->e.c1->elewa[0].ipwa =
                               (C1_IP_WA*)CCACALLOC(size_j,sizeof(C1_IP_WA));

    if (actele->e.c1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    }
    for (k=0; k<size_j; k++)
    {
      actele->e.c1->elewa[0].ipwa[k].epstn = 0.;
      actele->e.c1->elewa[0].ipwa[k].yip   = -1;

      for (j=0; j<6; j++) actele->e.c1->elewa[0].ipwa[k].sig[j] = 0.;
      for (j=0; j<9; j++) actele->e.c1->elewa[0].ipwa[k].eps[j] = 0.;
      for (j=0; j<3; j++) actele->e.c1->elewa[0].ipwa[k].sig[j] = 1.;
      for (j=0; j<3; j++) actele->e.c1->elewa[0].ipwa[k].eps[j] = 1.;
    }
   /*-------------------------------------------------------------------*/
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1init */
#endif
/*! @} (documentation module close)*/

