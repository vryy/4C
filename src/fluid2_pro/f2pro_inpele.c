/*!----------------------------------------------------------------------
\file
\brief read fluid2_pro element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO
#include "../headers/standardtypes.h"
#include "fluid2pro_prototypes.h"
#include "fluid2pro.h"
/*!---------------------------------------------------------------------
\brief read fluid2 element from input-file

<pre>                                                         genk 03/02
</pre>
\param  *ele	   ELEMENT	   (o)	   actual element
\return void

------------------------------------------------------------------------*/
void f2pro_inp(ELEMENT *ele)
{
  INT        i;             /* simply a counter                           */
  INT        ierr=0;        /* error flag                                 */
  char       buffer[50];

#ifdef DEBUG
  dstrc_enter("f2pro_inp");
#endif

  ele->e.f2pro = (FLUID2_PRO*)CCACALLOC(1,sizeof(FLUID2_PRO));

  /* read the element nodes */
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    ele->numnp=4;
    ele->distyp=quad4;
    ele->e.f2pro->ntyp=1;
    ele->e.f2pro->dm=dm_q1p0;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
    frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");

    /*
     * Lets not allocate too small pieces of memory. Two pointers
     * share one array peacefully. :) */
    ele->e.f2pro->press = (DOUBLE*)CCACALLOC(2, sizeof(DOUBLE));
    ele->e.f2pro->phi   = &(ele->e.f2pro->press[1]);
    ele->e.f2pro->dof   = (INT*)   CCACALLOC(2, sizeof(INT));
    ele->e.f2pro->ldof  = &(ele->e.f2pro->dof[1]);
  }
  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    ele->numnp=9;
    ele->distyp=quad9;
    ele->e.f2pro->ntyp=1;
    ele->e.f2pro->dm=dm_q2pm1;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
    frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");

    /*
     * Lets not allocate too small pieces of memory. Two pointers
     * share one array peacefully. :) */
    ele->e.f2pro->press = (DOUBLE*)CCACALLOC(6, sizeof(DOUBLE));
    ele->e.f2pro->phi   = &(ele->e.f2pro->press[3]);
    ele->e.f2pro->dof   = (INT*)   CCACALLOC(6, sizeof(INT));
    ele->e.f2pro->ldof  = &(ele->e.f2pro->dof[3]);
  }
  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    dserror("quad8 for FLUID2_PRO not implemented yet!\n");
/*    ele->numnp=8;
      ele->distyp=quad8;
      ele->e.f2->ntyp=1;
      ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
      if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
      frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
      if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");*/
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    dserror("tri3 for FLUID2_PRO not implemented yet!\n");
/*   ele->numnp=3;
     ele->distyp=tri3;
     ele->e.f2->ntyp=2;
     ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
     if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
     frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
     if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");*/
  }
  frchk("TRI6",&ierr);
  if (ierr==1)
  {
    dserror("tri6 for FLUID2_PRO not implemented yet!\n");
/*   ele->numnp=6;
     ele->distyp=tri6;
     ele->e.f2->ntyp=2;
     ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
     if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
     frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
     if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");*/
  }
/*------------------------------------------ reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1) dserror("Reading of FLUID2_PRO element failed\n");
  if (ele->mat==0) dserror("No material defined for FLUID2_PRO element\n");
/*-------------------------------------------- read the gaussian points */
  if (ele->numnp==4 || ele->numnp==8 || ele->numnp==9)
  {
    frint_n("GP",&(ele->e.f2pro->nGP[0]),2,&ierr);
    if (ierr!=1) dserror("Reading of FLUID2_PRO element failed: integration\n");
  }
/*-------------------------- read gaussian points for triangle elements */
  if (ele->numnp==3 || ele->numnp==6)
  {
    dserror("TRI elements not possble for FLUID2_PRO\n");
    frint("GP_TRI",&(ele->e.f2pro->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of FLUID2_PRO element failed: integration\n");
    frchar("GP_ALT",buffer,&ierr);
/*
  integration for TRI-elements is distinguished into different cases. This is
  necessary to get the right integration parameters from FLUID_DATA.
  The flag for the integration case is saved in nGP[1]. For detailed informations
  see /fluid2/f2_intg.c.
*/
    switch(ele->e.f2pro->nGP[0])
    {
    case 1:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=0;
      else
        dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 3:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=1;
      else if (strncmp(buffer,"gaussrad",8)==0)
        ele->e.f2pro->nGP[1]=2;
      else
        dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 4:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=3;
      else
        dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 6:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=4;
      else if (strncmp(buffer,"gaussrad",8)==0)
        ele->e.f2pro->nGP[1]=5;
      else
        dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 7:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=6;
      else if (strncmp(buffer,"gaussrad",8)==0)
        ele->e.f2pro->nGP[1]=7;
      else
        dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 9:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=8;
      else
        dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 12:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=9;
      else
        dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 13:
      if (strncmp(buffer,"standard",8)==0)
        ele->e.f2pro->nGP[1]=10;
      else
        dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    default:
      dserror("Reading of FLUID2 element failed: integration points\n");
    } /* end switch(ele->e.f2pro->nGP[0]) */
    if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
  } /* endif (ele->numnp==3 || ele->numnp==6) */

/*------------------------------------------------------ read net algo */
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
        strncmp(buffer,"ALE",3)==0 ||
        strncmp(buffer,"Ale",3)==0 )
      dserror("ALE not possible for FLUID2_PRO!\n");
    /*ele->e.f2pro->is_ale=1;*/
    else if (strncmp(buffer,"euler",5)==0 ||
             strncmp(buffer,"EULER",5)==0 ||
             strncmp(buffer,"Euler",5)==0 )
      ele->e.f2pro->is_ale=0;
    else
      dserror("Reading of FLUID2_PRO element failed: Euler/Ale");
  }
  else
    dserror("Reading of FLUID2_PRO element failed: NA");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void f2pro_dis(ELEMENT* vele, ELEMENT* pele, INT numele, INT nodeshift)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("f2pro_dis");
#endif

  pele->e.f2pro = (FLUID2_PRO*)CCACALLOC(1,sizeof(FLUID2_PRO));
  pele->Id = vele->Id + numele;
  pele->mat = vele->mat;
  pele->e.f2pro->is_ale = vele->e.f2pro->is_ale;
  pele->e.f2pro->dm = vele->e.f2pro->dm;
  pele->distyp = vele->distyp;
  pele->numnp = vele->numnp;

  /* Switch to Taylor-Hood here. This way quite some nodes remain
   * unused. Another cheat... */
  switch (pele->numnp)
  {
  case 9:
    pele->numnp = 4;
    pele->distyp = quad4;
    break;
  case 4:
    /* do nothing. this is not inf-sup stable, but it might work... */
    break;
  default:
    dserror("unsupported number of nodes %d",pele->numnp);
  }

  pele->lm = (INT*)CCACALLOC(pele->numnp,sizeof(INT));

  for (i=0; i<pele->numnp; i++)
  {
    pele->lm[i] = vele->lm[i] + nodeshift;
  }

  pele->e.f2pro->nGP[0] = vele->e.f2pro->nGP[0];
  pele->e.f2pro->nGP[1] = vele->e.f2pro->nGP[1];

  pele->e.f2pro->other = vele;
  vele->e.f2pro->other = pele;

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
