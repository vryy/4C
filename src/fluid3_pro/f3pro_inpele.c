/*!----------------------------------------------------------------------
\file
\brief read fluid3_pro element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID3_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3_PRO
#include "../headers/standardtypes.h"
#include "fluid3pro_prototypes.h"
#include "fluid3pro.h"
/*!---------------------------------------------------------------------
\brief read fluid3 element from input-file

<pre>                                                         genk 03/02
</pre>
\param  *ele	   ELEMENT	   (o)	   actual element
\return void

------------------------------------------------------------------------*/
void f3pro_inp(ELEMENT *ele)
{
  INT        i;             /* simply a counter                           */
  INT        ierr=0;        /* error flag                                 */
  char       buffer[50];

#ifdef DEBUG
  dstrc_enter("f3pro_inp");
#endif

/*------------------------------------------------ allocate the element */
  ele->e.f3pro = (FLUID3_PRO*)CCACALLOC(1,sizeof(FLUID3_PRO));

  /* read the element nodes */
  frchk("HEX8",&ierr);
  if (ierr==1)
  {
    ele->numnp=8;
    ele->distyp=hex8;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");

    ele->e.f3pro->dm=dm_q1p0;

    /*
     * Lets not allocate too small pieces of memory. Two pointers
     * share one array peacefully. :) */
    ele->e.f3pro->press = (DOUBLE*)CCACALLOC(3, sizeof(DOUBLE));
    ele->e.f3pro->pressm= &(ele->e.f3pro->press[1]);
    ele->e.f3pro->phi   = &(ele->e.f3pro->press[2]);
    ele->e.f3pro->dof   = (INT*)   CCACALLOC(2, sizeof(INT));
    ele->e.f3pro->ldof  = &(ele->e.f3pro->dof[1]);
  }


  frchk("HEX20",&ierr);
  if (ierr==1)
  {
    dserror("currently not supported");
    ele->numnp=20;
    ele->distyp=hex20;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
  }

  frchk("HEX27",&ierr);
  if (ierr==1)
  {
    ele->numnp=27;
    ele->distyp=hex27;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("HEX27",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");

    ele->e.f3pro->dm=dm_q2pm1;

    /*
     * Lets not allocate too small pieces of memory. Two pointers
     * share one array peacefully. :) */
    ele->e.f3pro->press = (DOUBLE*)CCACALLOC(12, sizeof(DOUBLE));
    ele->e.f3pro->pressm= &(ele->e.f3pro->press[4]);
    ele->e.f3pro->phi   = &(ele->e.f3pro->press[8]);
    ele->e.f3pro->dof   = (INT*)   CCACALLOC(8, sizeof(INT));
    ele->e.f3pro->ldof  = &(ele->e.f3pro->dof[4]);
  }

  frchk("TET4",&ierr);
  if (ierr==1)
  {
    dserror("currently not supported");
    ele->numnp=4;
    ele->distyp=tet4;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");

#if 0
    /* rearrange element node numbers for tet4 */
    lmtmp=ele->lm[0];
    ele->lm[0]=ele->lm[1];
    ele->lm[1]=lmtmp;
#endif
  }


  frchk("TET10",&ierr); /* rerrangement??????? */
  if (ierr==1)
  {
    dserror("currently not supported");
    ele->numnp=10;
    ele->distyp=tet10;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
    frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
  }



  /* reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;


  /* read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1) dserror("Reading of FLUID3 element failed");
  if (ele->mat==0) dserror("No material defined for FLUID3 element");
#if 0
  if (counter==0) cmat=ele->mat;
  else dsassert(ele->mat==cmat,"no different materials for fluid elements allowed!");
#endif

  /* read number of gaussian points */
  if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
  {
    frint_n("GP",&(ele->e.f3pro->nGP[0]),3,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP");
  }


  /* read number of gaussian points for tetrahedral elements */
  if (ele->numnp==4 || ele->numnp==10)
  {
    frint("GP_TET",&(ele->e.f3pro->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_TET");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_ALT");
    /*
     * integration for TET-elements is distinguished into different cases. This is
     * necessary to get the right integration parameters from FLUID_DATA.
     * The flag for the integration case is saved in nGP[1]. For detailed informations
     * see /fluid3/f3_intg.c.
     */
    switch(ele->e.f3pro->nGP[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3pro->nGP[1]=0;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3pro->nGP[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ele->e.f3pro->nGP[1]=2;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT");
        break;
      case 5:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3pro->nGP[1]=3;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!");
        break;
      default:
        dserror("Reading of FLUID3 element failed: integration points");
    } /* end switch(ele->e.f3->nGP[0]) */
  }


  /* read net algo */
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
        strncmp(buffer,"ALE",3)==0 ||
        strncmp(buffer,"Ale",3)==0 )
      ele->e.f3pro->is_ale=1;
    else if (strncmp(buffer,"euler",5)==0 ||
        strncmp(buffer,"EULER",5)==0 ||
        strncmp(buffer,"Euler",5)==0 )
      ele->e.f3pro->is_ale=0;
    else
      dserror("Reading of FLUID3 element failed: Euler/Ale");
  }
  else
    dserror("Reading of FLUID3 element failed: NA");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


void f3pro_dis(ELEMENT* vele, ELEMENT* pele, INT numele, INT nodeshift)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("f3pro_dis");
#endif

  pele->e.f3pro = (FLUID3_PRO*)CCACALLOC(1,sizeof(FLUID3_PRO));
  pele->Id = vele->Id + numele;
  pele->mat = vele->mat;
  pele->e.f3pro->is_ale = vele->e.f3pro->is_ale;
  pele->e.f3pro->dm = vele->e.f3pro->dm;
  pele->numnp = vele->numnp;
  pele->distyp = vele->distyp;

  /* Switch to "Taylor-Hood" here. This way quite some nodes remain
   * unused. Another cheat... */
  switch (pele->numnp)
  {
  case 27:
    pele->numnp = 8;
    pele->distyp = hex8;
    break;
  case 8:
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

  pele->e.f3pro->nGP[0] = vele->e.f3pro->nGP[0];
  pele->e.f3pro->nGP[1] = vele->e.f3pro->nGP[1];

  pele->e.f3pro->other = vele;
  vele->e.f3pro->other = pele;

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
#endif
