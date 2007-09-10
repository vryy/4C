/*
<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
 */
#ifndef CCADISCRET
#ifdef D_FLUID3_IS

#include "fluid3_is.h"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


void f3is_inp(ELEMENT* ele)
{

  INT         i;
  INT         ierr = 0;
  INT         create_ale;
  char        buffer[50];
  static INT  cmat;

#ifdef DEBUG
  dstrc_enter("f3is_inp");
#endif

  /* allocate the element */
  ele->e.f3is = (FLUID3_IS*)CCACALLOC(1,sizeof(FLUID3_IS));

  /* read the element nodes */
  frchk("HEX8",&ierr);
  if (ierr==1)
  {
    ele->numnp=8;
    ele->distyp=hex8;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("HEX20",&ierr);
  if (ierr==1)
  {
    ele->numnp=20;
    ele->distyp=hex20;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("HEX27",&ierr);
  if (ierr==1)
  {
    ele->numnp=27;
    ele->distyp=hex27;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("HEX27",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  frchk("TET4",&ierr);
  if (ierr==1)
  {
    ele->numnp=4;
    ele->distyp=tet4;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
    frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");

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
    dserror("TET10 element not yet tested!!!\n");
    ele->numnp=10;
    ele->distyp=tet10;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
    frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  /* reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;

  /* read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1)
    dserror("Reading of FLUID3_IS element failed");
  if (ele->mat==0)
    dserror("No material defined for FLUID3_IS element");
  if (cmat==0)
    cmat=ele->mat;
  else if(ele->mat!=cmat)
    dserror("no different materials for fluid elements allowed!");

  /* read number of gaussian points */
  if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
  {
    frint_n("GP",&(ele->e.f3is->nGP[0]),3,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP\n");
  }

  /* read number of gaussian points for tetrahedral elements */
  if (ele->numnp==4 || ele->numnp==10)
  {
    frint("GP_TET",&(ele->e.f3is->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_TET\n");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_ALT\n");
    /*
     * integration for TET-elements is distinguished into different cases. This is
     * necessary to get the right integration parameters from FLUID_DATA.
     * The flag for the integration case is saved in nGP[1]. For detailed informations
     * see /fluid3/f3_intg.c.
     */
    switch(ele->e.f3is->nGP[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3is->nGP[1]=0;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3is->nGP[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ele->e.f3is->nGP[1]=2;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT\n");
        break;
      case 5:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3is->nGP[1]=3;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      default:
        dserror("Reading of FLUID3 element failed: integration points\n");
    }
  }

  /* read net algo */
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
        strncmp(buffer,"ALE",3)==0 ||
        strncmp(buffer,"Ale",3)==0 )
      ele->e.f3is->is_ale=1;
    else if (strncmp(buffer,"euler",5)==0 ||
        strncmp(buffer,"EULER",5)==0 ||
        strncmp(buffer,"Euler",5)==0 )
      ele->e.f3is->is_ale=0;
    else
      dserror("Reading of FLUID3 element failed: Euler/Ale\n");
  }
  else
    dserror("Reading of FLUID3 element failed: NA\n");

  /* read type of ale elements, created or read */
  if (!frreadyes("CA",&create_ale))
    create_ale = 0;

  if (create_ale == 1 && ele->e.f3is->is_ale == 1)
  {
    genprob.create_ale    = 1;
    ele->e.f3is->create_ale = 1;
  }
  else
  {
    /*genprob.create_ale    = 0;*/
    ele->e.f3is->create_ale = 0;
  }

  /* set initial value for free surface flag */
  ele->e.f3is->fs_on=0;

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif
#endif
