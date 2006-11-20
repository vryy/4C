/*!----------------------------------------------------------------------
\file
\brief input of fluid2_is element

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifdef D_FLUID2_IS

#include "fluid2_is.h"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


void f2is_inp(ELEMENT* ele)
{
  INT        i;             /* simply a counter                           */
  INT        ierr=0;        /* error flag                                 */
  INT        create_ale;
  char      buffer[50];
  static INT cmat;

#ifdef DEBUG
  dstrc_enter("f2is_inp");
#endif

  ele->e.f2is = (FLUID2_IS*)CCACALLOC(1,sizeof(FLUID2_IS));

  /* read the element nodes */
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    ele->numnp=4;
    ele->distyp=quad4;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }
  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    ele->numnp=9;
    ele->distyp=quad9;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }
  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    ele->numnp=8;
    ele->distyp=quad8;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    ele->numnp=3;
    ele->distyp=tri3;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }
  frchk("TRI6",&ierr);
  if (ierr==1)
  {
    ele->numnp=6;
    ele->distyp=tri6;
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
  }

  /* reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;

  /* read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1)
    dserror("Reading of FLUID2_IS element failed");
  if (ele->mat==0)
    dserror("No material defined for FLUID2_IS element");
  if (cmat==0)
    cmat=ele->mat;
  else if(ele->mat!=cmat)
    dserror("no different materials for fluid elements allowed!");

  /*-------------------------------------------- read the gaussian points */
  if (ele->numnp==4 || ele->numnp==8 || ele->numnp==9)
  {
    frint_n("GP",&(ele->e.f2is->nGP[0]),2,&ierr);
    if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
  }

  /*-------------------------- read gaussian points for triangle elements */
  if (ele->numnp==3 || ele->numnp==6)
  {
    frint("GP_TRI",&(ele->e.f2is->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
    frchar("GP_ALT",buffer,&ierr);
/*
  integration for TRI-elements is distinguished into different cases. This is
  necessary to get the right integration parameters from FLUID_DATA.
  The flag for the integration case is saved in nGP[1]. For detailed informations
  see /fluid2/f2_intg.c.
*/
    switch(ele->e.f2is->nGP[0])
    {
    case 1:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=0;
      else
	dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 3:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=1;
      else if (strncmp(buffer,"gaussrad",8)==0)
	ele->e.f2is->nGP[1]=2;
      else
	dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 4:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=3;
      else
	dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 6:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=4;
      else if (strncmp(buffer,"gaussrad",8)==0)
	ele->e.f2is->nGP[1]=5;
      else
	dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 7:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=6;
      else if (strncmp(buffer,"gaussrad",8)==0)
	ele->e.f2is->nGP[1]=7;
      else
	dserror("Reading of FLUID2 element failed: GP_ALT\n");
      break;
    case 9:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=8;
      else
	dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 12:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=9;
      else
	dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    case 13:
      if (strncmp(buffer,"standard",8)==0)
	ele->e.f2is->nGP[1]=10;
      else
	dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
    default:
      dserror("Reading of FLUID2 element failed: integration points\n");
    } /* end switch(ele->e.f2is->nGP[0]) */
    if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
  } /* endif (ele->numnp==3 || ele->numnp==6) */

/*------------------------------------------------ set default velues: */
  ele->e.f2is->is_ale=0;

/*------------------------------------------------------ read net algo */
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
	strncmp(buffer,"ALE",3)==0 ||
	strncmp(buffer,"Ale",3)==0 )
      ele->e.f2is->is_ale=1;
    else if (strncmp(buffer,"euler",5)==0 ||
	     strncmp(buffer,"EULER",5)==0 ||
	     strncmp(buffer,"Euler",5)==0 )
      ele->e.f2is->is_ale=0;
    else
      dserror("Reading of FLUID2 element failed: Euler/Ale\n");
  }
  else
    dserror("Reading of FLUID2 element failed: NA\n");

  /* read type of ale elements, created or read */
  frint("CA",&(create_ale),&ierr);

  /*if (ierr!=1) dserror("Reading of FLUID3 element failed: flag CA not available\n");*/
  if (ierr!=1)
    create_ale = 0;

  if (create_ale == 1 && ele->e.f2is->is_ale == 1)
  {
    genprob.create_ale    = 1;
    ele->e.f2is->create_ale = 1;
  }
  else
  {
    genprob.create_ale    = 0;
    ele->e.f2is->create_ale = 0;
  }

/*----------------------------- set initial value for free surface flag */
  ele->e.f2is->fs_on=0;

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}

#endif
