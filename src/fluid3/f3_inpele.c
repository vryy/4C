/*-----------------------------------------------------------------------*/
/*!
\file
\brief read fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

 */
/*-----------------------------------------------------------------------*/


/*!
\addtogroup FLUID3
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3


#include "../headers/standardtypes.h"
#include "fluid3.h"
#include "fluid3_prototypes.h"
#include "../ale3/ale3.h"


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | general problem data                                               |
  | global variable GENPROB genprob is defined in global_control.c     |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*-----------------------------------------------------------------------*/
/*!
  \brief read fluid3 element from input-file

  \param  *ele      ELEMENT     (o)     actual element
  \param   counter  INT         (i)     how many elements have been read before

  \return void

  \author genk
  \date   05/02

 */
/*-----------------------------------------------------------------------*/
void f3inp(
    ELEMENT   *ele,
    INT        counter)
{

  INT         i;
  INT         ierr = 0;
#if 0
  INT         lmtmp;
#endif
  INT         create_ale;
  char        buffer[50];
  static INT  cmat;


#ifdef DEBUG
  dstrc_enter("f3inp");
#endif


  /* allocate the element */
  ele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));


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
  if (ierr!=1) dserror("Reading of FLUID3 element failed\n");
  if (ele->mat==0) dserror("No material defined for FLUID3 element\n");
  if (counter==0) cmat=ele->mat;
  else dsassert(ele->mat==cmat,"no different materials for fluid elements allowed!\n");


  /* read number of gaussian points */
  if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
  {
    frint_n("GP",&(ele->e.f3->nGP[0]),3,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP\n");
  }


  /* read number of gaussian points for tetrahedral elements */
  if (ele->numnp==4 || ele->numnp==10)
  {
    frint("GP_TET",&(ele->e.f3->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_TET\n");

    frchar("GP_ALT",buffer,&ierr);
    if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_ALT\n");
    /*
     * integration for TET-elements is distinguished into different cases. This is
     * necessary to get the right integration parameters from FLUID_DATA.
     * The flag for the integration case is saved in nGP[1]. For detailed informations
     * see /fluid3/f3_intg.c.
     */
    switch(ele->e.f3->nGP[0])
    {
      case 1:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3->nGP[1]=0;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      case 4:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3->nGP[1]=1;
        else if (strncmp(buffer,"gaussrad",8)==0)
          ele->e.f3->nGP[1]=2;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT\n");
        break;
      case 5:
        if (strncmp(buffer,"standard",8)==0)
          ele->e.f3->nGP[1]=3;
        else
          dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
        break;
      default:
        dserror("Reading of FLUID3 element failed: integration points\n");
    } /* end switch(ele->e.f3->nGP[0]) */
  }


  /* read net algo */
  frchar("NA",buffer,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"ale",3)==0 ||
        strncmp(buffer,"ALE",3)==0 ||
        strncmp(buffer,"Ale",3)==0 )
      ele->e.f3->is_ale=1;
    else if (strncmp(buffer,"euler",5)==0 ||
        strncmp(buffer,"EULER",5)==0 ||
        strncmp(buffer,"Euler",5)==0 )
      ele->e.f3->is_ale=0;
    else
      dserror("Reading of FLUID3 element failed: Euler/Ale\n");
  }
  else
    dserror("Reading of FLUID3 element failed: NA\n");


  /* read type of ale elements, created or read */
  if (!frreadyes("CA",&create_ale))
    create_ale = 0;

  if (create_ale == 1 && ele->e.f3->is_ale == 1)
  {
    genprob.create_ale    = 1;
    ele->e.f3->create_ale = 1;
  }
  else
  {
    /*genprob.create_ale    = 0;*/
    ele->e.f3->create_ale = 0;
  }



  /* set initial value for free surface flag */
  ele->e.f3->fs_on=0;


  /* submesh data */
#ifdef FLUID3_ML
  ele->e.f3->smisal = 0;
#endif


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
} /* end of f3inp */






/*-----------------------------------------------------------------------*/
/*!
  \brief create ale elements

  \param  *ele0       ELEMENT     (i) the source (fluid) element
  \parem  *ele1       ELEMENT     (o) the created ale element
  \param   numele     INT         (i) total number of elements in the problem
  \param   nodeshift  INT         (i) total number of nodes in the problem

  \return void

  \author mn
  \date   08/05

 */
/*-----------------------------------------------------------------------*/
void f3_createale(
    ELEMENT *ele0,
    ELEMENT *ele1,
    INT      numele,
    INT      nodeshift)
{

#ifdef D_ALE


  INT       j;             /* a counter */


#ifdef DEBUG
  dstrc_enter("f3_createale");
#endif

  /* allocate the element */
  ele1->e.ale3 = (ALE3*)CCACALLOC(1,sizeof(ALE3));
  if (ele1->e.ale3==NULL) dserror("Allocation of element ALE failed");


  /* set some general data */
  ele1->Id             = ele0->Id + numele;
  ele1->eltyp          = el_ale3;

  /* set data depending on dis type */
  switch (ele0->distyp)
  {
    case hex8:
      ele1->distyp  = hex8;
      ele1->numnp   = 8;
      ele1->lm      = (INT*)CCACALLOC(ele1->numnp,sizeof(INT));
      if (ele1->lm==NULL) dserror("Allocation of lm in ELEMENT failed");

      for (j=0;j<ele1->numnp;j++)
      {
        ele1->lm[j] = ele0->lm[j] + nodeshift;
      }

      ele1->mat            = 1;
      ele1->e.ale3->nGP[0] = 2;
      ele1->e.ale3->nGP[1] = 2;
      ele1->e.ale3->nGP[2] = 2;
      ele1->e.ale3->jacobi = 0;
      ele1->e.ale3->hex20_red = 0;
      break;

    case tet4:
      ele1->distyp  = tet4;
      ele1->numnp   = 4;
      ele1->lm      = (INT*)CCACALLOC(ele1->numnp,sizeof(INT));
      if (ele1->lm==NULL) dserror("Allocation of lm in ELEMENT failed");

      for (j=0;j<ele1->numnp;j++)
      {
        ele1->lm[j] = ele0->lm[j] + nodeshift;
      }

      ele1->mat            = 1;
      ele1->e.ale3->nGP[0] = 1;
      ele1->e.ale3->nGP[1] = 1;
      ele1->e.ale3->nGP[2] = 1;
      ele1->e.ale3->jacobi = 0;
      break;

    case hex20:
      ele1->distyp  = hex20;
      ele1->numnp   = 20;
      ele1->lm      = (INT*)CCACALLOC(ele1->numnp,sizeof(INT));
      if (ele1->lm==NULL) dserror("Allocation of lm in ELEMENT failed");

      for (j=0;j<ele1->numnp;j++)
      {
        ele1->lm[j] = ele0->lm[j] + nodeshift;
      }

      ele1->mat            = 1;
      ele1->e.ale3->nGP[0] = 3;
      ele1->e.ale3->nGP[1] = 3;
      ele1->e.ale3->nGP[2] = 3;
      ele1->e.ale3->jacobi = 0;
      ele1->e.ale3->hex20_red = 1;
      break;

    default:
      dserror("Creation of ale elements only for hex8 and hex20 elements!!\n");
  }  /* switch (ele0->distyp) */


  /* make ale element known to the fluid and vice versa */
  ele1->e.ale3->fluid_ele = ele0;
  ele0->e.f3->ale_ele     = ele1;


#ifdef DEBUG
  dstrc_exit();
#endif


#endif  /* ifdef D_ALE */

  return;
} /* end of f3_createale */



#endif



