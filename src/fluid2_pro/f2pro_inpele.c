/*!----------------------------------------------------------------------
\file
\brief read fluid2_pro element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
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

/*------------------------------------------------ allocate the element */
ele->e.f2pro = (FLUID2_PRO*)CCACALLOC(1,sizeof(FLUID2_PRO));
if (ele->e.f2pro==NULL) dserror("Allocation of element FLUID2_PRO failed\n");
/*---------------------------------------------- read the element nodes */
frchk("QUAD4",&ierr);
if (ierr==1)
{
   dserror("quad4 for FLUID2_PRO not implemented yet!\n");
/*   ele->numnp=4;
   ele->distyp=quad4;
   ele->e.f2->ntyp=1;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n"); */
}
frchk("QUAD9",&ierr);
if (ierr==1)
{
   ele->numnp=9;
   ele->distyp=quad9;
   ele->e.f2pro->ntyp=1;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
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
       dserror("Reading of FLUID2_PRO element failed: Euler/Ale\n");
}
else
   dserror("Reading of FLUID2_PRO element failed: NA\n");

frchar("DISMODE",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"Q2Q1",4)==0)
       ele->e.f2pro->dm=dm_q2q1;
   else
       dserror("Reading of FLUID2_PRO element failed: DISMODE\n");
}
else
   dserror("Reading of FLUID2_PRO element failed: DISMODE\n");
/*----------------------------- her now read all the other stuff needed */


/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2pro_inp */

/*!---------------------------------------------------------------------
\brief create elements of other discretisations

<pre>                                                         genk 08/02
</pre>
\param  *ele0    ELEMENT	   (i)
\parem  *ele1    ELEMENT           (o)
\return void

------------------------------------------------------------------------*/
void f2pro_dis(
    ELEMENT *ele0,
    ELEMENT *ele1,
    INT      numele,
    INT      numnode)
{
INT       j;             /* simply a counter                            */

#ifdef DEBUG
dstrc_enter("f2pro_dis");
#endif

/*------------------------------------------------ allocate the element */
ele1->e.f2pro = (FLUID2_PRO*)CCACALLOC(1,sizeof(FLUID2_PRO));
if (ele1->e.f2pro==NULL) dserror("Allocation of element FLUID2_PRO failed\n");

/*----------------------------------------------- set some general data */
ele1->Id=ele0->Id + numele;
ele1->mat=ele0->mat;
ele1->e.f2pro->is_ale=ele0->e.f2pro->is_ale;
ele1->e.f2pro->dm=ele0->e.f2pro->dm;

/*--------------------------------------- set data depending on dismode */
switch (ele1->e.f2pro->dm)
{
case dm_q2q1:
   ele1->numnp=4;
   ele1->distyp=quad4;
   ele1->e.f2pro->ntyp=1;
   ele1->lm = (INT*)CCACALLOC(ele1->numnp,sizeof(INT));
   if (ele1->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   for (j=0;j<ele1->numnp;j++)
   {
      ele1->lm[j] = ele0->lm[j] + numnode;
   } /* end of loop over all nodes */
   ele1->e.f2pro->nGP[0]=2;
   ele1->e.f2pro->nGP[1]=2;
break;
default:
   dserror("DISMODE unknown!\n");
} /* end switch(ele1->e.f2pro->dm) */

/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2pro_dis */

#endif
/*! @} (documentation module close)*/
