/*!----------------------------------------------------------------------
\file
\brief read fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID3 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3.h"
#include "fluid3_prototypes.h"
/*!---------------------------------------------------------------------
\brief read fluid3 element from input-file

<pre>                                                         genk 05/02
</pre>
\param  *ele	   ELEMENT	   (o)	   actual element
\return void                                                                       
\warning Node Numbers of TET4 are changed compered to the input
         file: node0=node1; node1=node0;
         This is necessary, since the GID-TET4 element is defined in
         a local left-system, which leads to a negative determinant 
         of the Jacobian matrix.
------------------------------------------------------------------------*/
void f3inp(ELEMENT *ele,INT counter)
{
INT  i;
INT  ierr=0;
INT  lmtmp;
char buffer[50];
static INT cmat;

#ifdef DEBUG 
dstrc_enter("f3inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.f3 = (FLUID3*)CCACALLOC(1,sizeof(FLUID3));
if (ele->e.f3==NULL) dserror("Allocation of element FLUID3 failed\n");
/*---------------------------------------------- read the element nodes */
frchk("HEX8",&ierr);
if (ierr==1)
{
   ele->numnp=8;
   ele->distyp=hex8;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("HEX20",&ierr);
if (ierr==1)
{
   dserror("HEX20 elements not implemented yet!!!\n"); 
   ele->numnp=20;
   ele->distyp=hex20;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("HEX27",&ierr);
if (ierr==1)
{
   dserror("HEX27 elements not implemented yet!!!\n");
   ele->numnp=27;
   ele->distyp=hex27;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
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
   /*-------------------------- rearrange element node numbers for tet4 */
   lmtmp=ele->lm[0];
   ele->lm[0]=ele->lm[1];
   ele->lm[1]=lmtmp;
}
frchk("TET10",&ierr); /* rerrangement??????? */
if (ierr==1)
{
   dserror("TET10 not implemented yet!!!\n");
   ele->numnp=10;
   ele->distyp=tet10;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
/*----------------------- reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of FLUID3 element failed\n");
if (ele->mat==0) dserror("No material defined for FLUID3 element\n");
if (counter==0) cmat=ele->mat;
else dsassert(ele->mat==cmat,"no different materials for fluid elements allowed!\n");
/*-------------------------------------------- read the gaussian points */
if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
{
   frint_n("GP",&(ele->e.f3->nGP[0]),3,&ierr);
   if (ierr!=1) dserror("Reading of FLUID3 element failed\n");
}   
/*----------------------- read gaussian points for tetrahedral elements */
if (ele->numnp==4 || ele->numnp==10)
{
/*   dserror("Tetrahedal Element not implemented yet!!!\n"); */
   frint("GP_TET",&(ele->e.f3->nGP[0]),&ierr);
   if (ierr!=1) dserror("Reading of FLUID3 element failed\n");
   frchar("GP_ALT",buffer,&ierr);
/*
integration for TET-elements is distinguished into different cases. This is
necessary to get the right integration parameters from FLUID_DATA. 
The flag for the integration case is saved in nGP[1]. For detailed informations
see /fluid3/f3_intg.c.
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
   if (ierr!=1) dserror("Reading of FLUID3 element failed: GP_ALT\n");
}  
/*------------------------------------------------------ read net algo */
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

/*----------------------------- set initial value for free surface flag */
ele->e.f3->fs_on=0;

/*-------------------------------------------------------- submesh data */
#ifdef FLUID3_ML
ele->e.f3->smisal = 0;   
#endif
/*--------------------------------------------------------------------- */
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3inp */

#endif
