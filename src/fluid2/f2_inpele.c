/*!----------------------------------------------------------------------
\file
\brief read fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h" 
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!--------------------------------------------------------------------- 
\brief read fluid2 element from input-file

<pre>                                                         genk 03/02
</pre>
\param  *ele	   ELEMENT	   (o)	   actual element
\param   counter   INT             (i)     counter for first element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_inp(ELEMENT *ele, INT counter)
{
INT        i;             /* simply a counter                           */
INT        ndum;          /* dummy value                                */
INT        ierr=0;        /* error flag                                 */
INT        itaumu;        /*                                            */ 
INT        itaump;	  /*						*/
INT        itauc;	  /* element flags                              */
char      buffer[50];
static INT cmat;

#ifdef DEBUG 
dstrc_enter("f2inp");
#endif

/*------------------------------------------------ allocate the element */      
ele->e.f2 = (FLUID2*)CCACALLOC(1,sizeof(FLUID2));
if (ele->e.f2==NULL) dserror("Allocation of element FLUID2 failed\n");
/*---------------------------------------------- read the element nodes */
frchk("QUAD4",&ierr);
if (ierr==1)
{
   ele->numnp=4;
   ele->distyp=quad4;
   ele->e.f2->ntyp=1;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("QUAD9",&ierr);
if (ierr==1)
{
   ele->numnp=9;
   ele->distyp=quad9;
   ele->e.f2->ntyp=1;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("QUAD8",&ierr);
if (ierr==1)
{
   ele->numnp=8;
   ele->distyp=quad8;
   ele->e.f2->ntyp=1;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("QUAD8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("TRI3",&ierr);
if (ierr==1)
{
   ele->numnp=3;
   ele->distyp=tri3;
   ele->e.f2->ntyp=2;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
#ifdef D_LS
   ele->e.f2->is_ls = 1;
#endif
}
frchk("TRI6",&ierr);
if (ierr==1)
{
   ele->numnp=6;
   ele->distyp=tri6;
   ele->e.f2->ntyp=2;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
} 
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of FLUID2 element failed\n");
if (ele->mat==0) dserror("No material defined for FLUID2 element\n");
if (counter==0) cmat=ele->mat;
else dsassert(ele->mat==cmat,"no different materials for fluid elements allowed!\n");
/*-------------------------------------------- read the gaussian points */
if (ele->numnp==4 || ele->numnp==8 || ele->numnp==9)
{
   frint_n("GP",&(ele->e.f2->nGP[0]),2,&ierr);
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");  
}   
/*-------------------------- read gaussian points for triangle elements */
if (ele->numnp==3 || ele->numnp==6)
{
   frint("GP_TRI",&(ele->e.f2->nGP[0]),&ierr);   
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
   frchar("GP_ALT",buffer,&ierr);
/*
integration for TRI-elements is distinguished into different cases. This is
necessary to get the right integration parameters from FLUID_DATA. 
The flag for the integration case is saved in nGP[1]. For detailed informations
see /fluid2/f2_intg.c.
*/   
   switch(ele->e.f2->nGP[0])
   {
   case 1:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=0;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");
      break;
   case 3:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=1;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=2;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT\n");	 
      break;
   case 4:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=3;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");  
      break;
   case 6:
       if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=4;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=5;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT\n");  
      break;
   case 7:
       if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=6;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=7;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT\n");  
      break;
   case 9:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=8;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");    
      break;
   case 12:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=9;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");   
      break;
   case 13:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=10;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!\n");   
      break;
   default:
      dserror("Reading of FLUID2 element failed: integration points\n");
   } /* end switch(ele->e.f2->nGP[0]) */
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration\n");
} /* endif (ele->numnp==3 || ele->numnp==6) */

/*------------------------------------------------------ read net algo */
frchar("NA",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"ale",3)==0 ||
       strncmp(buffer,"ALE",3)==0 ||
       strncmp(buffer,"Ale",3)==0 )
       ele->e.f2->is_ale=1;
   else if (strncmp(buffer,"euler",5)==0 ||
            strncmp(buffer,"EULER",5)==0 ||
            strncmp(buffer,"Euler",5)==0 )
       ele->e.f2->is_ale=0;
   else 
       dserror("Reading of FLUID2 element failed: Euler/Ale\n");
}
else
   dserror("Reading of FLUID2 element failed: NA\n");

/*-------------------------------------------- read flag for turbulence */
frchar("TURBULENCE",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"No",2)==0)
      ele->e.f2->turbu=0;
   else if (strncmp(buffer,"algebraic",9)==0)
      ele->e.f2->turbu=1;
   else if (strncmp(buffer,"kappa-eps",9)==0)
      ele->e.f2->turbu=2;
   else if (strncmp(buffer,"kappa-omega",11)==0)
      ele->e.f2->turbu=3;
   else
      dserror("Reading of FLUID2 element failed: TURBULENCE\n"); 
}
else 
  dserror("Reading of FLUID2 element failed: TURBULENCE\n");

/*----------------------------- set initial value for free surface flag */
ele->e.f2->fs_on=0;
   
/*-------------------------------------------------------- submesh data */
ele->e.f2->smisal = 0;   
/*--------------------------------------------------------------------- */
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2inp */


#endif
/*! @} (documentation module close)*/
