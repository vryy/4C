/*!----------------------------------------------------------------------
\file
\brief read fluid3 element

------------------------------------------------------------------------*/
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
void f3inp(ELEMENT *ele,int counter)
{
int  i;
int  ierr=0;
int  quad;
int  lmtmp;
long int  topology[100];
char *colpointer;
char buffer[50];
int  ndum;       /* dummy value */
int  ihelem;
int  itaumu;
int  itaump;
int  itauc;
static int cmat;

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
   ele->e.f3->ntyp=1;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
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
   ele->e.f3->ntyp=1;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
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
   ele->e.f3->ntyp=1;
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed\n");
}
frchk("TET4",&ierr);
if (ierr==1)
{
   ele->numnp=4;
   ele->distyp=tet4;
   ele->e.f3->ntyp=2;   
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
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
   ele->e.f3->ntyp=2; 
   ele->lm = (int*)CCACALLOC(ele->numnp,sizeof(int));
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

/*-------------------------------------------------- read stabilisation */
frchar("ISTABI",buffer,&ierr);
if (ierr==1)
{
  if (strncmp(buffer,"yes",3)==0)
     ele->e.f3->istabi=1;
  else if (strncmp(buffer,"no",2)==0)
     ele->e.f3->istabi=0;
  else
    dserror("Reading of FLUID3 element failed: ISTABI\n"); 
}    
else
  dserror("Reading of FLUID3 element failed: ISTABI\n"); 
  
/*----------------------------------------------- read stab flag iadvec */
frchar("IADVEC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f3->iadvec=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f3->iadvec=0;  
   else 
       dserror("Reading of FLUID3 element failed: IADVEC\n");                
}
else
   dserror("Reading of FLUID3 element failed: IADVEC\n");
/*----------------------------------------------- read stab flag ipres */
frchar("IPRES",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f3->ipres=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f3->ipres=0;
   else 
       dserror("Reading of FLUID3 element failed: IPRES\n");      
}
else
   dserror("Reading of FLUID3 element failed: IPRES\n");
/*----------------------------------------------- read stab flag ivisc */
frchar("IVISC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"GLS-",4)==0)
       ele->e.f3->ivisc=1;
   else if (strncmp(buffer,"GLS+",4)==0)
       ele->e.f3->ivisc=2;
   else if (strncmp(buffer,"no",2)==0) 
       ele->e.f3->ivisc=0;
   else 
       dserror("Reading of FLUID3 element failed: IVISC\n");       
}
else
   dserror("Reading of FLUID3 element failed: IVISC\n");
/*----------------------------------------------- read stab flag icont */
frchar("ICONT",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f3->icont=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f3->icont=0;
   else 
       dserror("Reading of FLUID3 element failed: ICONT\n");       
}
else
   dserror("Reading of FLUID3 element failed: ICONT\n");  
/*---------------------------------------------------- read stab norm */
frchar("NORM_P",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"L_2",3)==0)
       ele->e.f3->norm_p=2;
   else if (strncmp(buffer,"L_1",3)==0)
       ele->e.f3->norm_p=1;
   else 
       dserror("Reading of FLUID3 element failed: NORM_P\n");       
}
else
   dserror("Reading of FLUID3 element failed: NORM_P\n"); 
/*----------------------------------------------- read stab flag ninths */
frchar("NINTHS",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"at_center",9)==0)
      ele->e.f3->ninths=1;
   else if (strncmp(buffer,"every_intpt",11)==0)
      ele->e.f3->ninths=2;
   else
      dserror("Reading of FLUID3 element failed: NINTHS\n"); 
}
else 
  dserror("Reading of FLUID3 element failed: NINTHS\n");
/*----------------------------------------------- read stab flag istapc */
frchar("ISTAPC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"at_center",9)==0)
      ele->e.f3->istapc=1;
   else if (strncmp(buffer,"every_intpt",11)==0)
      ele->e.f3->istapc=2;
   else
      dserror("Reading of FLUID3 element failed: ISTAPC\n"); 
}
else 
  dserror("Reading of FLUID3 element failed: ISTAPC\n");      
/*----------------------------------------------- read stab flag istapa */
frint("ISTAPA",&(ele->e.f3->istapa),&ierr); 
if (ierr!=1) dserror("Reading of FLUID3 element failed: ISTAPA\n");
/*----------------------------------------------- read stab flag istapc */
frint("ISTAPC",&(ele->e.f3->istapc),&ierr); 
if (ierr!=1) dserror("Reading of FLUID3 element failed: ISTAPC\n");
/*--------------------------------------------------- read stab flag mk */
frint("MK",&(ele->e.f3->mk),&ierr); 
if (ierr!=1) dserror("Reading of FLUID3 element failed: MK\n");
/*----------------------------------------------- read stab flag ihelem */
frint("IHELEM",&ihelem,&ierr); 
if (ierr!=1) dserror("Reading of FLUID3 element failed: IHELEM\n");
/*---------------------------------------------- read stab const c_lamb */
frdouble("C_LAMB",&(ele->e.f3->clamb),&ierr); 
if (ierr!=1) dserror("Reading of FLUID3 element failed: C_LAMB\n");

/*-------------------------------------------- set some additional data */
/* initialisation */
ele->e.f3->istrle  = 0;
ele->e.f3->ivol   = 0;
ele->e.f3->iduring = 0;
ele->e.f3->idiaxy  = 0;
itaumu  = 0;
itaump  = 0;
itauc	= 0;

math_intextract(ihelem,&ndum,
               &(ele->e.f3->ihele[0]),&(ele->e.f3->ihele[1]),&(ele->e.f3->ihele[2]));

for(i=0;i<3;i++)
{
   if (ele->e.f3->ihele[i]==5)
   {
      ele->e.f3->istrle  = 1;
      if (ele->e.f3->iadvec!=0 && ele->e.f3->ninths==1)
         itaumu = -1;
      if (ele->e.f3->iadvec!=0 && ele->e.f3->ninths!=1)
         itaumu = 1;
      if (ele->e.f3->ipres!=0 && ele->e.f3->ninths==1)
         itaump = -1;
      if (ele->e.f3->ipres!=0 && ele->e.f3->ninths!=1)
         itaump = 1;
      if (ele->e.f3->icont!=0 && ele->e.f3->ninths==1)
         itauc = -1;
      if (ele->e.f3->icont!=0 && ele->e.f3->ninths!=1)
         itauc = 1;     
   }
   else if (ele->e.f3->ihele[i]!=0)
   {
      ele->e.f3->ivol = 1;
      if (ele->e.f3->iadvec!=0 && ele->e.f3->istapc==1)
         itaumu = -1;
      if (ele->e.f3->iadvec!=0 && ele->e.f3->istapc!=1)
         itaumu = 1;
      if (ele->e.f3->ipres!=0 && ele->e.f3->istapc==1)
         itaump = -1;
      if (ele->e.f3->ipres!=0 && ele->e.f3->istapc!=1)
         itaump = 1;
      if (ele->e.f3->icont!=0 && ele->e.f3->istapc==1)
         itauc = -1;
      if (ele->e.f3->icont!=0 && ele->e.f3->istapc!=1)
         itauc = 1;
	 
      if (ele->e.f3->ihele[i]==4)
         ele->e.f3->idiaxy = 1;         
   }
}

if (ele->e.f3->istrle==1 && ele->e.f3->ninths!=1)
   ele->e.f3->iduring = 1;
if (ele->e.f3->ivol==1 && ele->e.f3->istapc!=1)
   ele->e.f3->iduring = 1;

/*------------------------------------------- store data at the element */
ele->e.f3->itau[0] = itaumu;
ele->e.f3->itau[1] = itaump;
ele->e.f3->itau[2] = itauc;

/*----------------------------- set initial value for free surface flag */
ele->e.f3->fs_on=0;

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3inp */

#endif
