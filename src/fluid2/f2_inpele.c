#include "../headers/standardtypes.h" 
#include "fluid2_prototypes.h"
/*----------------------------------------------------------------------*
 | read fluid2 element                                    genk 3/02     |
 *----------------------------------------------------------------------*/
void f2_inp(ELEMENT *ele)
{
int  i;
int  ndum;       /* dummy value */
int  ihelem;
int  test;
int  ierr=0;
int  quad;
int  counter;
int  itaumu;
int  itaump;
int  itauc;
long int  topology[100];
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("f2inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.f2 = (FLUID2*)CALLOC(1,sizeof(FLUID2));
if (ele->e.f2==NULL) dserror("Allocation of element FLUID2 failed");
/*---------------------------------------------- read the element nodes */
frchk("QUAD4",&ierr);
if (ierr==1)
{
   ele->numnp=4;
   ele->distyp=quad4;
   ele->e.f2->ntyp=1;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/* ############### NOT CHECKED!!! ############################ */
frchk("QUAD9",&ierr);
if (ierr==1)
{
   ele->numnp=9;
   ele->distyp=quad9;
   ele->e.f2->ntyp=1;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("QUAD9",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI3",&ierr);
if (ierr==1)
{
   ele->numnp=3;
   ele->distyp=tri3;
   ele->e.f2->ntyp=2;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TRI6",&ierr);
if (ierr==1)
{
   ele->numnp=6;
   ele->distyp=tri6;
   ele->e.f2->ntyp=2;
   ele->lm = (int*)CALLOC(ele->numnp,sizeof(int));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TRI6",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
} 
/* ^^^^^^^^^^^^^^^                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
/* ############### NOT CHECKED!!! ############################ */
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of FLUID2 element failed");
if (ele->mat==0) dserror("No material defined for FLUID2 element");
/*-------------------------------------------- read the gaussian points */
if (ele->numnp==4 || ele->numnp==9)
{
   frint_n("GP",&(ele->e.f2->nGP[0]),2,&ierr);
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration");  
}   
/*-------------------------- read gaussian points for triangle elements */
if (ele->numnp==3 || ele->numnp==6)
{
   frint("GP_TRI",&(ele->e.f2->nGP[0]),&ierr);   
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration");
   frchar("GP_ALT",buffer,&ierr);
/*
integration for TRI-elements is distinguished into different cases. This is
necessary to get the right integration parameters from F2_DATA. 
The flag for the integration case is saved in nGP[1]. For detailed informations
see /fluid2/f2_intg.c.
*/   
   switch(ele->e.f2->nGP[0])
   {
   case 1:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=0;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!");
   case 3:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=1;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=2;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT");	 
   case 4:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=3;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!");  
   case 6:
       if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=4;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=5;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT");  
   case 7:
       if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=6;
      else if (strncmp(buffer,"gaussrad",8)==0)
         ele->e.f2->nGP[1]=7;
      else
         dserror("Reading of FLUID2 element failed: GP_ALT");  
   case 9:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=8;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!");    
   case 12:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=9;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!");   
   case 13:
      if (strncmp(buffer,"standard",8)==0)
         ele->e.f2->nGP[1]=10;
      else
         dserror("Reading of FLUID2 element failed: gauss-radau not possible!");   
   default:
      dserror("Reading of FLUID2 element failed: integration");
   }
   if (ierr!=1) dserror("Reading of FLUID2 element failed: integration");
}
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
       dserror("Reading of FLUID2 element failed: Euler/Ale");
}
else
   dserror("Reading of FLUID2 element failed: NA");
/*-------------------------------------------------- read stabilisation */
frchar("ISTABI",buffer,&ierr);
if (ierr==1)
{
  if (strncmp(buffer,"yes",3)==0)
     ele->e.f2->istabi=1;
  else if (strncmp(buffer,"no",2)==0)
     ele->e.f2->istabi=0;
  else
    dserror("Reading of FLUID2 element failed: ISTABI"); 
}

/*----------------------------------------------- read stab flag iadvec */
frchar("IADVEC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f2->iadvec=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f2->iadvec=0;  
   else 
       dserror("Reading of FLUID2 element failed: IADVEC");                
}
else
   dserror("Reading of FLUID2 element failed: IADVEC");
/*----------------------------------------------- read stab flag ipres */
frchar("IPRES",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f2->ipres=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f2->ipres=0;
   else 
       dserror("Reading of FLUID2 element failed: IPRES");      
}
else
   dserror("Reading of FLUID2 element failed: IPRES");
/*----------------------------------------------- read stab flag ivisc */
frchar("IVISC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"GLS-",4)==0)
       ele->e.f2->ivisc=1;
   else if (strncmp(buffer,"GLS+",4)==0)
       ele->e.f2->ivisc=2;
   else if (strncmp(buffer,"no",2)==0) 
       ele->e.f2->ivisc=0;
   else 
       dserror("Reading of FLUID2 element failed: IVISC");       
}
else
   dserror("Reading of FLUID2 element failed: IVISC");
/*----------------------------------------------- read stab flag icont */
frchar("ICONT",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"yes",3)==0)
       ele->e.f2->icont=1;
   else if (strncmp(buffer,"no",2)==0)
       ele->e.f2->icont=0;
   else 
       dserror("Reading of FLUID2 element failed: ICONT");       
}
else
   dserror("Reading of FLUID2 element failed: ICONT");  
/*---------------------------------------------------- read stab norm */
frchar("NORM_P",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"L_2",3)==0)
       ele->e.f2->norm_p=2;
   else if (strncmp(buffer,"L_1",3)==0)
       ele->e.f2->norm_p=1;
   else 
       dserror("Reading of FLUID2 element failed: NORM_P");       
}
else
   dserror("Reading of FLUID2 element failed: NORM_P"); 
/*----------------------------------------------- read stab flag ninths */
frchar("NINTHS",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"at_center",9)==0)
      ele->e.f2->ninths=1;
   else if (strncmp(buffer,"every_intpt",11)==0)
      ele->e.f2->ninths=2;
   else
      dserror("Reading of FLUID2 element failed: NINTHS"); 
}
else 
  dserror("Reading of FLUID2 element failed: NINTHS");
/*----------------------------------------------- read stab flag istapc */
frchar("ISTAPC",buffer,&ierr);
if (ierr==1)
{
   if (strncmp(buffer,"at_center",9)==0)
      ele->e.f2->istapc=1;
   else if (strncmp(buffer,"every_intpt",11)==0)
      ele->e.f2->istapc=2;
   else
      dserror("Reading of FLUID2 element failed: ISTAPC"); 
}
else 
  dserror("Reading of FLUID2 element failed: ISTAPC");      
/*----------------------------------------------- read stab flag istapa */
frint("ISTAPA",&(ele->e.f2->istapa),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: ISTAPA");
/*----------------------------------------------- read stab flag istapc */
frint("ISTAPC",&(ele->e.f2->istapc),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: ISTAPC");
/*--------------------------------------------------- read stab flag mk */
frint("MK",&(ele->e.f2->mk),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: MK");
/*----------------------------------------------- read stab flag ihelem */
frint("IHELEM",&ihelem,&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: IHELEM");
/*---------------------------------------------- read stab const c_lamb */
frdouble("C_LAMB",&(ele->e.f2->clamb),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: C_LAMB");
/*-------------------------------------------------- read stab const c1 */
frdouble("C1",&(ele->e.f2->c1),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: C1");
/*-------------------------------------------------- read stab const c2 */
frdouble("C2",&(ele->e.f2->c2),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: C2");
/*-------------------------------------------------- read stab const c3 */
frdouble("C3",&(ele->e.f2->c3),&ierr); 
if (ierr!=1) dserror("Reading of FLUID2 element failed: C3");
/*----------------------------- her no read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */
/*----------------------------- her now read all the other stuff needed */

/*-------------------------------------------- set some additional data */
/* initialisation */
ele->e.f2->istrle  = 0;
ele->e.f2->iarea   = 0;
ele->e.f2->iduring = 0;
ele->e.f2->idiaxy  = 0;
itaumu  = 0;
itaump  = 0;
itauc	= 0;

intextract(ihelem,&ndum,
           &(ele->e.f2->ihele[0]),&(ele->e.f2->ihele[1]),&(ele->e.f2->ihele[2]));

for(i=0;i<3;i++)
{
   if (ele->e.f2->ihele[i]==5)
   {
      ele->e.f2->istrle  = 1;
      if (ele->e.f2->iadvec!=0 && ele->e.f2->ninths==1)
         itaumu = -1;
      if (ele->e.f2->iadvec!=0 && ele->e.f2->ninths!=1)
         itaumu = 1;
      if (ele->e.f2->ipres!=0 && ele->e.f2->ninths==1)
         itaump = -1;
      if (ele->e.f2->ipres!=0 && ele->e.f2->ninths!=1)
         itaump = 1;
      if (ele->e.f2->icont!=0 && ele->e.f2->ninths==1)
         itauc = -1;
      if (ele->e.f2->icont!=0 && ele->e.f2->ninths!=1)
         itauc = 1;     
   }
   else if (ele->e.f2->ihele[i]!=0)
   {
      ele->e.f2->iarea = 1;
      if (ele->e.f2->iadvec!=0 && ele->e.f2->istapc==1)
         itaumu = -1;
      if (ele->e.f2->iadvec!=0 && ele->e.f2->istapc!=1)
         itaumu = 1;
      if (ele->e.f2->ipres!=0 && ele->e.f2->istapc==1)
         itaump = -1;
      if (ele->e.f2->ipres!=0 && ele->e.f2->istapc!=1)
         itaump = 1;
      if (ele->e.f2->icont!=0 && ele->e.f2->istapc==1)
         itauc = -1;
      if (ele->e.f2->icont!=0 && ele->e.f2->istapc!=1)
         itauc = 1;
	 
      if (ele->e.f2->ihele[i]==4)
         ele->e.f2->idiaxy = 1;         
   }
}

if (ele->e.f2->istrle==1 && ele->e.f2->ninths!=1)
   ele->e.f2->iduring = 1;
if (ele->e.f2->iarea==1 && ele->e.f2->istapc!=1)
   ele->e.f2->iduring = 1;

/*------------------------------------------- store data at the element */
ele->e.f2->itau[0] = itaumu;
ele->e.f2->itau[1] = itaump;
ele->e.f2->itau[2] = itauc;
   
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2inp */
