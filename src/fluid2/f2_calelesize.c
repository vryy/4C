#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#define ZERO  (0.0)
#define ONE   (1.0)
#define TWO   (2.0)
#define THREE (3.0)
#define FOUR  (4.0)
#define EIGHT (8.0)

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | routine to calculate elment size and stabilisation parameter         |
 | (at center) for one element                             genk 04/02   |
 |                                                                      |
 | ele->e.f2->iadvec: adevction stab.					|
 |    0 = no								|
 |    1 = yes								|
 | ele->e.f2->ipres: pressure stab.					|
 |    0 = no								|
 |    1 = yes								|
 | ele->e.f2->ivisc: diffusion stab.					|
 |    0 = no								|
 |    1 = GLS-								|
 |    2 = GLS+								|
 | ele->e.f2->icont: continuity stab.					|
 |    0 = no								|
 |    1 = yes								|
 | ele->e.f2->istapa: version of stab. parameter			|
 |    35 = diss wall instationary					|
 |    36 = diss wall stationanary					|
 | ele->e.f2->norm_P: p-norm						|
 |    p = 1<=p<=oo							|
 |    0 = max.-norm (p=oo)						|
 | ele->e.f2->mk: higher order elements control flag                    |
 |    0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)		|
 |    1 = min(1/3,2*C)							|
 |   -1 = mk=1/3  (--> element order via approx. nodal-distance)	|
 | ele->e.f2->ihele[]:                                                  |
 |    x/y/z = length-def. for velocity/pressure/continuity stab         |
 |    0 = don't compute 						|
 |    1 = sqrt(area)							|
 |    2 = area equivalent diameter					|
 |    3 = diameter/sqrt(2)						|
 |    4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle) 	|
 |    5 = streamlength (element length in flow direction		|
 | ele->e.f2->ninths: number of integration points for streamlength	|
 |    1 = at center of element						|
 |    2 = at every int pt used for element.-stab.-matrices		|
 | ele->e.f2->istapc: flag for stabilisation parameter calculation	|
 |    1 = at center of element						|
 |    2 = at every integration point					|
 | ele->e.f2->clamb \							|
 | ele->e.f2->c1     |_>> stabilisation constants (input)		|
 | ele->e.f2->c2     |							|
 | ele->e.f2->c3    /							|
 | ele->e.f2->istrle: has streamlength to be computed			|
 | ele->e.f2->iarea: calculation of area length 			|
 | ele->e.f2->iduring: calculation during int.-pt.loop			|
 | ele->e.f2->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	|
 | ele->e.f2->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	|
 | ele->e.f2->itau[2]: flag for tau_c calc. (-1: before, 1:during) 	|
 | ele->e.f2->hk[i]: "element sizes" (vel / pre / cont)                 | 
 | ele->e.f2->idiaxy: has diagonals to be computed			|
 | dynvar->tau[0]: stability parameter momentum / velocity (tau_mu)	|
 | dynvar->tau[1]: stability parameter momentum / pressure (tau_mp)	|
 | dynvar->tau[2]: stability parameter continuity (tau_c)	        |
 |									|
 *----------------------------------------------------------------------*/
void f2_calelesize(
	           ELEMENT         *ele,
		   F2_DATA         *data,
		   FLUID_DYN_CALC  *dynvar,
	           double          *funct,
	           double         **deriv,
	           double         **deriv2,	       
		   double         **xjm,
		   double         **evel,	       
		   double          *velint,
		   double         **cutp		
                  )
{
int ieval = 0;
int igc = 0;
int i,ilen;
int istrnint;
int isharea;
int ntyp;
int actmat;
int iel;
double visc;
double area, det;
double strle;
double e1,e2;
double fac,facr,facs;
double dia,dia1,dia2;
double dx,dy;
double gcoor[2];
DIS_TYP typ;

#ifdef DEBUG 
dstrc_enter("f2_calelesize");
#endif		

/*---------------------------------------------------------- initialise */
ntyp   = ele->e.f2->ntyp;
iel    = ele->numnp;
typ    = ele->distyp;

istrnint = ele->e.f2->istrle * ele->e.f2->ninths;
isharea  = dynvar->ishape * ele->e.f2->iarea;

/*----------------------------------------------------------------------*
 | calculations at element center: area & streamlength                  |
 | NOTE:                                                                |
 |    area is always calculated using only 1 integrationpoint           |     
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,facr,facs with their constant values in the calls of   |
 |         f2_rec / f2_tri!!!!!!                                        |
 *----------------------------------------------------------------------*/

if (isharea==1)
{
   area  = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(ntyp)
   {
   case 1:    /* --> quad - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,2);
      break;
   case 2:       /* --> tri - element */              
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      f2_tri(funct,deriv,deriv2,e1,e2,typ,2);
      break;
   default:
      dserror("ntyp unknown!");      
   }
   ieval++;
/* ------------------------------------------- compute jacobian matrix */      
   f2_jaco(funct,deriv,xjm,&det,ele,iel);
   fac=facr*facs*det;
   area += fac;
   if (istrnint==1)    /* compute streamlength */
   {
      f2_veli(velint,funct,evel,iel);
      ieval++;
      f2_gcoor(funct,ele,iel,gcoor);
      igc++;
      f2_calstrlen(&strle,velint,ele,gcoor,cutp,ntyp);            
   }
   if (ele->e.f2->idiaxy==1)    /* compute diagonal based diameter */
   {
      switch(ntyp)
      {
      case 1:
         dx = ele->node[1]->x[0] - ele->node[3]->x[0];
	 dy = ele->node[1]->x[1] - ele->node[3]->x[1];
	 dia1 = sqrt(dx*dx+dy*dy);
	 dx = ele->node[2]->x[0] - ele->node[4]->x[0];
	 dy = ele->node[2]->x[1] - ele->node[4]->x[1];
	 dia2 = sqrt(dx*dx+dy*dy);
/*------ dia=sqrt(2)*area/(1/2*(dia1+dia2))=sqrt(8)*area/(dia1+dia2) ---*/
	 dia = sqrt(EIGHT)*area/(dia1+dia2); 
         break;
      case 2:    /* get global coordinate of element center */
         if (igc==0)
	    f2_gcoor(funct,ele,iel,&gcoor);
	 dia = ZERO;
	 for (i=0;i<3;i++)
	 {
	    dx = gcoor[0] - ele->node[0]->x[i];
	    dy = gcoor[1] - ele->node[2]->x[i];
	    dia += dx*dx + dy*dy;
	 }
	 dia = FOUR*area/sqrt(THREE*dia);
         break;
      default:
          dserror("ntyp unknown!");
      }
   }
/*--------------------------------------------------- set element sizes *
  ----loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for(ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f2->ihele[ilen]==1)
         ele->e.f2->hk[ilen] = sqrt(area);
      else if (ele->e.f2->ihele[ilen]==2)
         ele->e.f2->hk[ilen] = TWO*sqrt(area/PI);
      else if (ele->e.f2->ihele[ilen]==3)
         ele->e.f2->hk[ilen] = sqrt(TWO*area/PI);
      else if (ele->e.f2->ihele[ilen]==4)
         ele->e.f2->hk[ilen] = dia;
      else if (ele->e.f2->ninths==1)
         ele->e.f2->hk[ilen] = strle;  
   }
}   
/*----------------------------------------------------------------------*
 | calculations at element center: only streamlength                    |
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,facr,facs with their constant values in the calls of   |
 |         f2_rec / f2_tri!!!!!!                                        |
 *----------------------------------------------------------------------*/
else if (istrnint==1 && isharea !=1) 
{
   area  = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(ntyp)
   {
   case 1:    /* --> quad - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,2);
      break;
   case 2:       /* --> tri - element */              
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      f2_tri(funct,deriv,deriv2,e1,e2,typ,2);
      break;
   default:
      dserror("ntyp unknown!");
   }
   ieval++;
/* ------------------------------------------- compute jacobian matrix */      
   f2_jaco(funct,deriv,xjm,&det,ele,iel);
/*----------------------------------------------- compute streamlength */
   f2_veli(velint,funct,evel,iel);
   ieval++;
   f2_gcoor(funct,ele,iel,gcoor);
   igc++;
   f2_calstrlen(&strle,velint,ele,gcoor,cutp,ntyp);       
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f2->ihele[ilen]==5)
         ele->e.f2->hk[ilen] = strle;   
   }
}

/*----------------------------------------------------------------------*
  calculate stabilisation parameter
 *----------------------------------------------------------------------*/
if(ele->e.f2->istapc==1 || istrnint==1)
{
   switch(ieval) /* ival>2: vel at intpoint already available! */
   {
   case 0:
/*------ get only values of integration parameters and shape functions
        no derivatives -------------------------------------------------*/
      switch(ntyp)
      {
      case 1:    /* --> quad - element */
         e1   = data->qxg[0][0];
         facr = data->qwgt[0][0];
         e2   = data->qxg[0][0];
         facs = data->qwgt[0][0];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,1);
         break;
      case 2:       /* --> tri - element */              
         e1   = data->txgr[0][0];
         facr = data->twgt[0][0];
         e2   = data->txgs[0][0];
         facs = ONE;
         f2_tri(funct,deriv,deriv2,e1,e2,typ,1);
         break;      
      default:
         dserror("ntyp unknown!");
      }
      f2_veli(velint,funct,evel,iel);
      break;
   case 1:            
      f2_veli(velint,funct,evel,iel);
      break;
   case 2:
      break;
   default:
      dserror("wrong value for ieval");
   }
/*----------------------------------- calculate stabilisation parameter */               
   actmat=ele->mat-1;
   visc = mat[actmat].m.fluid->viscosity;
   f2_calstabpar(ele,dynvar,velint,visc,iel,ntyp,-1); 
}


#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calelesize */

/*----------------------------------------------------------------------*
 | routine to calculate elment size and stabilisation parameter     |
 | (at center) for one element during integration loop     genk 04/02   |
 *----------------------------------------------------------------------*/
 void f2_calelesize2(
	             ELEMENT         *ele,
	  	     FLUID_DYN_CALC  *dynvar,
	             double          *funct,	       	       
		     double          *velint,
		     double         **cutp,
		     double           visc,
		     int              iel,
		     int              ntyp	
                    )
{
int    ilen;
int    istrnint;
double strle;
double val;
double gcoor[2];

#ifdef DEBUG 
dstrc_enter("f2_calelesize2");
#endif

/*---------------------------------------------------------- initialise */
istrnint = ele->e.f2->istrle * ele->e.f2->ninths;
val = ZERO;

if (istrnint==2)
{
/*------------------------------------------------ compute streamlength */
   f2_gcoor(funct,ele,iel,gcoor);
   f2_calstrlen(&strle,velint,ele,gcoor,cutp,ntyp);
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f2->ihele[ilen]==5)
         ele->e.f2->hk[ilen] = strle;   
   }
}
/*----------------------------------- calculate stabilisation parameter */               
   f2_calstabpar(ele,dynvar,velint,visc,iel,ntyp,1); 

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calelesize2 */
	
/*----------------------------------------------------------------------*
 | routine to calculate streamlength                                    |
 | (only approxmation for higher order elements - boundaries are        |
 | assumed to be straight)                                   genk 04/02 |
 *----------------------------------------------------------------------*/
void f2_calstrlen(
                   double   *strle,    /* streamlength  */
		   double   *velint,   /* velocities at integr. point */
		   ELEMENT  *ele,      /* actual element */
                   double   *gcoor,    /* global coordinates of int. point */
		   double  **cutp,     /* cutting points */       
		   int       ntyp      /* flag for element type */
                 )
{
int nodcut=-1;
int nodmax;
int inod;
double dl,dx,dy,dxh,dyh;
double dsub,dval;

#ifdef DEBUG 
dstrc_enter("f2_calstrlen");
#endif

dval = FABS(velint[0])+FABS(velint[1]);
if (dval == ZERO)  /* no flow at this point - take some arbitr. measure for streamlength */
{
   dx = ele->node[1]->x[0] - ele->node[3]->x[0];
   dy = ele->node[1]->x[1] - ele->node[3]->x[1];
   goto calc2;   
}

/*----------------------------------------------------------------------*
   streamlength is calculated via cutting points of velocity vector
   with straight boundaries                                             */
switch(ntyp)
{
case 1: /* max number of nodes for quad: 4 --> C-numbering nodmax = 4-1 */
   nodmax = 3;
   break;
case 2:  /* max number of nodes for tri: 3 --> C-numbering nodmax = 3-1 */
   nodmax = 2;
   break;
default:
   dserror("ntyp unknown!");   
}        
 /*------------------------------------------------- get cutting points */
for (inod=0;inod<nodmax;inod++)
{
   dxh = ele->node[inod+1]->x[0] - ele->node[inod]->x[0];
   dyh = ele->node[inod+1]->x[1] - ele->node[inod]->x[1];
   dsub = dxh*velint[1]-dyh*velint[0];
   if (dsub==ZERO)  /* check for parallel vectors */
      continue;
   dl = ((ele->node[inod]->x[1]-gcoor[1])*velint[0] -	\
	 (ele->node[inod]->x[0]-gcoor[0])*velint[1])/dsub;
   if (dl>=ZERO && dl<=ONE)
   {
      nodcut++;
      cutp[0][nodcut]=ele->node[inod]->x[0]+dl*dxh;
      cutp[1][nodcut]=ele->node[inod]->x[1]+dl*dyh;
      if (nodcut==1)
	 goto calc1;
   } 
}
/*------------------------------------------------- the last boundary */
dxh = ele->node[0]->x[0]-ele->node[nodmax]->x[0];
dyh = ele->node[0]->x[1]-ele->node[nodmax]->x[1];

dsub = dxh*velint[1] - dyh*velint[0];
if (dsub==ZERO)
   dserror("Couldn't find two cutting points!");
dl = ((ele->node[nodmax]->x[1]-gcoor[1])*velint[0] -	\
      (ele->node[nodmax]->x[0]-gcoor[0])*velint[1])/dsub;
if (dl>=ZERO && dl <= ONE)
{
   nodcut++;
   cutp[0][nodcut]=ele->node[nodmax]->x[0]+dl*dxh;
   cutp[1][nodcut]=ele->node[nodmax]->x[1]+dl*dyh; 
   if(nodcut==1)
      goto calc1;
}      

dserror("Couldn't find two cutting points!");

calc1:
dx = cutp[0][1]-cutp[0][0];
dy = cutp[1][1]-cutp[1][0];

calc2:
*strle = sqrt(dx*dx+dy*dy);
		
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_calstrlen */		  
