#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#define Q12  (0.5)
#define Q14  (0.25)
#define Q16  (0.1666666666667)
#define ZERO (0.0)
#define ONE  (1.0)
#define TWO  (2.0)

/*----------------------------------------------------------------------*
 | routine for evaluation of shpae functions and their natural first    |
 | and secend derivatives with respect to r/s for                       |
 | R E C T A N G L E S                                     genk 04/02   |
 *----------------------------------------------------------------------*/
void f2_rec(
            double     *funct,    /* shape functions */ 
            double    **deriv,    /* first natural deriv. of shape funct. */
            double    **deriv2,   /* second natural deriv. of shape funct. */
	    double      r,        
            double      s,        /* coordinates */
            DIS_TYP     typ,      /* element type */
            int         icode     /* evaluation flag */
	   )
{
int    i,ii;
double rp,rm,sp,sm,r2,s2;
double rh,sh,rs;
double rhm,shm,rhp,shp;

#ifdef DEBUG 
dstrc_enter("f2_rec");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
case quad4: /* LINEAR shape functions and their natural derivatives ----*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   
   funct[0]=Q14*rp*sp;
   funct[1]=Q14*rm*sp;
   funct[2]=Q14*rm*sm;
   funct[3]=Q14*rp*sm;

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= Q14*sp;
      deriv[1][0]= Q14*rp;
      deriv[0][1]=-Q14*sp;
      deriv[1][1]= Q14*rm;
      deriv[0][2]=-Q14*sm;
      deriv[1][2]=-Q14*rm;
      deriv[0][3]= Q14*sm;
      deriv[1][3]=-Q14*rp;
   }

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= ZERO;
      deriv2[1][0]= ZERO;
      deriv2[2][0]= Q14;
      deriv2[0][1]= ZERO;
      deriv2[1][1]= ZERO;
      deriv2[2][1]=-Q14;
      deriv2[0][2]= ZERO;
      deriv2[1][2]= ZERO;
      deriv2[2][2]= Q14;
      deriv2[0][3]= ZERO;
      deriv2[1][3]=ZERO;
      deriv2[2][3]=-Q14;      
   }
   break;

/*----------------------------------------------------------------------*/
case quad8: /* QUADRATIC shape functions and their natural derivatives
               without central node ------------------------------------*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   r2=ONE-r*r;
   s2=ONE-s*s;
   
   funct[0]=Q14*rp*sp;
   funct[1]=Q14*rm*sp;
   funct[2]=Q14*rm*sm;
   funct[3]=Q14*rp*sm;
   funct[4]=Q12*r2*sp;
   funct[5]=Q12*rm*s2;
   funct[6]=Q12*r2*sm;
   funct[7]=Q12*rm*s2;
   funct[0]=funct[0] - Q12*(funct[4]+funct[7]);   

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= Q14*sp;
      deriv[1][0]= Q14*rp;
      deriv[0][1]=-Q14*sp;
      deriv[1][1]= Q14*rm;
      deriv[0][2]=-Q14*sm;
      deriv[1][2]=-Q14*rm;
      deriv[0][3]= Q14*sm;
      deriv[1][3]=-Q14*rp;
      deriv[0][4]=-ONE*r*sp;
      deriv[1][4]= Q12*r2;
      deriv[0][5]=-Q12*  s2;
      deriv[1][5]=-ONE*rm*s;
      deriv[0][6]=-ONE*r*sm;
      deriv[1][6]=-Q12*r2;
      deriv[0][7]= Q12*  s2;
      deriv[1][7]=-ONE*rp*s;

      deriv[0][0]=deriv[0][0] - Q12*(deriv[0][4]+deriv[0][7]);
      deriv[1][0]=deriv[1][0] - Q12*(deriv[1][4]+deriv[1][7]);      
   }

   for(i=2;i<=4;i++)
   {
      ii=i+3;
      funct[i]=funct[i] - Q12*(funct[ii]+funct[ii+1]);
      if (icode>1)
      {
         deriv[0][i]=deriv[0][i] - Q12*(deriv[0][ii]+deriv[0][ii+1]);
	 deriv[1][i]=deriv[1][i] - Q12*(deriv[1][ii]+deriv[1][ii+1]);
      }
   }

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= ZERO;
      deriv2[1][0]= ZERO;
      deriv2[2][0]= Q14;
      deriv2[0][1]= ZERO;
      deriv2[1][1]= ZERO;
      deriv2[2][1]=-Q14;
      deriv2[0][2]= ZERO;
      deriv2[1][2]= ZERO;
      deriv2[2][2]= Q14;
      deriv2[0][3]= ZERO;
      deriv2[1][3]= ZERO;
      deriv2[2][3]=-Q14;
      deriv2[0][4]=-(ONE+s);
      deriv2[1][4]= ZERO;
      deriv2[2][4]=-r;
      deriv2[0][5]= ZERO;
      deriv2[1][5]=-(ONE-r);
      deriv2[2][5]= s;
      deriv2[0][6]=-(ONE-s);
      deriv2[1][6]= ZERO;
      deriv2[2][6]= r;
      deriv2[0][7]= ZERO;
      deriv2[1][7]=-(ONE+r);
      deriv2[2][7]=-s;
      
      deriv2[0][0]=deriv2[0][0] - Q12*(deriv2[0][4]+deriv2[0][7]);
      deriv2[1][0]=deriv2[1][0] - Q12*(deriv2[1][4]+deriv2[1][7]);
      deriv2[2][0]=deriv2[2][0] - Q12*(deriv2[2][4]+deriv2[2][7]); 
      
      for(i=2;i<=4;i++)
      {
         ii=i+3;
         deriv2[0][i]=deriv2[0][i] - Q12*(deriv2[0][ii]+deriv2[0][ii+1]);
	 deriv2[1][i]=deriv2[1][i] - Q12*(deriv2[1][ii]+deriv2[1][ii+1]);
	 deriv2[2][i]=deriv2[2][i] - Q12*(deriv2[2][ii]+deriv2[2][ii+1]);
      }           
   }
   break;

/*----------------------------------------------------------------------*/
case quad9: /* QUADRATIC  shape functions and their natural derivatives
               with central node ---------------------------------------*/
/*--------------------------------------------------- form basic values */
   rp=ONE+r;
   rm=ONE-r;
   sp=ONE+s;
   sm=ONE-s;
   r2=ONE-r*r;
   s2=ONE-s*s;
   rh=Q12*r;
   sh=Q12*s;
   rs=rh*sh;
   rhp=r+Q12;
   rhm=r-Q12;
   shp=s+Q12;
   shm=s-Q12;
   
   funct[0]= rs*rp*sp;
   funct[1]=-rs*rm*sp;
   funct[2]= rs*rm*sm;
   funct[3]=-rs*rp*sm;
   funct[4]= sh*sp*r2;
   funct[5]=-rh*rm*s2;
   funct[6]=-sh*sm*r2;
   funct[7]= rh*rp*s2;
   funct[8]= r2*s2;   	       
      
   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]= rhp*sh*sp;
      deriv[1][0]= shp*rh*rp;
      deriv[0][1]= rhm*sh*sp;
      deriv[1][1]=-shp*rh*rm;
      deriv[0][2]=-rhm*sh*sm;
      deriv[1][2]=-shm*rh*rm;
      deriv[0][3]=-rhp*sh*sm;
      deriv[1][3]= shm*rh*rp;
      deriv[0][4]=-TWO*r*sh*sp;
      deriv[1][4]= shp*r2;
      deriv[0][5]= rhm*s2;
      deriv[1][5]= TWO*s*rh*rm;
      deriv[0][6]= TWO*r*sh*sm;
      deriv[1][6]= shm*r2;
      deriv[0][7]= rhp*s2;
      deriv[1][7]=-TWO*s*rh*rp; 
      deriv[0][8]=-TWO*r*s2;
      deriv[1][8]=-TWO*s*r2;         
   }
   
   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= sh*sp;
      deriv2[1][0]= rh*rp;
      deriv2[2][0]= shp*rhp;
      deriv2[0][1]= sh*sp;
      deriv2[1][1]=-rh*rm;
      deriv2[2][1]= shp*rhm;
      deriv2[0][2]=-sh*sm;
      deriv2[1][2]=-rh*rm;
      deriv2[2][2]= shm*rhm;
      deriv2[0][3]=-sh*sm;
      deriv2[1][3]= rh*rp;
      deriv2[2][3]= shm*rhp;
      deriv2[0][4]=-TWO*sh*sp;
      deriv2[1][4]= r2;
      deriv2[2][4]=-TWO*r*shp;
      deriv2[0][5]= s2;
      deriv2[1][5]= TWO*rh*rm;
      deriv2[2][5]=-TWO*s*rhm;
      deriv2[0][6]= TWO*sh*sm;
      deriv2[1][6]= r2;
      deriv2[2][6]=-TWO*r*shm;
      deriv2[0][7]= s2;
      deriv2[1][7]=-TWO*rh*rp;
      deriv2[2][7]=-TWO*s*rhp;
      deriv2[0][8]=-TWO*s2;
      deriv2[1][8]=-TWO*r2;
      deriv2[2][8]= TWO*s*TWO*r;         
   }
   break;

/*----------------------------------------------------------------------*/   
default:
   dserror("distyp unknown");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_rec */

/*----------------------------------------------------------------------*
 | routine for evaluation of shpae functions and their natural first    |
 | and secend derivatives with respect to r/s for                       |
 | T R I A N G L E S                                       genk 06/02   |
 *----------------------------------------------------------------------*/
void f2_tri(
            double     *funct,    /* shape functions */ 
            double    **deriv,    /* first natural deriv. of shape funct. */
            double    **deriv2,   /* second natural deriv. of shape funct. */
	    double      r, 
            double      s,	  /* coordinates */
            DIS_TYP     typ,	  /* element type */
            int         icode	  /* evaluation flag */
	   )
{

#ifdef DEBUG 
dstrc_enter("f2_tri");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
case tri3:
   dserror("shape functions for TRI3 not implemented yet!");
   break;

case tri6:
   dserror("shape functions for TRI6 not implemented yet!"); 
default:
   dserror("distyp unknown");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_rec */

/*----------------------------------------------------------------------*
 | routine to calculate the globael derivatives wrt x,y at point r,s    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_gder(
              double   **derxy,    /* global derivatives wrt. x/y       */
	      double   **deriv,    /* derivatives of shape functions    */
	      double   **xjm,      /* jacobian matrix                   */
	      double     det,      /* jacobian determinant              */
	      int        iel       /* number of nodes in actual element */
	    )
{
int    k;

#ifdef DEBUG 
dstrc_enter("f2_tri");
#endif

/*------------------------------------------------------- initialistion */
for(k=0;k<iel;k++)
{
   derxy[0][k]=ZERO;
   derxy[1][k]=ZERO;
}

/* NOTE: calculation of inversion is skipped 
         --> change of subscripts and signs in loop!! */

for (k=0;k<iel;k++)
{
   derxy[0][k] +=(  xjm[1][1]*deriv[0][k]  + (-xjm[0][1]*deriv[1][k]))/det;
   derxy[1][k] +=((-xjm[1][0]*deriv[0][k]) +   xjm[0][0]*deriv[1][k] )/det;
} 	  

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_gder */

/*----------------------------------------------------------------------*
 | routine to get global coordinates for given shape function values    |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_gcoor(
              double     *funct,    /* shape functions */ 
              ELEMENT    *ele,
	      int         iel,      /* number of nodes in actual element */
	      double     *gcoor     /* global coordinates */
	     )
{
int i;

#ifdef DEBUG 
dstrc_enter("f2_gder2");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;

for(i=0;i<iel;i++)
{
   gcoor[0] += funct[i] * ele->node[i]->x[0];
   gcoor[1] += funct[i] * ele->node[i]->x[1];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_gcoor */

/*----------------------------------------------------------------------*
 | routine to calculate second global derivatives wrt x/y at point r,s  |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void f2_gder2(
               ELEMENT     *ele,
	       double     **xjm,            
               double     **bi,
	       double     **xder2,
	       double     **derxy,
	       double     **derxy2,
               double     **deriv2,
	       int          iel
	     )
{
int i;
double x00,x01,x02,x10,x11,x12,x20,x21,x22;
double det,dum;
double r0,r1,r2;

#ifdef DEBUG 
dstrc_enter("f2_gder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
x00 = xjm[0][0]*xjm[0][0];
x01 = xjm[0][1]*xjm[0][1];
x02 = TWO*xjm[0][0]*xjm[0][1];
x10 = xjm[1][0]*xjm[1][0];
x11 = xjm[1][1]*xjm[1][1];
x12 = TWO*xjm[1][0]*xjm[1][1];
x20 = xjm[0][0]*xjm[1][0];
x21 = xjm[0][1]*xjm[1][1];
x22 = xjm[0][0]*xjm[1][1] + xjm[0][1]*xjm[1][0];

/*-------------------------------------- inverse of jacobian_bar matrix */
det =   x00*x11*x22 + x01*x12*x21 + x10*x21*x02 \
      - x20*x11*x02 - x00*x12*x21 - x01*x10*x22 ;
dum = ONE/det;         
bi[0][0] =   dum*(x11*x22 - x12*x21);
bi[1][0] =  -dum*(x10*x22 - x20*x12);
bi[2][0] =   dum*(x10*x21 - x20*x11);
bi[0][1] =  -dum*(x01*x22 - x21*x02);
bi[1][1] =   dum*(x00*x22 - x02*x20);
bi[2][1] =  -dum*(x00*x21 - x20*x01);
bi[0][2] =   dum*(x01*x12 - x11*x02);
bi[1][2] =  -dum*(x00*x12 - x10*x02);
bi[2][2] =   dum*(x00*x11 - x01*x10);

/*---------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   xder2[i][0]=ZERO;
   xder2[i][1]=ZERO;
}
for (i=0;i<iel;i++)
{
   derxy2[0][i]=ZERO;
   derxy2[1][i]=ZERO;
   derxy2[2][i]=ZERO;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++)
{
   xder2[0][0] = xder2[0][0] + deriv2[0][i] * ele->node[i]->x[0];
   xder2[1][0] = xder2[1][0] + deriv2[1][i] * ele->node[i]->x[0];
   xder2[2][0] = xder2[2][0] + deriv2[2][i] * ele->node[i]->x[0];
   xder2[0][1] = xder2[0][1] + deriv2[0][i] * ele->node[i]->x[1];
   xder2[1][1] = xder2[1][1] + deriv2[1][i] * ele->node[i]->x[1];
   xder2[2][1] = xder2[2][1] + deriv2[2][i] * ele->node[i]->x[1];
}

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++)
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i];
   
   derxy2[0][i] += bi[0][0]*r0 + bi[0][1]*r1 + bi[0][2]*r2;
   derxy2[1][i] += bi[1][0]*r0 + bi[1][1]*r1 + bi[1][2]*r2;
   derxy2[2][i] += bi[2][0]*r0 + bi[2][1]*r1 + bi[2][2]*r2;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_gcoor */
