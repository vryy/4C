/*!----------------------------------------------------------------------
\file
\brief shape functions and their natural derivatives for fluid2 element

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
static DOUBLE Q12 = ONE/TWO;
static DOUBLE Q14 = ONE/FOUR;
static DOUBLE Q16 = ONE/SIX;
static DOUBLE Q52 = FIVE/TWO;
static DOUBLE C23 = FIVE*FOUR+THREE;
/*----------------------------------------------------------- prototype */
void fluid_mf(INT mctrl);
void dyn_fsi(INT mctrl);
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for rectangles

<pre>                                                         genk 04/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s are evaluated for 
 R E C T A N G L E S 
		     
Numbering of the nodes:


                    ^ s
                    |
              1     |4    0
              o-----o-----o
              |           |
              |           |7
            5 o     o     o -----> r
              |     8     |
              |           |
              o-----o-----o
              2     6     3

</pre>
\param  *funct     DOUBLE   (o)    shape functions
\param **deriv     DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2    DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r 	   DOUBLE   (i)    coordinate
\param   s 	   DOUBLE   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   INT	    (i)    evaluation flag
\return void                                                                       

------------------------------------------------------------------------*/
void f2_rec(
            DOUBLE     *funct,     
            DOUBLE    **deriv,    
            DOUBLE    **deriv2,   
	    DOUBLE      r,        
            DOUBLE      s,        
            DIS_TYP     typ,      
            INT         icode     
	   )
{
INT    i,ii;
DOUBLE rp,rm,sp,sm,r2,s2;
DOUBLE rh,sh,rs;
DOUBLE rhm,shm,rhp,shp;

#ifdef DEBUG 
dstrc_enter("f2_rec");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
/*----------------------------------------------------------------------*/
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
   } /* endif (icode>1) */

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
   } /* endif (icode==3) */
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
   
   funct[4]=Q12*r2*sp;
   funct[5]=Q12*rm*s2;
   funct[6]=Q12*r2*sm;
   funct[7]=Q12*rp*s2;
   funct[0]=Q14*rp*sp-Q12*(funct[4]+funct[7]);
   funct[1]=Q14*rm*sp-Q12*(funct[4]+funct[5]); 
   funct[2]=Q14*rm*sm-Q12*(funct[5]+funct[6]); 
   funct[3]=Q14*rp*sm-Q12*(funct[6]+funct[7]); 
   
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

      deriv[0][0]-= Q12*(deriv[0][4]+deriv[0][7]);
      deriv[1][0]-= Q12*(deriv[1][4]+deriv[1][7]);      

      for(i=1;i<4;i++)
      {
         ii=i+3;
         deriv[0][i] -= Q12*(deriv[0][ii]+deriv[0][ii+1]);
	 deriv[1][i] -= Q12*(deriv[1][ii]+deriv[1][ii+1]);     
      } /* end loop over i */
   } /* endif (icode>1) */

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
      
      deriv2[0][0] -= Q12*(deriv2[0][4]+deriv2[0][7]);
      deriv2[1][0] -= Q12*(deriv2[1][4]+deriv2[1][7]);
      deriv2[2][0] -= Q12*(deriv2[2][4]+deriv2[2][7]); 
      
      for(i=1;i<4;i++)
      {
         ii=i+3;
         deriv2[0][i] -= Q12*(deriv2[0][ii]+deriv2[0][ii+1]);
	 deriv2[1][i] -= Q12*(deriv2[1][ii]+deriv2[1][ii+1]);
	 deriv2[2][i] -= Q12*(deriv2[2][ii]+deriv2[2][ii+1]);
      } /* end loop over i */
   } /* endif (icode==3) */
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
   } /* endif (icode>1) */
   
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
   } /* endif (icode==3) */
break;

/*----------------------------------------------------------------------*/   
default:
   dserror("distyp unknown");
} /* end switch(typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_rec */

/*!--------------------------------------------------------------------- 
\brief extroplation for  rectangles

<pre>                                                         genk 10/02

In this routine function values at the gauss points are extrapolated to 
the element nodes.
   icode = 0: only evaluation of function parameters
   icode = 1: evaluationn of function parameters and extrapolation
   icode = 2: only extrapolation
		     
</pre>
\param  *funval      DOUBLE   (o)    function value at the element nodes 
\param  *fpar        DOUBLE   (i/o)  function parameters
\param   r           DOUBLE   (i)    local coord. in r direction 
\param   s           DOUBLE   (i)    local coord. in s direction 
\param   fval        DOUBLE   (i)    function values at gauss points
\param   igauss      INT      (i)    number of gauss points
\param   icode       INT      (i)    evaluation flag
\return  void                                                                      
------------------------------------------------------------------------*/
void f2_recex(
	      DOUBLE           *funval,     
	      DOUBLE           *fpar,
	      DOUBLE            r,    
	      DOUBLE            s,
              DOUBLE           *fval,         
	      INT               igauss,
	      INT               icode
            )
{
DOUBLE p1,p2,p3,p4,p5,p6,p7,p8,p9;
DOUBLE f1,f2,f3,f4,f5;

#ifdef DEBUG 
dstrc_enter("f2recex");
#endif

/*------------------------evaluate function parameters-------------------*/
if (icode==0 || icode==1)
{
   if (igauss >= 1) 
      p1 = fval[0];
   if (igauss >= 4) 
   {
      p2 = fval[1];
      p3 = fval[2];
      p4 = fval[3];
   }/*--end of if(igauss >=4)--*/
   if (igauss==8) dserror("extrapolation for 8-noded quad not implemented yet!\n");
   if (igauss >= 9)
   {
      p5 = fval[4];
      p6 = fval[5];
      p7 = fval[6];
      p8 = fval[7];
      p9 = fval[8];
   }/*--end of if(igauss >=9)--*/

   if (igauss == 1)
      fpar[0] = p1;
   else if (igauss == 4)
   {
      f1 = ONE/FOUR;
      f2 = f1*sqrt(THREE);
      f3 = f1*THREE;
      fpar[0] = f2 * (-p1+p3-p2+p4);
      fpar[1] = f3 * ( p1-p3-p2+p4);
      fpar[2] = f2 * (-p1-p3+p2+p4);
      fpar[3] = f1 * ( p1+p3+p2+p4);
   }
   else if (igauss == 8)
      dserror("extrapolation for 8-noded quad not implemented yet!\n"); 
   else if (igauss == 9)
   {
      f3 = FIVE/SIX;
      f4 = f3/TWO;
      f1 = f3*f3;
      f2 = f1*sqrt(THREE/FIVE);
      f5 = f3*sqrt(THREE/FIVE);
      fpar[0] = f1*( p1+p7+p3+p9-TWO*(p4+p2+p8+p6)+FOUR*p5);
      fpar[1] = f2*(-p1-p7+p3+p9+TWO*(p4-p6));
      fpar[2] = f2*(-p1+p7-p3+p9+TWO*(p2-p8));
      fpar[3] = f3*( p2+p8-TWO*p5);
      fpar[4] = f4*( p1-p7-p3+p9);
      fpar[5] = f3*( p4+p6-TWO*p5);
      fpar[6] = f5*(-p2+p8);
      fpar[7] = f5*(-p4+p6);
      fpar[8] = p5;
   }/*--end of if(igauss ==1)--*/
   else
     dserror ("extropolation for number of gauss point not implemented yet!\n");
   icode++;
} /* endif (icode==0) || (icode==1) */

/*---------calculate function value f(r,s)------*/
if (icode==2)
{
   if (igauss==1)
   { 
      *funval = fpar[0];	
   }
   else if (igauss==4)
   {
      p1 = fpar[0];
      p2 = fpar[1];
      p3 = fpar[2];
      p4 = fpar[3];
      *funval = p1*r + p2*r*s + p3*s +p4;
   }
   else if (igauss==8) 
      dserror("extrapolation of 8-noded quads not implemented yet!\n");
   else if (igauss==9)
   {
      p1 = fpar[0];
      p2 = fpar[1];
      p3 = fpar[2];
      p4 = fpar[3];
      p5 = fpar[4];
      p6 = fpar[5];
      p7 = fpar[6];
      p8 = fpar[7];
      p9 = fpar[8];
      *funval = p1*r*s*r*s+p2*r*s*r+p3*r*s*s+p4*r*r+p5*r*s+p6*s*s+p7*r+p8*s+p9;
   }/*--end of if (igauss ==1)--*/ 
   else
     dserror ("extropolation for number of gauss point not implemented yet!\n");   
} /* endif (icode==2) */
	
#ifdef DEBUG 
dstrc_exit();
#endif

return ;
} /* end of f2recex */

/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for triangles

<pre>                                                         genk 06/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s are evaluated for 
T R I A N G L E S 
		     
</pre>
\param  *funct     DOUBLE   (o)    shape functions
\param **deriv     DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2    DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r 	   DOUBLE   (i)    coordinate
\param   s 	   DOUBLE   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   INT	    (i)    evaluation flag
\return void                                                                       

------------------------------------------------------------------------*/
void f2_tri(
            DOUBLE     *funct,      
            DOUBLE    **deriv,    
            DOUBLE    **deriv2,   
	    DOUBLE      r,        
            DOUBLE      s,	  
            DIS_TYP     typ,	  
            INT         icode	  
	   )
{

DOUBLE rr,rs,ss;

#ifdef DEBUG 
dstrc_enter("f2_tri");
#endif

/*------------------------------- selection of polynomial interpolation */
switch (typ)
{
/*----------------------------------------------------------------------*/
case tri3: /* LINEAR shape functions and their natural derivatives -----*/
/*----------------------------------------------------------------------*/

   funct[0]=ONE-r-s;
   funct[1]=r;
   funct[2]=s; 

   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]=-ONE;
      deriv[1][0]=-ONE;
      deriv[0][1]= ONE;
      deriv[1][1]=ZERO;
      deriv[0][2]=ZERO;
      deriv[1][2]= ONE;
   } /* endif (icode>1) */
break;

/*----------------------------------------------------------------------*/
case tri6: /* QUADRATIC shape functions and ther natural derivatives ---*/
/*--------------------------------------------------- form basic values */
   rr=r*r;
   ss=s*s;
   rs=r*s;

   funct[0]=(ONE-TWO*r-TWO*s)*(ONE-r-s);
   funct[1]=TWO*rr-r;
   funct[2]=TWO*ss-s;
   funct[3]=FOUR*(r-rr-rs);
   funct[4]=FOUR*rs;
   funct[5]=FOUR*(s-rs-ss);
   
   if(icode>1) /* --> first derivative evaluation */
   {
      deriv[0][0]=-THREE+FOUR*(r+s);
      deriv[1][0]= deriv[0][0];

      deriv[0][1]= FOUR*r-ONE;
      deriv[1][1]= ZERO;

      deriv[0][2]= ZERO;
      deriv[1][2]= FOUR*s-ONE;

      deriv[0][3]= FOUR*(ONE-TWO*r-s);
      deriv[1][3]=-FOUR*r;

      deriv[0][4]= FOUR*s;
      deriv[1][4]= FOUR*r;

      deriv[0][5]=-FOUR*s;
      deriv[1][5]= FOUR*(ONE-r-TWO*s);         
   } /* endif (icode>1) */

   if(icode==3) /* --> second derivative evaluation */
   {
      deriv2[0][0]= FOUR;
      deriv2[1][0]= FOUR;
      deriv2[2][0]= FOUR;

      deriv2[0][1]= FOUR;
      deriv2[1][1]= ZERO;
      deriv2[2][1]= ZERO;

      deriv2[0][2]= ZERO;
      deriv2[1][2]= FOUR;
      deriv2[2][2]= ZERO;

      deriv2[0][3]=-EIGHT;
      deriv2[1][3]= ZERO;
      deriv2[2][3]=-FOUR;

      deriv2[0][4]= ZERO;
      deriv2[1][4]= ZERO;
      deriv2[2][4]= FOUR;

      deriv2[0][5]= ZERO;
      deriv2[1][5]=-EIGHT;
      deriv2[2][5]=-FOUR;	  
   } /* endif (icode==3) */
break;   

/*----------------------------------------------------------------------*/
default:
   dserror("distyp unknown");
} /* end swtich(typ) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_tri */

/*!---------------------------------------------------------------------
\brief extroplation for  triangles

<pre>                                                         genk 10/02

In this routine function values at the gauss points are extrapolated to 
the element nodes.
   icode = 0: only evaluation of function parameters
   icode = 1: evaluationn of function parameters and extrapolation
   icode = 2: only extrapolation
		     
</pre>
\param  *funval      DOUBLE   (o)    function value at the element nodes 
\param  *fpar        DOUBLE   (i/o)  function parameters
\param   r           DOUBLE   (i)    local coord. in r direction 
\param   s           DOUBLE   (i)    local coord. in s direction 
\param   fval        DOUBLE   (i)    function values at gauss points
\param   igauss      INT      (i)    number of gauss points
\param   icode       INT      (i)    evaluation flag
\return  void                                                                      
------------------------------------------------------------------------*/
void f2_triex(
	      DOUBLE           *funval,     
	      DOUBLE           *fpar,
	      DOUBLE            r,    
	      DOUBLE            s,
              DOUBLE           *fval,         
	      INT               igauss,
	      INT               icode
            )
{

#ifdef DEBUG 
dstrc_enter("f2_triex");
#endif

/*------------------------evaluate function parameters-------------------*/
if (icode==0 || icode==1)
{
   if (igauss==1) 
      fpar[0] = fval[0];
   else if (igauss==3) 
   {
      fpar[0] = TWO*(fval[1]-fval[2]);
      fpar[1] = TWO*(fval[1]-fval[0]);
      fpar[2] = fval[0]-fval[1]+fval[2];
   }
   else if (igauss==4)
   {
      fpar[0] = -Q52*(fval[0]-fval[1]);
      fpar[1] = -Q52*(fval[0]-fval[2]); 
      fpar[2] = C23*fval[0]-SEVEN*(fval[1]+fval[2])+THREE*fval[3];
      fpar[2] = fpar[2]/TWELVE;            
   }/*--end of if(igauss >=9)--*/
   else if (igauss==6)
   {
      fpar[0] = 0.372483342138412E+01*fval[0]				    
   	      + 0.372483342138411E+01*fval[1]				    
   	      - 0.165973973290167E+00*fval[2]				    
   	      - 0.799630368687822E+01*fval[3]				    
   	      + 0.356305408700079E+00*fval[4]				    
   	      + 0.356305408700074E+00*fval[5];				    
      fpar[1] = 0.761564081605845E+01*fval[0]				    
   	      - 0.165973973290160E+00*fval[1]				    
   	      - 0.165973973290163E+00*fval[2]				    
   	      - 0.799630368687823E+01*fval[3]				    
   	      + 0.870891450427838E+01*fval[4]				    
   	      - 0.799630368687823E+01*fval[5];				    
      fpar[2] = 0.372483342138411E+01*fval[0]				    
   	      - 0.165973973290164E+00*fval[1]				    
   	      + 0.372483342138411E+01*fval[2]				    
   	      + 0.356305408700076E+00*fval[3]				    
   	      + 0.356305408700076E+00*fval[4]				    
   	      - 0.799630368687821E+01*fval[5];				    
      fpar[3] =-0.545993310748378E+01*fval[0]				    
   	      - 0.198973373528444E+01*fval[1]				    
   	      + 0.165973973290165E+00*fval[2]				    
   	      + 0.799630368687822E+01*fval[3]				    
   	      - 0.112120572260041E+01*fval[4]				    
   	      + 0.408594905200256E+00*fval[5];				    
      fpar[4] =-0.545993310748379E+01*fval[0]				    
   	      + 0.165973973290164E+00*fval[1]				    
   	      - 0.198973373528444E+01*fval[2]				    
   	      + 0.408594905200259E+00*fval[3]				    
   	      - 0.112120572260041E+01*fval[4]				    
   	      + 0.799630368687821E+01*fval[5];				    
      fpar[5] = 0.187365927351161E+01*fval[0]				    
   	      + 0.138559587411937E+00*fval[1]				    
   	      + 0.138559587411936E+00*fval[2]				    
   	      - 0.638559587411938E+00*fval[3]				    
   	      + 0.126340726488395E+00*fval[4]				    
   	      - 0.638559587411938E+00*fval[5];	  
   }
   else if (igauss==7)
   {
      fpar[0] = 3.843171599376112E+00*fval[0]				    
    	      - 5.788952456939242E+00*fval[1]				    
    	      + 3.843171599376112E+00*fval[2]				    
    	      + 8.165917142333541E-01*fval[3]				    
    	      - 5.128422945129115E-02*fval[4]				    
    	      + 8.165917142333529E-01*fval[5]				    
    	      - 3.479289940828401E+00*fval[6];
      fpar[1] = 7.737627428203517E+00*fval[0]
    	      - 5.788952456939243E+00*fval[1]				    
    	      - 5.128422945129109E-02*fval[2]				    
    	      + 7.422135885405950E+00*fval[3]				    
    	      - 5.128422945129107E-02*fval[4]				    
    	      - 5.788952456939244E+00*fval[5]
    	      - 3.479289940828401E+00*fval[6];				    
      fpar[2] = 3.843171599376112E+00*fval[0]
    	      + 8.165917142333533E-01*fval[1]				    
    	      - 5.128422945129116E-02*fval[2]				    
    	      + 8.165917142333544E-01*fval[3]				    
    	      + 3.843171599376113E+00*fval[4]
    	      - 5.788952456939242E+00*fval[5]				    
    	      - 3.479289940828402E+00*fval[6];				    
      fpar[3] =-5.674119101307225E+00*fval[0]
    	      + 5.788952456939242E+00*fval[1]				    
    	      - 2.012224097445000E+00*fval[2]				    
    	      - 1.485644212302241E+00*fval[3]
    	      + 5.128422945129113E-02*fval[4]				    
    	      - 1.475392161644655E-01*fval[5]				    
    	      + 3.479289940828402E+00*fval[6];				    
      fpar[4] =-5.674119101307225E+00*fval[0]
    	      - 1.475392161644663E-01*fval[1]				    
    	      + 5.128422945129116E-02*fval[2]
    	      - 1.485644212302242E+00*fval[3]				    
    	      - 2.012224097445000E+00*fval[4]				    
    	      + 5.788952456939242E+00*fval[5]				    
    	      + 3.479289940828402E+00*fval[6];				    
      fpar[5] = 1.974392460120363E+00*fval[0]
    	      - 4.126757274200198E-01*fval[1]				    
    	      + 1.434449581892504E-01*fval[2]				    
    	      + 2.563767706488678E-01*fval[3]				    
    	      + 1.434449581892504E-01*fval[4]				    
    	      - 4.126757274200199E-01*fval[5]				    
    	      - 6.923076923076923E-01*fval[6];	 
   }
   else
      dserror("tri-extrapolation for number of igauss>7 not implemented yet!\n");
   icode++;
} /* endif (icode==0) || (icode==1) */

/*---------calculate function value f(r,s)------*/
if (icode==2)
{
   if (igauss==1)   
      *funval = fpar[0];	   
   else if (igauss==3)   
      *funval = fpar[0]*r + fpar[1]*s + fpar[2];   
   else if (igauss==4) 
      *funval = fpar[0]*r + fpar[1]*s + fpar[2];
   else if (igauss==6) 
      *funval = fpar[0]*r*r + fpar[1]*r*s + fpar[2]*s*s
               +fpar[3]*r   + fpar[4]*s   + fpar[5];
   else if (igauss==7)
      *funval = fpar[0]*r*r + fpar[1]*r*s + fpar[2]*s*s
               +fpar[3]*r   + fpar[4]*s   + fpar[5];     
   else
     dserror ("extropolation for number of gauss point not implemented yet!\n");   
} /* endif (icode==2) */
	
#ifdef DEBUG 
dstrc_exit();
#endif

return ;
} /* end of f2_triex */

/*!---------------------------------------------------------------------
\brief degnenerateed shape functions and their natural derivatives 

<pre>                                                         genk 06/02

In this routine the degnerated shape function and their natural first 
derivative with respect to r is evaluated for T R I A N G L E S 
and  RECTANGLES. They are used for the integrals over the element 
edges.
		     
</pre>
\param  *funct     DOUBLE   (o)    shape functions
\param **deriv     DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2    DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r 	   DOUBLE   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   INT	    (i)    evaluation flag
\return void                                                                       

------------------------------------------------------------------------*/
void f2_degrectri(
                    DOUBLE     *funct, 
                    DOUBLE    **deriv, 
                    DOUBLE      r, 
                    DIS_TYP     typ,
                    INT         option
		 )
{
DOUBLE         rr,rp,rm,r2;

#ifdef DEBUG 
dstrc_enter("f2_degrectri");
#endif

/*----------------------------------------------------------------------*/
/* if option ==0 only funtion evaluation, if option==1 also derivatives */
/*----------------------------------------------------------------------*/
rr = r*r;
rp = ONE+r;
rm = ONE-r;
r2 = ONE-rr;

switch(typ)
{
/*------------------------------------------------ rectangular elements */
case quad4: case tri3:/*-------------- linear rectangular interpolation */
   funct[0] = Q12*rm;
   funct[1] = Q12*rp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -Q12;
      deriv[0][1]=  Q12;
   }
break;
case quad8: case quad9: case tri6:/*---------- quadratic interpolation  */
   funct[0] = -Q12*r*rm;
   funct[1] = r2;
   funct[2] = Q12*r*rp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -Q12+r;
      deriv[0][1]= -TWO*r;
      deriv[0][2]=  Q12+r;                                              
   }
break;
default:
   dserror("unknown typ of interpolation");
break;
} /* end of switch typ */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_degrectri */

/*!---------------------------------------------------------------------
\brief jacobian matrix

<pre>                                                         genk 04/02

In this routine the jacobian matrix and its determinant is calculated
		     
</pre>
\param **xyze      DOUBLE   (i)    nodal coordinates
\param  *funct     DOUBLE   (i)    natural shape functions
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param   iel       INT      (i)    num. of nodes of act. ele
\param  *ele       ELEMENT  (i)    actual element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_jaco(DOUBLE    **xyze,
             DOUBLE     *funct,    
             DOUBLE    **deriv,   
             DOUBLE    **xjm,     
             DOUBLE     *det,          
             INT         iel,        
             ELEMENT    *ele
	    )

{
INT k;

#ifdef DEBUG 
dstrc_enter("f2_jaco");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;
xjm[1][0] = ZERO ;
xjm[1][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     xjm[0][0] += deriv[0][k] * xyze[0][k] ;
     xjm[0][1] += deriv[0][k] * xyze[1][k] ;
     xjm[1][0] += deriv[1][k] * xyze[0][k] ;
     xjm[1][1] += deriv[1][k] * xyze[1][k] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/        
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
    
if(*det<ZERO)
{   
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   printf("NEGATIVE JACOBIAN DETERMINANT\n");
   printf("STOPPING ...");
#ifdef PARALLEL
   dserror("not regulary!\n");
#else
   if (genprob.probtyp==prb_fluid && genprob.numfld==2) 
   {
      fluid_mf(99);
      dserror("regulary!\n");
   }
#ifdef D_FSI
   else if (genprob.probtyp==prb_fsi) 
   {
      dyn_fsi(99);
      dserror("regulary!\n");
   }
#endif
   else dserror("not regulary!\n");
#endif
}
   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_jaco */

/*!---------------------------------------------------------------------
\brief jacobian matrix

<pre>                                                         genk 10/02

In this routine the jacobian matrix and its determinant is calculated;
no check if det<ZERO!!!
		     
</pre>
\param **xyze      DOUBLE   (i)    nodal coordinates
\param  *funct     DOUBLE   (i)    natural shape functions
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param   iel       INT      (i)    num. of nodes of act. ele
\param  *ele       ELEMENT  (i)    actual element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_jaco2(DOUBLE    **xyze,
             DOUBLE     *funct,    
             DOUBLE    **deriv,   
             DOUBLE    **xjm,     
             DOUBLE     *det,          
             INT         iel,        
             ELEMENT    *ele
	    )

{
INT k;

#ifdef DEBUG 
dstrc_enter("f2_jaco");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;
xjm[1][0] = ZERO ;
xjm[1][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     xjm[0][0] += deriv[0][k] * xyze[0][k] ;
     xjm[0][1] += deriv[0][k] * xyze[1][k] ;
     xjm[1][0] += deriv[1][k] * xyze[0][k] ;
     xjm[1][1] += deriv[1][k] * xyze[1][k] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/        
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
      
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_jaco */

/*!---------------------------------------------------------------------
\brief jacobian matrix for integration over the element edges

<pre>                                                         genk 04/02

In this routine the jacobian matrix and its determinant is calculated
		     
</pre>
\param **xyze      DOUBLE   (i)    nodal coordinates
\param  *funct     DOUBLE   (i)    natural shape functions
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param   iel       INT      (i)    num. of nodes of act. ele
\param  *iedgnod   INT      (i)    edgenodes
\return void                                                                       

------------------------------------------------------------------------*/
void f2_edgejaco(
                 DOUBLE    **xyze,
                 DOUBLE     *funct,    
                 DOUBLE    **deriv,   
                 DOUBLE    **xjm,     
                 DOUBLE     *det,          
                 INT         iel,
                 INT        *iedgnod
	        )

{
INT k;
INT node;

#ifdef DEBUG 
dstrc_enter("f2_edgejaco");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     node=iedgnod[k];
     xjm[0][0] += deriv[0][k] * xyze[0][node] ;
     xjm[0][1] += deriv[0][k] * xyze[1][node] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/        
*det = sqrt(xjm[0][0]* xjm[0][0] + xjm[0][1]* xjm[0][1]);
    
if(*det<ZERO)
dserror("NEGATIVE JACOBIAN DETERMINANT ALONG EDGE\n");
   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_edgejaco */

/*!---------------------------------------------------------------------
\brief global derivates

<pre>                                                         genk 04/02

In this routine the global derivatives w.r.t. x,y at point r,s are
calculated.
		     
</pre>
\param **derxy     DOUBLE   (o)    global derivatives wrt. x/y
\param **deriv     DOUBLE   (i)    derivatives of shape functions
\param **xjm       DOUBLE   (i)    jacobian matrix
\param   det       DOUBLE   (i)    jacobian determinant
\param   iel       INT      (i)    number of nodes in actual element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_gder(
              DOUBLE   **derxy,     
	      DOUBLE   **deriv,    
	      DOUBLE   **xjm,      
	      DOUBLE     det,      
	      INT        iel       
	    )
{
INT    k;

#ifdef DEBUG 
dstrc_enter("f2_gder");
#endif

/*------------------------------------------------------- initialistion */
for(k=0;k<iel;k++) /* loop all nodes of the element */
{
   derxy[0][k]=ZERO;
   derxy[1][k]=ZERO;
} /* end loop over iel */

/* NOTE: calculation of inversion is skipped 
         --> change of subscripts and signs in loop!! -----------------*/

for (k=0;k<iel;k++) /* loop all nodes of the element */
{
   derxy[0][k] +=(  xjm[1][1]*deriv[0][k]  + (-xjm[0][1]*deriv[1][k]))/det;
   derxy[1][k] +=((-xjm[1][0]*deriv[0][k]) +   xjm[0][0]*deriv[1][k] )/det;
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_gder */

/*!---------------------------------------------------------------------
\brief global coordinates

<pre>                                                         genk 04/02

In this routine the global coordinates for given shape function values
are set.
		     
</pre>
\param **xyze       DOUBLE   (i)    nodal coordinates
\param  *funct      DOUBLE   (i)    shape functions
\param   iel        DOUBLE   (i)    number of nodes in act. element
\param  *gcoor      DOUBLE   (o)    global coordinates
\return void                                                                       

------------------------------------------------------------------------*/
void f2_gcoor(
              DOUBLE    **xyze,
	      DOUBLE     *funct,            
	      INT         iel,      
	      DOUBLE     *gcoor     
	     )
{
INT i;

#ifdef DEBUG 
dstrc_enter("f2_gcoor");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;

for(i=0;i<iel;i++) /* loop all nodes of the element */
{
   gcoor[0] += funct[i] * xyze[0][i];
   gcoor[1] += funct[i] * xyze[1][i];
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_gcoor */

/*!--------------------------------------------------------------------- 
\brief second global derivatives

<pre>                                                         genk 04/02

In this routine the second global derivatives w.r.t x/y at point r,s
are calculated.
		     
</pre>
\param  *ele 	   ELEMENT  (i)    actual element
\param **xyze      DOUBLE   (i)    element coordinates
\param **xjm 	   DOUBLE   (i)    jacobian matrix
\param **bi 	   DOUBLE   (-)    working array
\param **xder2     DOUBLE   (-)    working array
\param **derxy     DOUBLE   (i)    glob. coord.. deriv.
\param **derxy2    DOUBLE   (o)    2nd. glob. coord. deriv.
\param **deriv2    DOUBLE   (i)    2nd. nat. deriv. of shape funcs
\param   iel	   INT	    (i)    number of nodes of actual ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_gder2(
               ELEMENT     *ele,
	       DOUBLE     **xyze,     
	       DOUBLE     **xjm,      
               DOUBLE     **bi,     
	       DOUBLE     **xder2,  
	       DOUBLE     **derxy,  
	       DOUBLE     **derxy2, 
               DOUBLE     **deriv2, 
	       INT          iel     
	     )
{
INT i;
DOUBLE x00,x01,x02,x10,x11,x12,x20,x21,x22;
DOUBLE det,dum;
DOUBLE r0,r1,r2;

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
det =   x00*x11*x22 + x01*x12*x20 + x10*x21*x02 \
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
} /* end loop over i */
for (i=0;i<iel;i++)
{
   derxy2[0][i]=ZERO;
   derxy2[1][i]=ZERO;
   derxy2[2][i]=ZERO;
} /* end loop over iel */

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   xder2[0][0] += deriv2[0][i] * xyze[0][i];
   xder2[1][0] += deriv2[1][i] * xyze[0][i];
   xder2[2][0] += deriv2[2][i] * xyze[0][i];
   xder2[0][1] += deriv2[0][i] * xyze[1][i];
   xder2[1][1] += deriv2[1][i] * xyze[1][i];
   xder2[2][1] += deriv2[2][i] * xyze[1][i];    
} /* end loop over iel */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i];
   
   derxy2[0][i] += bi[0][0]*r0 + bi[0][1]*r1 + bi[0][2]*r2;
   derxy2[1][i] += bi[1][0]*r0 + bi[1][1]*r1 + bi[1][2]*r2;
   derxy2[2][i] += bi[2][0]*r0 + bi[2][1]*r1 + bi[2][2]*r2;
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_gder2 */



#endif
/*! @} (documentation module close)*/
