/*!----------------------------------------------------------------------
\file
\brief shape functions and their natural derivatives for fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
static double Q12 = ONE/TWO;
static double Q14 = ONE/FOUR;
static double Q16 = ONE/SIX;
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
\param  *funct     double   (o)    shape functions
\param **deriv     double   (o)    1st natural deriv. of shape funct.
\param **deriv2    double   (o)    2nd natural deriv. of shape funct.
\param   r 	   double   (i)    coordinate
\param   s 	   double   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   int	    (i)    evaluation flag
\return void                                                                       

------------------------------------------------------------------------*/
void f2_rec(
            double     *funct,     
            double    **deriv,    
            double    **deriv2,   
	    double      r,        
            double      s,        
            DIS_TYP     typ,      
            int         icode     
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
\brief shape functions and their natural derivatives for triangles

<pre>                                                         genk 06/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s are evaluated for 
T R I A N G L E S 
		     
</pre>
\param  *funct     double   (o)    shape functions
\param **deriv     double   (o)    1st natural deriv. of shape funct.
\param **deriv2    double   (o)    2nd natural deriv. of shape funct.
\param   r 	   double   (i)    coordinate
\param   s 	   double   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   int	    (i)    evaluation flag
\return void                                                                       

------------------------------------------------------------------------*/
void f2_tri(
            double     *funct,      
            double    **deriv,    
            double    **deriv2,   
	    double      r,        
            double      s,	  
            DIS_TYP     typ,	  
            int         icode	  
	   )
{

double rr,rs,rrr,ss,sss,rrs,rss;

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
   rrr=rr*r;
   rrs=rr*s;
   rss=r*ss;
   sss=ss*s;

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
} /* end of f2_rec */

/*!---------------------------------------------------------------------                                         
\brief jacobian matrix

<pre>                                                         genk 04/02

In this routine the jacobian matrix and its determinant is calculated
		     
</pre>
\param  *funct     double   (i)    natural shape functions
\param **deriv     double   (i)    natural deriv. of shape funcs
\param **xjm       double   (o)    jacobian matrix
\param  *det       double   (o)    determinant of jacobian matrix
\param  *ele       ELEMENT  (i)    actual element
\param   iel       int      (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_jaco(double     *funct,    
             double    **deriv,   
             double    **xjm,     
             double     *det,     
             ELEMENT    *ele,     
             int         iel        
	    )

{
int k;

#ifdef DEBUG 
dstrc_enter("f2_jaco");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = 0.0 ;
xjm[0][1] = 0.0 ;
xjm[1][0] = 0.0 ;
xjm[1][1] = 0.0 ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     xjm[0][0] += deriv[0][k] * ele->node[k]->x[0] ;
     xjm[0][1] += deriv[0][k] * ele->node[k]->x[1] ;
     xjm[1][0] += deriv[1][k] * ele->node[k]->x[0] ;
     xjm[1][1] += deriv[1][k] * ele->node[k]->x[1] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/        
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
    
if(*det<0.0)
{   
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}
   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_jaco */

/*!---------------------------------------------------------------------                                         
\brief global derivates

<pre>                                                         genk 04/02

In this routine the global derivatives w.r.t. x,y at point r,s are
calculated.
		     
</pre>
\param **derxy     double   (o)    global derivatives wrt. x/y
\param **deriv     double   (i)    derivatives of shape functions
\param **xjm       double   (i)    jacobian matrix
\param   det       double   (i)    jacobian determinant
\param   iel       int      (i)    number of nodes in actual element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_gder(
              double   **derxy,     
	      double   **deriv,    
	      double   **xjm,      
	      double     det,      
	      int        iel       
	    )
{
int    k;

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
\param *funct      double   (i)    shape functions
\param *ele        double   (i)    actual element
\param  iel        double   (i)    number of nodes in act. element
\param *gcoor      double   (o)    global coordinates
\return void                                                                       

------------------------------------------------------------------------*/
void f2_gcoor(
              double     *funct,      
              ELEMENT    *ele,      
	      int         iel,      
	      double     *gcoor     
	     )
{
int i;

#ifdef DEBUG 
dstrc_enter("f2_gcoor");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;

for(i=0;i<iel;i++) /* loop all nodes of the element */
{
   gcoor[0] += funct[i] * ele->node[i]->x[0];
   gcoor[1] += funct[i] * ele->node[i]->x[1];
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
\param **xjm 	   double   (i)    jacobian matrix
\param **bi 	   double   (-)    working array
\param **xder2     double   (-)    working array
\param **derxy     double   (i)    glob. coord.. deriv.
\param **derxy2    double   (o)    2nd. glob. coord. deriv.
\param **deriv2    double   (i)    2nd. nat. deriv. of shape funcs
\param   iel	   int	    (i)    number of nodes of actual ele
\return void                                                                       

------------------------------------------------------------------------*/
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
   xder2[0][0] += deriv2[0][i] * ele->node[i]->x[0];
   xder2[1][0] += deriv2[1][i] * ele->node[i]->x[0];
   xder2[2][0] += deriv2[2][i] * ele->node[i]->x[0];
   xder2[0][1] += deriv2[0][i] * ele->node[i]->x[1];
   xder2[1][1] += deriv2[1][i] * ele->node[i]->x[1];
   xder2[2][1] += deriv2[2][i] * ele->node[i]->x[1];
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
