/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief routine to calculate element size and stabilisation parameter

<pre>                                                         genk 05/02

   ele->e.f3->iadvec: adevction stab.					 
      0 = no								
      1 = yes								
   ele->e.f3->ipres: pressure stab.					
      0 = no								
      1 = yes								
   ele->e.f3->ivisc: diffusion stab.					
      0 = no								
      1 = GLS-  							
      2 = GLS+  							
   ele->e.f3->icont: continuity stab.					
      0 = no								
      1 = yes								
   ele->e.f3->istapa: version of stab. parameter			
      35 = diss wall instationary					
      36 = diss wall stationanary					
   ele->e.f3->norm_P: p-norm						
      p = 1<=p<=oo							
      0 = max.-norm (p=oo)						
   ele->e.f3->mk: higher order elements control flag			
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)		
      1 = min(1/3,2*C)  						
     -1 = mk=1/3  (--> element order via approx. nodal-distance)	
   ele->e.f3->ihele[]:  						
      x/y/z = length-def. for velocity/pressure/continuity stab 	
      0 = don't compute 						
      1 = sqrt(area)							
      2 = area equivalent diameter					
      3 = diameter/sqrt(2)						
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle) 	
      5 = streamlength (element length in flow direction		
   ele->e.f3->ninths: number of integration points for streamlength	
      1 = at center of element  					
      2 = at every int pt used for element.-stab.-matrices		
   ele->e.f3->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f3->clamb \							
   ele->e.f3->c1     |_>> stabilisation constants (input)		
   ele->e.f3->c2     |  						
   ele->e.f3->c3    /							
   ele->e.f3->istrle: has streamlength to be computed			
   ele->e.f3->iarea: calculation of area length 			
   ele->e.f3->iduring: calculation during int.-pt.loop  		
   ele->e.f3->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f3->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f3->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f3->hk[i]: "element sizes" (vel / pre / cont) 		  
   ele->e.f3->idiaxy: has diagonals to be computed			
   dynvar->tau[0]: stability parameter momentum / velocity (tau_mu)	
   dynvar->tau[1]: stability parameter momentum / pressure (tau_mp)	
   dynvar->tau[2]: stability parameter continuity (tau_c)
</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *data    FLUID_DATA	       (i)
\param  *dynvar  FLUID_DYN_CALC        (i/o)
\param  *funct   double 	       (-)   shape functions
\param **deriv   double 	       (-)   deriv. of shape funcs
\param **deriv2  double 	       (-)   2nd deriv. of sh. funcs
\param **derxy   double 	       (-)   global derivatives
\param **xjm     double 	       (-)   jacobian matrix
\param **evel    double 	       (i)   element velocities
\param  *velint  double 	       (-)   vel. at integr. point
\param **cutp    double 	       (-)   cutting points
\return void             

------------------------------------------------------------------------*/			     
void f3_calelesize(
	           ELEMENT         *ele,
		   FLUID_DATA      *data,
		   FLUID_DYN_CALC  *dynvar,
	           double          *funct,
	           double         **deriv,
	           double         **deriv2,	       
                   double         **derxy,
		   double         **xjm,
		   double         **evel,	       
		   double          *velint,
		   double         **wa1		
                  )
{
int ieval = 0;       /* evaluation flag			                */
int i,ilen, inod;    /* simply a counter	        		*/
int istrnint;        /* evaluation flag		     	                */
int ishvol;          /* evaluation flag		        	        */
int ntyp;            /* element type (TET or HEX)  		        */
int actmat;          /* number of actual material		        */
int iel;             /* number of nodes of actual element               */
double visc;         /* fluid viscosity                                 */
double det;          /* determinant of jacobian                         */
double vol;          /* element volume                                  */
double val;          /* temporary calculation value                     */
double velno;        /* velocity norm                                   */
double strle;        /* streamlength                                    */
double e1,e2,e3;     /* natural coordinates of inegration point         */
double fac,facr;     /* factors                                         */
double facs,fact;    /* factors                                         */
double velino[3];    /* normed velocity vector at integration point     */
DIS_TYP typ;

#ifdef DEBUG 
dstrc_enter("f3_calelesize");
#endif		

/*---------------------------------------------------------- initialise */
ntyp   = ele->e.f3->ntyp;
iel    = ele->numnp;
typ    = ele->distyp;

istrnint = ele->e.f3->istrle * ele->e.f3->ninths;
ishvol  = dynvar->ishape * ele->e.f3->ivol;

/*----------------------------------------------------------------------*
 | calculations at element center: area & streamlength                  |
 | NOTE:                                                                |
 |    volume is always calculated using only 1 integrationpoint         |  
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,e3, facr,facs,fact with their constant values in the   |
 |         calls of f3_hex / f3_tet!!!!!!                               |
 *----------------------------------------------------------------------*/

if (ishvol==1)
{
   vol   = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(ntyp)
   {
   case 1:   /* --> hex - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      e3   = data->qxg[0][0];
      fact = data->qwgt[0][0]; 
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,2);
   break;        
   case 2:  /* --> tet - element */
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      e3   = data->txgs[0][0]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,2);      
   break;
   default:
      dserror("ntyp unknown!"); 
   } /*end switch(ntyp) */
   ieval++;
/* ------------------------------------------- compute jacobian matrix */        
   f3_jaco(funct,deriv,xjm,&det,ele,iel);
   fac=facr*facs*fact*det;
   vol += fac;
   if (istrnint==1)    /* compute streamlength */
   {
      f3_veli(velint,funct,evel,iel);      
      f3_gder(derxy,deriv,xjm,wa1,det,iel); 
      ieval++;
      val = ZERO;
      velno=sqrt( velint[0]*velint[0] \
                + velint[1]*velint[1] \
		+ velint[2]*velint[2]);
      if(velno>=EPS6)
      {
	 velino[0] = velint[0]/velno;
         velino[1] = velint[1]/velno;
         velino[2] = velint[2]/velno;
      }
      else
      {
         velino[0] = ONE;
	 velino[1] = ZERO;
	 velino[2] = ZERO;
      }         
      for (inod=0;inod<iel;inod++) /* loop element nodes */
      {
         val += FABS(velino[0]*derxy[0][inod] \
	            +velino[1]*derxy[1][inod] \
		    +velino[2]*derxy[2][inod]);
      } /* end of loop over elements */
      strle=TWO/val;      
   } /* endif (istrnint==1) */
/*--------------------------------------------------- set element sizes *
  ----loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for(ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f3->ihele[ilen]==1)
         ele->e.f3->hk[ilen] = pow(vol,(ONE/THREE));
      else if (ele->e.f3->ihele[ilen]==2)
         ele->e.f3->hk[ilen] = pow((SIX*vol/PI),(ONE/THREE));
      else if (ele->e.f3->ihele[ilen]==3)
         ele->e.f3->hk[ilen] = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);
      else if (ele->e.f3->ihele[ilen]==4) 
         dserror("ihele[i] = 4: calculation of element size not possible!!!");
         else if (ele->e.f3->ninths==1)   
         ele->e.f3->hk[ilen] = strle; 
   } /* end of loop over ilen */
} /* endif (ishvol==1) */
/*----------------------------------------------------------------------*
 | calculations at element center: only streamlength                    |
 | NOTE:                                                                |
 |    volume is always calculated using only 1 integrationpoint         |  
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,e3, facr,facs,fact with their constant values in the   |
 |         calls of f3_hex / f3_tet!!!!!!                               |
 *----------------------------------------------------------------------*/
else if (istrnint==1 && ishvol !=1) 
{
   vol   = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(ntyp)
   {
   case 1:   /* --> hex - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      e3   = data->qxg[0][0];
      fact = data->qwgt[0][0]; 
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,2);
   break;        
   case 2:  /* --> tet - element */
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      e3   = data->txgs[0][0]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,2);
   break;
   default:
      dserror("ntyp unknown!"); 
   } /*end switch(ntyp) */ 
/* ------------------------------------------- compute jacobian matrix */        
   f3_jaco(funct,deriv,xjm,&det,ele,iel);
/* --------------------------------------------- compute stream length */    
   f3_veli(velint,funct,evel,iel);	
   f3_gder(derxy,deriv,xjm,wa1,det,iel); 
   ieval++;
   val = ZERO;
   velno=sqrt( velint[0]*velint[0] \
   	     + velint[1]*velint[1] \
     	     + velint[2]*velint[2]);
   if(velno>=EPS6)
   {
      velino[0] = velint[0]/velno;
      velino[1] = velint[1]/velno;
      velino[2] = velint[2]/velno;
   }
   else
   {
      velino[0] = ONE;
      velino[1] = ZERO;
      velino[2] = ZERO;
   }	     
   for (inod=0;inod<iel;inod++) /* loop element nodes */
   { 
      val += FABS(velino[0]*derxy[0][inod] \
     		 +velino[1]*derxy[1][inod] \
     		 +velino[2]*derxy[2][inod]);      
   } /* end of loop over element nodes */
   strle=TWO/val;	  
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f3->ihele[ilen]==5)
         ele->e.f3->hk[ilen] = strle;   
   } /* end of loop over ilen */
} /* endif (istrnint==1 && ishvol !=1) */

/*----------------------------------------------------------------------*
  calculate stabilisation parameter
 *----------------------------------------------------------------------*/
if(ele->e.f3->istapc==1 || istrnint==1)
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
         e3   = data->qxg[0][0];
         fact = data->qwgt[0][0]; 
         f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,2);
      break;
      case 2:       /* --> tet - element */              
         e1   = data->txgr[0][0];
         facr = data->twgt[0][0];
         e2   = data->txgs[0][0];
         facs = ONE;
         e3   = data->txgs[0][0]; 
         fact = ONE;
         f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,2);
      break;      
      default:
         dserror("ntyp unknown!");
      } /* end switch (ntyp) */
      f3_veli(velint,funct,evel,iel);
   break;
   case 1:            
      f3_veli(velint,funct,evel,iel);
   break;
   case 2:
   break;
   default:
      dserror("wrong value for ieval");
   } /* end switch (ieval) */
/*----------------------------------- calculate stabilisation parameter */               
   actmat=ele->mat-1;
   visc = mat[actmat].m.fluid->viscosity;
   f3_calstabpar(ele,dynvar,velint,visc,iel,ntyp,-1);    
} /* endif (ele->e.f3->istapc==1 || istrnint==1) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calelesize */

/*!---------------------------------------------------------------------                                         
\brief routine to calculate element size and stabilisation parameter

<pre>                                                         genk 05/02

in this routine the element size and the stabilisation parameter 
is calculated for one element during the integration loop
		     
</pre>
\param  *ele     ELEMENT	        (i)    actual element
\param  *dynvar  FLUID_DYN_CALC         (i/o)
\param  *velint  double 		(-)    vel at intpoint
\param  *derxy   double 		(-)    global derivatives
\param   visc    double 		(i)    fluid viscosity
\param   iel     int		        (i)    act. num. of ele nodes
\param   ntyp    int		        (i)    element type
\return void                                               
\sa f3_calelesize()                               

------------------------------------------------------------------------*/
void f3_calelesize2(
	             ELEMENT         *ele,
		     FLUID_DYN_CALC  *dynvar,
                     double          *velint,	       
                     double         **derxy,	       	
		     double           visc,
		     int              iel,
		     int              ntyp
		  ) 
{
int    ilen, inod; /* simply a counter                                  */
int    istrnint;   /* evaluation flag                                   */
double strle;      /* stream length                                     */
double val;	   /* temporary calculation value                       */
double velno;	   /* velocity norm                                     */
double velino[3];  /* normed velocity vector at integration point       */

#ifdef DEBUG 
dstrc_enter("f3_calelesize2");
#endif

/*---------------------------------------------------------- initialise */
istrnint = ele->e.f3->istrle * ele->e.f3->ninths;
val = ZERO;

if (istrnint==2)
{
/*------------------------------------------------ compute streamlength */
   velno=sqrt( velint[0]*velint[0] \
   	     + velint[1]*velint[1] \
     	     + velint[2]*velint[2]);
   if(velno>=EPS6)
   {
      velino[0] = velint[0]/velno;
      velino[1] = velint[1]/velno;
      velino[2] = velint[2]/velno;
   }
   else
   {
      velino[0] = ONE;
      velino[1] = ZERO;
      velino[2] = ZERO;
   }	     
   for (inod=0;inod<iel;inod++)
   {
      val += FABS(velino[0]*derxy[0][inod] \
     		 +velino[1]*derxy[1][inod] \
     		 +velino[2]*derxy[2][inod]);      
   }
   strle=TWO/val;
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (ele->e.f3->ihele[ilen]==5)
         ele->e.f3->hk[ilen] = strle;   
   } /* end of loop over ilen */
} /* endif (istrnint==2) */
/*----------------------------------- calculate stabilisation parameter */               
f3_calstabpar(ele,dynvar,velint,visc,iel,ntyp,1); 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_calelesize2 */
#endif
