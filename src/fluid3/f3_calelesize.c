/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

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

   ele->e.f3->stabi.gls->iadvec: adevction stab.					 
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->ipres: pressure stab.					
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->ivisc: diffusion stab.					
      0 = no								
      1 = GLS-  							
      2 = GLS+  							
   ele->e.f3->stabi.gls->icont: continuity stab.					
      0 = no								
      1 = yes								
   ele->e.f3->stabi.gls->istapa: version of stab. parameter			
      35 = diss wall instationary					
      36 = diss wall stationanary					
   ele->e.f3->stabi.gls->norm_P: p-norm						
      p = 1<=p<=oo							
      0 = max.-norm (p=oo)						
   ele->e.f3->stabi.gls->mk: higher order elements control flag			
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)		
      1 = min(1/3,2*C)  						
     -1 = mk=1/3  (--> element order via approx. nodal-distance)	
   ele->e.f3->stabi.gls->ihele[]:  						
      x/y/z = length-def. for velocity/pressure/continuity stab 	
      0 = don't compute 						
      1 = sqrt(area)							
      2 = area equivalent diameter					
      3 = diameter/sqrt(2)						
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle) 	
      5 = streamlength (element length in flow direction		
   ele->e.f3->stabi.gls->ninths: number of integration points for streamlength	
      1 = at center of element  					
      2 = at every INT pt used for element.-stab.-matrices		
   ele->e.f3->stabi.gls->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f3->stabi.gls->clamb \							
   ele->e.f3->c1               |_>> stabilisation constants (input)		
   ele->e.f3->c2               |  						
   ele->e.f3->c3              /							
   ele->e.f3->stabi.gls->istrle: has streamlength to be computed			
   ele->e.f3->stabi.gls->iarea: calculation of area length 			
   ele->e.f3->stabi.gls->iduring: calculation during INT.-pt.loop  		
   ele->e.f3->stabi.gls->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f3->stabi.gls->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f3->stabi.gls->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f3->hk[i]: "element sizes" (vel / pre / cont) 		  
   ele->e.f3->stabi.gls->idiaxy: has diagonals to be computed			
   dynvar->tau[0]: stability parameter momentum / velocity (tau_mu)	
   dynvar->tau[1]: stability parameter momentum / pressure (tau_mp)	
   dynvar->tau[2]: stability parameter continuity (tau_c)
</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *data    FLUID_DATA	       (i)
\param  *dynvar  FLUID_DYN_CALC        (i/o)
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of sh. funcs
\param **derxy   DOUBLE 	       (-)   global derivatives
\param **xjm     DOUBLE 	       (-)   jacobian matrix
\param **evel    DOUBLE 	       (i)   element velocities
\param  *velint  DOUBLE 	       (-)   vel. at integr. point
\param **cutp    DOUBLE 	       (-)   cutting points
\return void             

------------------------------------------------------------------------*/			     
void f3_calelesize(
	           ELEMENT         *ele,
		   FLUID_DATA      *data,
		   FLUID_DYN_CALC  *dynvar,
	           DOUBLE          *funct,
	           DOUBLE         **deriv,
	           DOUBLE         **deriv2,	       
                   DOUBLE         **derxy,
		   DOUBLE         **xjm,
		   DOUBLE         **evel,	       
		   DOUBLE          *velint,
		   DOUBLE         **wa1		
                  )
{
INT ieval = 0;       /* evaluation flag			                */
INT ilen, inod;    /* simply a counter	        		*/
INT istrnint;        /* evaluation flag		     	                */
INT ishvol;          /* evaluation flag		        	        */
INT ntyp;            /* element type (TET or HEX)  		        */
INT actmat;          /* number of actual material		        */
INT iel;             /* number of nodes of actual element               */
DOUBLE visc;         /* fluid viscosity                                 */
DOUBLE det;          /* determinant of jacobian                         */
DOUBLE vol;          /* element volume                                  */
DOUBLE val;          /* temporary calculation value                     */
DOUBLE velno;        /* velocity norm                                   */
DOUBLE strle;        /* streamlength                                    */
DOUBLE e1,e2,e3;     /* natural coordinates of inegration point         */
DOUBLE fac,facr;     /* factors                                         */
DOUBLE facs,fact;    /* factors                                         */
DOUBLE velino[3];    /* normed velocity vector at integration point     */
DIS_TYP typ;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calelesize");
#endif		

/*---------------------------------------------------------- initialise */
ntyp   = ele->e.f3->ntyp;
iel    = ele->numnp;
typ    = ele->distyp;
gls    = ele->e.f3->stabi.gls;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

istrnint = gls->istrle * gls->ninths;
ishvol  = dynvar->ishape * gls->iareavol;

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
      if (gls->ihele[ilen]==1)
         ele->e.f3->hk[ilen] = pow(vol,(ONE/THREE));
      else if (gls->ihele[ilen]==2)
         ele->e.f3->hk[ilen] = pow((SIX*vol/PI),(ONE/THREE));
      else if (gls->ihele[ilen]==3)
         ele->e.f3->hk[ilen] = pow((SIX*vol/PI),(ONE/THREE))/sqrt(THREE);
      else if (gls->ihele[ilen]==4) 
         dserror("ihele[i] = 4: calculation of element size not possible!!!");
         else if (gls->ninths==1)   
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
      if (gls->ihele[ilen]==5)
         ele->e.f3->hk[ilen] = strle;   
   } /* end of loop over ilen */
} /* endif (istrnint==1 && ishvol !=1) */

/*----------------------------------------------------------------------*
  calculate stabilisation parameter
 *----------------------------------------------------------------------*/
if(gls->istapc==1 || istrnint==1)
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
\param  *velint  DOUBLE 		(-)    vel at intpoint
\param  *derxy   DOUBLE 		(-)    global derivatives
\param   visc    DOUBLE 		(i)    fluid viscosity
\param   iel     INT		        (i)    act. num. of ele nodes
\param   ntyp    INT		        (i)    element type
\return void                                               
\sa f3_calelesize()                               

------------------------------------------------------------------------*/
void f3_calelesize2(
	             ELEMENT         *ele,
		     FLUID_DYN_CALC  *dynvar,
                     DOUBLE          *velint,	       
                     DOUBLE         **derxy,	       	
		     DOUBLE           visc,
		     INT              iel,
		     INT              ntyp
		  ) 
{
INT    ilen, inod; /* simply a counter                                  */
INT    istrnint;   /* evaluation flag                                   */
DOUBLE strle;      /* stream length                                     */
DOUBLE val;	   /* temporary calculation value                       */
DOUBLE velno;	   /* velocity norm                                     */
DOUBLE velino[3];  /* normed velocity vector at integration point       */
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calelesize2");
#endif		

/*---------------------------------------------------------- initialise */
gls      = ele->e.f3->stabi.gls;
istrnint = gls->istrle * gls->ninths;
val = ZERO;

if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");

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
      if (gls->ihele[ilen]==5)
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
