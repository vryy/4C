/*!---------------------------------------------------------------------
\file
\brief fluid3 prototypes

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/************************************************************************
 | f3_calele.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief control routine for element integration of fluid3

<pre>                                                         genk 05/02

This routine controls the element evaluation:
-actual vel. and pres. variables are set
-stabilisation parameters are calculated
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
-stiffness matrix and load vectors are permuted for assembling
-element load vector due to dirichlet conditions is calculated				      
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\param   init	         INT	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_calele(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,
                ARRAY          *estif_global,
                ARRAY          *emass_global, 
	        ARRAY          *etforce_global,
	        ARRAY          *eiforce_global,
		ARRAY          *edforce_global,	
		INT            *hasdirich,
		INT            *hasext,		
		INT             init
	       );

/************************************************************************
 | f3_calelesize.c                                                      |
 ************************************************************************/
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
      2 = at every INT pt used for element.-stab.-matrices		
   ele->e.f3->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f3->clamb \							
   ele->e.f3->c1     |_>> stabilisation constants (input)		
   ele->e.f3->c2     |  						
   ele->e.f3->c3    /							
   ele->e.f3->istrle: has streamlength to be computed			
   ele->e.f3->iarea: calculation of area length 			
   ele->e.f3->iduring: calculation during INT.-pt.loop  		
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
                  );
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
		  );

/************************************************************************
 | f3_calextrhs.c                                                       |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief galerkin part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

                   /
   + (1-THETA)*dt |  v * b_old   d_omega
                 /  

               /
   + THETA*dt |  v * b   d_omega
             /      	  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param   *edeadn      DOUBLE          (i)    ele dead load at n
\param   *edeadng     DOUBLE          (i)    ele dead load at n+1
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgalexfv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,     
		  DOUBLE          *funct,       
                  DOUBLE          *edeadn,
		  DOUBLE          *edeadng,
		  DOUBLE           fac,      
		  INT              iel       
              );
/*!--------------------------------------------------------------------- 
\brief stabilisation part of external forces for vel dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:

                     /
   + thetas(l,r)*dt |  taum_mu * u * grad(v) * b^   d_omega
                   /  

                       /
   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b^  d_omega
                     /      	  		      

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn    

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *ele         ELEMENT	      (i)    actual element
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives
\param  **derxy2,     DOUBLE	      (i)    2nd. global derivatives    
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  ihoel       INT	      (i)    flag for higer ord. ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabexfv(
                    FLUID_DYN_CALC  *dynvar, 
                    ELEMENT         *ele,  
                    DOUBLE          *eforce,     
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,      
                    DOUBLE          *edead,
	 	    DOUBLE          *velint,  
		    DOUBLE           fac,      
                    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel,
		    INT              flag      
                   );
/*!--------------------------------------------------------------------- 
\brief stabilisation part of external forces for pre dofs

<pre>                                                         genk 09/02

In this routine the stabilisation part of the time forces for pre dofs
is calculated:

                      /
   - thetas(l,r)*dt  |  tau_mp * grad(q) * b^  d_omega
                    /	  		      

This routine is called twice with different values:
1. values are added to Iteration RHS (evaluation at n+1):
    thetas(l,r) = THETA*dt
    b^ = b = deadng
2. values are added to Time RHS (evaluation at n):
   thetas(l,r) = (1-THETA)*dt
   b^ = b_old = deadn    

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry 
    eforce[2*iel]     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param  **derxy,      DOUBLE	      (i)    global derivatives    
\param   *edead       DOUBLE          (i)    ele dead load at n or n+1
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\param	  flag	      INT	      (i)    flag for n or n+1
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabexfp(
                    FLUID_DYN_CALC  *dynvar, 
                    DOUBLE          *eforce,     
		    DOUBLE         **derxy,       
                    DOUBLE          *edead,  
		    DOUBLE           fac,      
		    INT              iel,
		    INT              flag      
                   ) ;		   	      

/************************************************************************
 | f3_calfuncderiv.c                                                    |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief shape functions and their natural derivatives for hexaeder

<pre>                                                         genk 05/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for 
 H E X A H E D E R
 
   Numbering of the nodes:

                           ^ t
                           |
                           |
                           |
                    8      |  15        7
                    o---------o---------o
                   /|                  /|
                  / |                 / |
                 /  |                /	|
              16o   |     o       14o   |
               /    o20       o    /	o19
              /     |             /     |
             /      |  13      6 /	|
          5 o---------o---------o	|
            |   o   |     o   	|   o   |  ---------->
            |       o---------o-|-------o           s
            |      / 4       11 |      /3
            |     /             |     /
          17o    /    o         o18  /
            | 12o         o     |   o10
            |  /                |  /
            | /                 | /
            |/	                |/
            o---------o---------o
	    1	/     9         2
	       /
	      /
	     /
	    r  

   PROBLEM: GID has a different numbering of the element nodes than this one.
            So either the shape functions for hex20 and hex27 (see drawing)
	    has to be adapted or during the input phase the numbering has to 
	    be adapted to the shape functions. 
	    This is all in progress and should be done for fluid3 and
	    brick1 the same way!!!!	   

</pre>
\param  *funct     DOUBLE   (o)    shape functions
\param **deriv     DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2    DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r 	   DOUBLE   (i)    coordinate
\param   s 	   DOUBLE   (i)    coordinate
\param   t 	   DOUBLE   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   INT	    (i)    evaluation flag
\return void                                                                       
\warning shape functions for hex20/hex27/tet10 not implemented yet!!!

------------------------------------------------------------------------*/
 void f3_hex(
            DOUBLE     *funct,     
            DOUBLE    **deriv,    
            DOUBLE    **deriv2,   
	    DOUBLE      r,        
            DOUBLE      s,        
            DOUBLE      t,        
            DIS_TYP     typ,      
            INT         icode     
	   );
/*!--------------------------------------------------------------------- 
\brief shape functions and their natural derivatives for tetraeder

<pre>                                                         genk 08/02

In this routine the shape functions and their natural first and second
derivatives with respect to r/s/t are evaluated for 
T E T R A E D E R  
		     
</pre>
\param  *funct     DOUBLE   (o)    shape functions
\param **deriv     DOUBLE   (o)    1st natural deriv. of shape funct.
\param **deriv2    DOUBLE   (o)    2nd natural deriv. of shape funct.
\param   r 	   DOUBLE   (i)    coordinate
\param   s 	   DOUBLE   (i)    coordinate
\param   t 	   DOUBLE   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   INT	    (i)    evaluation flag
\return void                                                                       
\warning shape functions for TET10 not implemented yet!!!

------------------------------------------------------------------------*/
void f3_tet(
            DOUBLE     *funct,     
            DOUBLE    **deriv,    
            DOUBLE    **deriv2,   
	    DOUBLE      r,        
            DOUBLE      s,        
            DOUBLE      t,        
            DIS_TYP     typ,      
            INT         icode     
	   );	   
/*!--------------------------------------------------------------------- 
\brief jacobian matrix

<pre>                                                         genk 05/02

In this routine the jacobian matrix and its determinant is calculated
		     
</pre>
\param  *funct     DOUBLE   (i)    natural shape functions
\param **deriv     DOUBLE   (i)    natural deriv. of shape funcs
\param **xjm       DOUBLE   (o)    jacobian matrix
\param  *det       DOUBLE   (o)    determinant of jacobian matrix
\param  *ele       ELEMENT  (i)    actual element
\param   iel       INT      (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_jaco(DOUBLE     *funct,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             ELEMENT    *ele,
             INT         iel);
/*!--------------------------------------------------------------------- 
\brief global derivates

<pre>                                                         genk 05/02

In this routine the global derivatives w.r.t. x,y,z at point r,s,t are
calculated.
		     
</pre>
\param **derxy     DOUBLE   (o)    global derivatives wrt. x/y/z
\param **deriv     DOUBLE   (i)    derivatives of shape functions
\param **xjm       DOUBLE   (i)    jacobian matrix
\param **xji       DOUBLE   (-)    inverse of jacobian
\param   det       DOUBLE   (i)    jacobian determinant
\param   iel       INT      (i)    number of nodes in actual element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_gder(
              DOUBLE   **derxy,     
	      DOUBLE   **deriv,    
	      DOUBLE   **xjm,      
	      DOUBLE   **xji,      
	      DOUBLE     det,      
	      INT        iel       
	    );
/*!--------------------------------------------------------------------- 
\brief global coordinates

<pre>                                                         genk 05/02

In this routine the global coordinates for given shape function values
are set.
		     
</pre>
\param *funct      DOUBLE   (i)    shape functions
\param *ele        DOUBLE   (i)    actual element
\param  iel        DOUBLE   (i)    number of nodes in act. element
\param *gcoor      DOUBLE   (o)    global coordinates
\return void                      

------------------------------------------------------------------------*/
void f3_gcoor(
              DOUBLE     *funct,     
              ELEMENT    *ele,
	      INT         iel,      
	      DOUBLE     *gcoor      
	     );
/*!--------------------------------------------------------------------- 
\brief second global derivatives

<pre>                                                         genk 05/02

In this routine the second global derivatives w.r.t x/y/z at point r,s,t
are calculated.
		     
</pre>
\param  *ele 	   ELEMENT  (i)    actual element
\param **xjm 	   DOUBLE   (i)    jacobian matrix
\param **bm 	   DOUBLE   (-)    working array
\param **xder2     DOUBLE   (-)    working array
\param **derxy     DOUBLE   (i)    glob. coord.. deriv.
\param **derxy2    DOUBLE   (o)    2nd. glob. coord. deriv.
\param **deriv2    DOUBLE   (i)    2nd. nat. deriv. of shape funcs
\param   iel	   INT	    (i)    number of nodes of actual ele
\return void        

------------------------------------------------------------------------*/
void f3_gder2(
               ELEMENT     *ele,
	       DOUBLE     **xjm,            
               DOUBLE     **bm,
	       DOUBLE     **xder2,
	       DOUBLE     **derxy,
	       DOUBLE     **derxy2,
               DOUBLE     **deriv2,
	       INT          iel
	     );
	     	     	    	     
/************************************************************************
 | f3_calgalmat.c                                                       |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief evaluate galerkin part of Kvv

<pre>                                                         genk 05/02

In this routine the galerkin part of matrix Kvv is calculated:

    /
   |  2 * nue * eps(v) : eps(u)   d_omega
  /

    /
   |  v * u_old * grad(u)     d_omega
  /

    /
   |  v * u * grad(u_old)     d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Kvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param  *velint    DOUBLE	   (i)    vel at INT point
\param **vderxy    DOUBLE	   (i)    global vel derivatives
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   iel	   INT  	   (i)	  number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calkvv(
                FLUID_DYN_CALC  *dynvar,
		DOUBLE         **estif,
		DOUBLE          *velint,
		DOUBLE         **vderxy,
		DOUBLE          *funct,
		DOUBLE         **derxy,
		DOUBLE           fac,
		DOUBLE           visc,
		INT              iel		 
              );
/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Kvp

<pre>                                                         genk 05/02

In this routine the galerkin part of matrix Kvp is calculated:

    /
   |  - div(v) * p     d_omega
  /

    /
   | - q * div(u)      d_omega
  / 

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				    
      --> Kvp is stored in estif[(0..(3*iel-1)][(4*iel)..(4*iel-1)] 
      --> Kpv is stored in estif[((3*iel)..(4*iel-1)][0..(3*iel-1)] 
      
</pre>
\param **estif     DOUBLE	   (i/o)  ele stiffness matrix
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param **derxy     DOUBLE	   (i)    global coord. deriv.
\param   fac       DOUBLE	   (i)    weighting factor
\param   iel       INT  	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calkvp(
		DOUBLE         **estif,
		DOUBLE          *funct,
		DOUBLE         **derxy,
		DOUBLE           fac,
		INT              iel                
              );
/*!---------------------------------------------------------------------
\brief evaluate galerkin part of Mvv

<pre>                                                         genk 05/02

In this routine the galerkin part of matrix Mvv is calculated:

    /
   |  v * u    d_omega
  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Mvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)] 
      
</pre>
\param **emass     DOUBLE	   (i/o)  ele mass matrix
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac       DOUBLE	   (i)    weighting factor
\param   iel       INT   	   (i)    number of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calmvv(
		DOUBLE         **estif,
		DOUBLE          *funct,
		DOUBLE           fac,
		INT              iel  
              );
	      	      	      

/************************************************************************
 | f3_calint.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element

<pre>                                                         genk 05/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadn    DOUBLE	   (-)    ele dead load (selfweight) at n 
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *vel2int   DOUBLE	   (-)    vel at integration point
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_calint(
               FLUID_DATA      *data, 
	       ELEMENT         *ele,
	       FLUID_DYN_CALC  *dynvar,
               INT             *hasext,
               DOUBLE         **estif,
	       DOUBLE         **emass,
	       DOUBLE          *etforce,
	       DOUBLE          *eiforce,
	       DOUBLE          *funct,
	       DOUBLE         **deriv,
	       DOUBLE         **deriv2,
	       DOUBLE         **xjm,
	       DOUBLE         **derxy,
	       DOUBLE         **derxy2,
	       DOUBLE         **eveln,
	       DOUBLE         **evelng,
	       DOUBLE          *epren,
	       DOUBLE          *edeadn,
	       DOUBLE          *edeadng,
	       DOUBLE          *velint,
	       DOUBLE          *vel2int,
	       DOUBLE          *covint,
	       DOUBLE         **vderxy,
	       DOUBLE          *pderxy,
	       DOUBLE         **vderxy2,
	       DOUBLE         **wa1,
	       DOUBLE         **wa2
	      );      	      	      	     	     	    	   	   

/************************************************************************
 | f3_caliterrhs.c                                                      |
 ************************************************************************/ 
/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for vel dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the iteration forces for vel dofs
is calculated:

                   /
   (+/-) THETA*dt |  v * u * grad(u)  d_omega
                 /  


see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *eforce    DOUBLE	   (i/o)  element force vector
\param  *covint    DOUBLE	   (i)	  conv. vels at INT. point
\param  *funct     DOUBLE	   (i)    nat. shape funcs
\param   fac 	   DOUBLE	   (i)    weighting factor
\param	 iel	   INT		   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgalifv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		  DOUBLE          *covint,
		  DOUBLE          *funct,
		  DOUBLE           fac,
		  INT              iel
                 );
/*!---------------------------------------------------------------------  
\brief stabilisation part of iteration forces for vel dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the iteration forces for vel dofs
is calculated:

                   /
   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
                 /    

                       /
   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) *  * grad(u)  d_omega
                     /



see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param   *covint   DOUBLE	   (i)    conv. vels at INT. point
\param   *velint   DOUBLE	   (i)    vel at integr. point
\param	 *funct    DOUBLE	   (i)	  nat. shape funcs
\param  **derxy    DOUBLE	   (i)	  global derivative
\param  **derxy2   DOUBLE	   (i)    2nd global derivative
\param    fac 	   DOUBLE	   (i)    weighting factor
\param    visc     DOUBLE	   (i)    fluid viscosity
\param    ihoel    INT		   (i)    flag for higher order ele
\param	  iel	   INT		   (i)    num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifv(
                  FLUID_DYN_CALC  *dynvar, 
                  ELEMENT         *ele,
		  DOUBLE          *eforce,   /* element force vector */
		  DOUBLE          *covint,
		  DOUBLE          *velint,
		  DOUBLE          *funct,
		  DOUBLE         **derxy,
		  DOUBLE         **derxy2,
		  DOUBLE           fac,
		  DOUBLE           visc,
		  INT              ihoel,
		  INT              iel
                 );
/*!--------------------------------------------------------------------- 
\brief stabilisation part of iteration forces for pre dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the iteration forces for pre 
dofs is calculated:

                   /
   (-/+) THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
                 /  

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:							
    there's only one full element force vector  	
    for pre-dofs the pointer eforce points to the entry  
    eforce[3*iel]					
      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *eforce   DOUBLE	   (i/o)  element force vector
\param   *covint   DOUBLE	   (i)    conv. vels at INT. point
\param  **derxy    DOUBLE	   (i)    global derivative
\param    fac      DOUBLE	   (i)    weighting factor
\param	  iel	   INT  	   (i)	  num. of nodes of act. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabifp(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		  DOUBLE          *covint,
		  DOUBLE          **derxy,
		  DOUBLE           fac,
		  INT              iel
                 );
		 		 		 
/************************************************************************
 | f3_calservice.c                                                      |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         genk 05/02

get the element velocities and the pressure at different times 

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the	 
      transformation in every time step 			 
				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param  **eveln    DOUBLE	   (o)    ele vels at time n
\param  **evelng   DOUBLE	   (o)    ele vels at time n+g
\param   *epren    DOUBLE	   (o)    ele pres at time n
\param   *edeadn   DOUBLE          (o)    ele dead load at n (selfweight)
\param   *edeadng  DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *hasext   INT             (o)    flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,
                DOUBLE         **eveln,
	        DOUBLE         **evelng,
	        DOUBLE          *epren,
		DOUBLE          *edeadn,
		DOUBLE          *edeadng,
		INT             *hasext		
	      );
/*!--------------------------------------------------------------------- 
\brief routine to calculate velocities at integration point

<pre>                                                         genk 05/02		 
				      
</pre>
\param   *velint   DOUBLE        (o)   velocities at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_veli(
             DOUBLE  *velint,    
             DOUBLE  *funct,    
	     DOUBLE **evel,     
	     INT      iel       
	    );
/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure at integration point

<pre>                                                         genk 05/02
				      
</pre>
\param  *preint    DOUBLE        (o)   pressure at integration point
\param  *funct     DOUBLE        (i)   shape functions
\param  *epre      DOUBLE        (i)   pressure at element nodes
\param   iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_prei(
             DOUBLE  *preint,    
             DOUBLE  *funct,    
	     DOUBLE  *epre,     
	     INT      iel       
	    );
/*!--------------------------------------------------------------------- 
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 05/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
vderxy[0][2] = Ux,z  
				      
</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_vder(
             DOUBLE **vderxy,     
             DOUBLE **derxy,    
	     DOUBLE **evel,     
	     INT      iel       
	    ) ;
/*!--------------------------------------------------------------------- 
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 05/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
vderxy[0][2] = Ux,z  
				      
</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_vder(
             DOUBLE **vderxy,     
             DOUBLE **derxy,    
	     DOUBLE **evel,     
	     INT      iel       
	    );
/*!---------------------------------------------------------------------  
\brief routine to calculate 2nd velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the 2nd derivatives of the velocity
w.r.t x/y/z are calculated
   vderxy2[0][0] = Ux,xx 
   vderxy2[0][3] = Ux,xy 
   vderxy2[1][4] = Ux,xz 
   vderxy2[2][5] = Ux,yz 
				      
</pre>
\param  **vderxy2  DOUBLE        (o)   2nd velocity derivativs
\param  **derxy2   DOUBLE        (i)   2nd global derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_vder2(
             DOUBLE **vderxy2,   
             DOUBLE **derxy2,   
	     DOUBLE **evel,     
	     INT      iel       
	    );
void f3_pder(
             DOUBLE  *pderxy,    
             DOUBLE **derxy,    
	     DOUBLE  *epre,     
	     INT      iel       
	    );
/*!--------------------------------------------------------------------- 
\brief convective velocities 

<pre>                                                         genk 05/02		 

in this routine the convective velocity is calculated at the 
integration point:
 u * grad(u)
 e.g. 3D: COVx = Ux*Ux,x + Uy*Ux,y + Uz*Ux,z

</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param   *velint   DOUBLE        (i)   velocity at integration point
\param   *covint   DOUBLE        (i)   convective velocity at INT point
\return void                                                                       

------------------------------------------------------------------------*/
void f3_covi(
             DOUBLE **vderxy,    
             DOUBLE  *velint,   
	     DOUBLE  *covint    
	    );
/*!--------------------------------------------------------------------- 
\brief permutation of element force vector 

<pre>                                                         genk 05/02		 

routine to rearrange the entries of the element force vector	  
this is necessary since we would like to use the existing assembly
routines for the RHS						  
hence a splitting of vel- and pre dof in the element force vector 
is not possible any more!!!!					  


</pre>
\param   *eforce   DOUBLE        (i/o) element force vector
\param  **tmp      DOUBLE        (i)   working array
\param    iel	   DOUBLE        (i)   number of nodes in this ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_permeforce(
		   DOUBLE    *eforce,
		   DOUBLE   **tmp,
		   INT        iel		   		   
	          );
/*!--------------------------------------------------------------------- 
\brief permutation of element stiffness matrix

<pre>                                                         genk 05/02		 

routine to add galerkin and stabilisation parts of the elment	   
stiffness matrix and to rearrange its entries!  		   
this is necessary since we would like to use the existing assembly 
routines for the stiffness matrix				   
hence a splitting of vel- and pre dofs is not possible any more!!!!				  

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array		
\param	  iel	  INT		 (i)   number of nodes in ele
\param	 *dynvar  FLUID_DYN_CALC
\return void                                                                       

------------------------------------------------------------------------*/
void f3_permestif(                  
		   DOUBLE         **estif,
		   DOUBLE         **emass,
		   DOUBLE         **tmp,
		   INT              iel,
		   FLUID_DYN_CALC  *dynvar		   		   
	          );
		  		  	    	    	    	    	    	    	      
/************************************************************************
 | f3_calstabmat.c                                                      |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_c * div(v) * div(u)   d_omega
  /

    /
   |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
  /
  
    /
   |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
  /  

    /
   |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
  /

    /
   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) d_omega
  /  

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
  /

    /
   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
  /

  
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Kvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)	   weighting factor	   
\param   visc	   DOUBLE	   (i)	   fluid viscosity
\param   iel	   INT		   (i)	   num. of nodes in ele
\param   ihoel	   INT		   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkvv(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE         **vderxy,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,	 
                    INT              ihoel
                   );
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Kvp

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kvv is calculated:

    /
   |  tau_mu * u_old * grad(v) * grad(p)   d_omega
  /

    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				   
      --> Kvp is stored in estif[(0..(3*iel-1)][(3*iel)..(4*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkvp(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
                   );
/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Mvv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Mvv is calculated:

    /
   |  -/+ tau_mu * u_old * grad(v) * u d_omega
  /
  
    /
   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
  /  
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			    
      --> Mvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *ele	   ELEMENT	   (i)	   actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabmvv(
                    ELEMENT         *ele,
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
    		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
                   );
/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kpv is calculated:

    /
   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
  /
  
    /
   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
  /  

    /
   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				    
      --> Kpv is stored in estif[((3*iel)..(4*iel-1)][0..(3*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)      
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param  *velint    DOUBLE	   (i)     vel. at integr. point
\param **vderxy    DOUBLE	   (i)     global vel. deriv.
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param **derxy2    DOUBLE	   (i)     2nd global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   visc      DOUBLE	   (i)     fluid viscosity
\param   iel	   INT  	   (i)	   num. of nodes in ele
\param   ihoel     INT  	   (i)	   flag for higer ord. ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkpv(
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE          *velint,
		    DOUBLE         **vderxy,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE         **derxy2,
		    DOUBLE           fac,
		    DOUBLE           visc,
		    INT              iel,
		    INT              ihoel		 
                   );
/*!--------------------------------------------------------------------- 
\brief evaluate stabilisaton part of Kpp

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Kpp is calculated:

    /
   |  - tau_mp * grad(q) *grad(p) d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  			     
      --> Kpp is stored in				     
	      estif[((3*iel)..(4*iel-1)][((3*iel)..(4*iel-1)] 
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)   ele stiffness matrix
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT  	   (i)     num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabkpp(
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel		 
                   );
/*!---------------------------------------------------------------------
\brief evaluate stabilisaton part of Mpv

<pre>                                                         genk 05/02

In this routine the stabilisation part of matrix Mpv is calculated:

    /
   |  - tau_mp * grad(q) * u d_omega
  /
        
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
  
NOTE: there's only one elestif  				   
      --> Mpv is stored in estif[((3*iel)..(4*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **emass     DOUBLE	   (i/o)   ele mass matrix
\param  *funct     DOUBLE	   (i)     nat. shape functions
\param **derxy     DOUBLE	   (i)     global derivatives
\param   fac	   DOUBLE	   (i)     weighting factor
\param   iel	   INT		   (i)	   num. of nodes in ele

\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabmpv(
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE         **estif,
		    DOUBLE          *funct,
		    DOUBLE         **derxy,
		    DOUBLE           fac,
		    INT              iel		 
                   );

/************************************************************************
 | f3_calstabpar.c                                                      |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief routine to calculate stability parameter                

<pre>                                                         genk 05/02   
  									 
   ele->e.f3->iadvec: advection stab.					
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
      2 = at every INT pt used for element.-stab.-matrices		
   ele->e.f3->istapc: flag for stabilisation parameter calculation	
      1 = at center of element  					
      2 = at every integration point					
   ele->e.f3->clamb \							
   ele->e.f3->c1     |_>> stabilisation constants (input)		
   ele->e.f3->c2     |  						
   ele->e.f3->c3    /							
   ele->e.f3->istrle: has streamlength to be computed			
   ele->e.f3->iarea: calculation of area length 			
   ele->e.f3->iduring: calculation during INT.-pt.loop  		
   ele->e.f3->itau[0]: flag for tau_mu calc. (-1: before, 1:during)	
   ele->e.f3->itau[1]: flag for tau_mp calc. (-1: before, 1:during)	
   ele->e.f3->itau[2]: flag for tau_c calc. (-1: before, 1:during)	
   ele->e.f3->hk[i]: element sizes (vel / pre / cont)			
   ele->e.f3->idiaxy: has diagonals to be computed			
   dynvar->tau[0]: stability parameter momentum / velocity (tau_mu)	
   dynvar->tau[1]: stability parameter momentum / pressure (tau_mp)	
   dynvar->tau[2]: stability parameter continuity (tau_c)		

</pre>

\param   *ele,        ELEMENT	      (i)    actual element
\param   *dynvar,     FLUID_DYN_CALC  (i/o)
\param   *velint,     DOUBLE	      (i)    vel at center
\param    visc,       DOUBLE	      (i)    viscosity
\param    iel,        INT	      (i)    number of nodes	     
\param	  ntyp,       INT	      (i)    element type
\param	  iflag       INT	      (i)    flag for evaluation
\return void                                                                       

------------------------------------------------------------------------*/ 
void f3_calstabpar(
	            ELEMENT         *ele,      
		    FLUID_DYN_CALC  *dynvar,
		    DOUBLE          *velint,  
		    DOUBLE           visc,    
		    INT              iel,     
		    INT              ntyp,    
		    INT              iflag    
                  );
		  
/************************************************************************
 | f3_caltimerhs.c                                                      |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for vel dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the time forces for vel dofs
is calculated:

      /
   + |  v * u     d_omega
    /  
    
                      /
   (-) (1-THETA)*dt  |  v * u * grad(u)     d_omega
                    /
		  
                      /
   (-) (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
                    /  
		  
                    /
   + (1-THETA)*dt  |  div(v) * p  d_omega
                  /		  		      

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:					     
   in ONESTEP methods: velint  = vel2int = U(n)
   in TWOSTEP methods: velint  = U(n+gamma)  
   in TWOSTEP methods: vel2int = U(n)	     
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *velint      DOUBLE	      (i)    vel. at integr. point
\param   *vel2int     DOUBLE	      (i)    vel. at integr. point
\param   *covint      DOUBLE	      (i)    conv. vel. at integr. p.
\param	 *funct       DOUBLE	      (i)    nat. shape functions      
\param	**derxy       DOUBLE	      (i)    global derivatives
\param	**vderxy      DOUBLE	      (i/    global vel. deriv.
\param	  preint      DOUBLE	      (i)    pres. at integr. point
\param	  visc	      DOUBLE	      (i)    fluid viscosity
\param	  fac	      DOUBLE	      (i)    weighting factor
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgaltfv(
                  FLUID_DYN_CALC  *dynvar, 
                  DOUBLE          *eforce,    
		  DOUBLE          *velint,
		  DOUBLE          *vel2int,
		  DOUBLE          *covint,
		  DOUBLE          *funct,
		  DOUBLE         **derxy,
		  DOUBLE         **vderxy,
		  DOUBLE           preint,
		  DOUBLE           visc,
		  DOUBLE           fac,
		  INT              iel
              )  ;
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for pre dofs

<pre>                                                         genk 05/02

In this routine the galerkin part of the time forces for pre dofs
is calculated:

                    /
   + (1-THETA)*dt  |  q * div(u)  d_omega
                  /

see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:						       
    there's only one full element force vector         
    for pre-dofs the pointer eforce points to the entry
    eforce[3*iel]				        	    
      
</pre>
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *eforce      DOUBLE	      (i/o)  element force vector
\param   *funct       DOUBLE	      (i)    nat. shape functions
\param  **vderxy      DOUBLE	      (i)    global vel. deriv.
\param    fac         DOUBLE	      (i)    weighting factor	     
\param	  iel	      INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calgaltfp(
                   FLUID_DYN_CALC  *dynvar, 
                   DOUBLE          *eforce,    
		   DOUBLE          *funct,
		   DOUBLE         **vderxy,
		   DOUBLE           fac,
		   INT              iel
                  ) ;
/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for vel dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:
		 		                   	  		      
      /
   + |  tau_mu * u * grad(v) * u  d_omega
    /
    
        /
   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
      /
      
                     /
   (-) (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
                   /          
    
                     /
   +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
                   /

                   /
   + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
                 /

                     /
   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
                   /		 
		 	
       /
  (-) |  tau_c * div(v) * div(u)  d_omega
     /
   
                    /
  (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
                  /
		
                    /
  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
                  /
		  
                  /
  + (1-THETA)*dt | tau_mu * u * grad(v) * b    d_omega
                /
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'
      
</pre>
\param   *dynvar     FLUID_DYN_CALC   (i)
\param   *ele        ELEMENT	      (i)    actual element
\param   *eforce     DOUBLE	      (i/o)  element force vector
\param   *velint     DOUBLE	      (i)    vel. at integr. point
\param   *vel2int    DOUBLE	      (i)    vel. at integr. point
\param	 *covint     DOUBLE	      (i)    conv. vel. at integr. p.
\param	**derxy      DOUBLE	      (i)    global derivatives
\param	**derxy2     DOUBLE	      (i)    2nd global derivatives
\param	**vderxy     DOUBLE	      (i)    global vel. deriv.
\param	**vderxy2    DOUBLE	      (i)    2nd global vel. deriv.
\param	 *pderxy     DOUBLE	      (i)    global pressure deriv.
\param	  fac	     DOUBLE	      (i)    weighting factor
\param	  visc	     DOUBLE	      (i)    fluid viscosity
\param	  ihoel      INT	      (i)    flag for higer ord. ele
\param	  iel	     INT	      (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabtfv(
                   FLUID_DYN_CALC  *dynvar, 
                   ELEMENT         *ele,
	           DOUBLE          *eforce,    
	 	   DOUBLE          *velint,
		   DOUBLE          *vel2int,
		   DOUBLE          *covint,
		   DOUBLE         **derxy,
		   DOUBLE         **derxy2,
		   DOUBLE         **vderxy,
		   DOUBLE         **vderxy2,
		   DOUBLE          *pderxy,
		   DOUBLE           fac,
		   DOUBLE           visc,
		   INT              ihoel,
		   INT              iel
                  );
/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for pre dofs

<pre>                                                         genk 05/02

In this routine the stabilisation part of the time forces for vel dofs
is calculated:
		 		                   	  		      
          /
   (-)   |  tau_mp * grad(q) * u  d_omega
        /
      
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
                   /
		   
                     /
   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
                   /		         
				       
see also dissertation of W.A. Wall chapter 4.4 'Navier-Stokes Loeser'

NOTE:						       
    there's only one full element force vector         
    for pre-dofs the pointer eforce points to the entry
    eforce[3*iel]				       
      
</pre>
\param   *dynvar,     FLUID_DYN_CALC   (i)
\param   *eforce,     DOUBLE	       (i/o)  element force vector
\param  **derxy,      DOUBLE	       (i)    global derivatives
\param  **vderxy2,    DOUBLE	       (i)    2nd global vel. deriv.
\param   *velint,     DOUBLE	       (i)    vel. at integr. point
\param	 *covint,     DOUBLE	       (i)    conv. vel. at integr. p.
\param	 *pderxy,     DOUBLE	       (i)    global pressure deriv.
\param	  visc,       DOUBLE	       (i)    fluid viscosity
\param	  fac,        DOUBLE	       (i)    weighting factor
\param	  ihoel,      INT	       (i)    flag for higer ord. ele
\param	  iel	      INT	       (i)    num. of nodes in ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabtfp(
                   FLUID_DYN_CALC  *dynvar, 
                   DOUBLE          *eforce,    
     		   DOUBLE         **derxy,
		   DOUBLE         **vderxy2,
		   DOUBLE          *velint,
		   DOUBLE          *covint,
		   DOUBLE          *pderxy,
		   DOUBLE           visc,
		   DOUBLE           fac,
		   INT              ihoel,
		   INT              iel
                  ) ;
		  		  		  	      
/************************************************************************
 | f3_inpele.c                                                          |
 ************************************************************************/
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
void f3inp(ELEMENT *ele, INT counter);

/************************************************************************
 | f3_intg.c                                                            |
 ************************************************************************/
void f3_intg(FLUID_DATA      *data,
             INT              option);

/************************************************************************
 | f3_main.c                                                            |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief main fluid3 control routine

<pre>                                                         genk 05/02
</pre>
\param  *actpart	 PARTITION    (i)	    
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param	*etforce_global  ARRAY        (o)    element time force vector
\param  *eiforce_global  ARRAY        (o)    element iter force vecotr
\param	*edforce_global  ARRAY        (o)    ele dirichl. force vector
\param	*action	         CALC_ACTION  (i)
\param	*hasdirich	 INT          (o)    flag
\param  *hasext          INT          (o)    flag
\return void

------------------------------------------------------------------------*/
void fluid3(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global, 
	    ARRAY       *etforce_global,
	    ARRAY       *eiforce_global,
	    ARRAY       *edforce_global,
            CALC_ACTION *action,
	    INT         *hasdirich,
	    INT         *hasext,     	    
	    CONTAINER   *container
           );
	   
	    
