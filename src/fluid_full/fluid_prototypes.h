/*!----------------------------------------------------------------------
\file
\brief fluid_prototypes

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/   
/************************************************************************
 | fluid_dirich.c                                                       |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief routine to initialise the dirichlet boundary conditions

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are initialised:
- real pressure from input-file is transformed to kinematic pressure:
  actgnode->dirich->dirich_val.a.dv[predof] /= dens

- if the initial field is a zero-field, then set the dirichlet values in
  the solution history of the node:
    'actnode->sol.a.da[0][j] = initval*acttimefac' (--> output)
    'actnode->sol_increment.a.da[1][j] = initval*acttimefac'
                                 |        |         |
                            time (n)      |         |               
                          initial value from input  |               
                                       factor from timecurve (T=0.0)			     
</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)    

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_initdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn);
/*!---------------------------------------------------------------------                                         
\brief routine to set dirichlet boundary conditions at time <time>

<pre>                                                         genk 04/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->time>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[3][j] = initval*acttimefac'
                                 |        |         |
                            time (n+1)    |         |               
                          initial value from input  |               
                                       factor from timecurve			     
</pre>
\param *actfield FIELD         (i)  actual field (fluid)   
\param *fdyn	 FLUID_DYNAMIC (i)  

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_setdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn);
/*!---------------------------------------------------------------------                                         
\brief routine to calculate the element dirichlet load vector

<pre>                                                         genk 04/02

in this routine the element load vector due to dirichlet conditions
is calcluated. The prescribed values are taken from the node solution
history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[3][j]'.
the element load vector 'dforce' is calculated by eveluating
</pre>
\code
      dforces[i] -= estif[i][j] * dirich[j];
\endcode			     

\param  *actele    ELEMENT   (i)   actual element	  
\param  *dforces   double    (o)   dirichlet force vector
\param **estif     double    (i)   element stiffness matrix
\param  *hasdirich int       (o)   flag if s.th. was written to dforces

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_caldirich(
                     ELEMENT   *actele,  
		     double    *dforces, 
                     double   **estif,   
		     int       *hasdirich
		    ); 

/************************************************************************
 | fluid_dyn.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief routine to control fluid dynamic analyis

<pre>                                                         genk 03/02

In this routine the different control programs for fluid-problems are
called. This depends on the paremeter TIMEINTEGR, which is stored in
fdyn->iop
iop=0: Stationary Solution
iop=1: Predictor Multicorrector scheme
iop=2: Semi-Implicit-One-Step Method
iop=3: Semi-Implicit-Two-Step Method
iop=4: One-Step-Theta Scheme
iop=5: Fractional-Step-Theta Scheme 

see dissertation of W.A. WALL, chapter 4.2 Zeitdiskretisierung;		     
</pre>


\return void                                                                             

------------------------------------------------------------------------*/
void dyn_fluid(void);

/************************************************************************
 | fluid_imp_semimp.c                                                   |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief implicit and semi-implicit algorithms for fluid problems

<pre>                                                         genk 03/02

this routine conrols the implicit and semi-implicit algorithms for fluid
problems combined with different nonlinear iteration schemes.

time-discretisation:
fdyn->iop=2: Semi-Implicit-One-Step Method
fdyn->iop=3: Semi-Implicit-Two-Step Method
fdyn->iop=4: One-Step-Theta Scheme
fdyn->iop=5: Fractional-Step-Theta Scheme 

see dissertation of W.A. Wall chapter 4.2 'Zeitdiskretisierung'

nonlinear iteration scheme:
fdyn->ite=0: no nonlinear iteration
fdyn->ite=1: fixed-point-like iteration
fdyn->ite=2: Newton iteration
fdyn->ite=3: fixed-point iteration

see dissertation chapter 4.3 'Linearisierung und Iteratonsverfahren'.
			     
</pre>   
\param *fdyn	 FLUID_DYNAMIC (i)    

\return void 
\warning up to now only the One-Step-Theta scheme combined with a
         fixed-point-like iteration scheme is tested! 

------------------------------------------------------------------------*/
void fluid_isi(FLUID_DYNAMIC *fdyn);

/************************************************************************
 | fluid_pr_mcorr.c                                                     |
 ************************************************************************/
 /*---------------------------------------------------------------------*
 | routine to control predictor-multicorrector algorithm for fluid      |
 | problems                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_pm(void);

/************************************************************************
 | fluid_service.c                                                      |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief routine to check starting algorithm

<pre>                                                         genk 04/02

this routine conrols the starting algorithms schemes. For the first 'nums'
iteration steps a different time integration scheme may be used then 
during the rest of the simulation.
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC (i)  
\param *nfrastep 	int           (o)  number of fract. steps

\return void 
\warning this routine is not completely tested yet!

------------------------------------------------------------------------*/
void fluid_startproc(
                     FLUID_DYNAMIC  *fdyn,
		     int            *nfrastep 
		    );
/*!---------------------------------------------------------------------                                         
\brief calculating time integration constants

<pre>                                                         genk 03/02

in this routine the constants for the time integration algorithms are 
calculated
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC   (i)  
\param *dynvar	        FLUID_DYN_CALC  (i/o)  

\return void 
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_tcons(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar);
/*!---------------------------------------------------------------------                                         
\brief setting flags for nonlinear iteration schemes

<pre>                                                         genk 03/02

in this routine the flags for the different nonlinear iteration schemes
are set. Depending on the iteration schemes the evaluation of some 
termes in the LHS or the RHS have to turned on or off:

nik <->  EVALUATION OF LHS-MATRICES (w/o NONLINEAR TERM)    
nic <->  EVALUATION OF NONLINEAR LHS N-CONVECTIVE	    
nir <->  EVALUATION OF NONLINEAR LHS N-REACTION 	    
nie <->  EVALUATE ONLY LHS-TERMS FOR EXPLICIT VELOCITY      
nil <->  EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)      
nif <->  EVALUATION OF "TIME - RHS" (F-hat)		    
nii <->  EVALUATION OF "ITERATION - RHS"		    
nis <->  STATIONARY CASE (NO TIMEDEPENDENT TERMS)	   
			     
</pre>   
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param *dynvar	   FLUID_DYN_CALC  (i/o)  
\param  itnum      int  	   (i)     actual number of iterations
\return void 
\warning up to now, only fixed-point like iteration checked!!!

------------------------------------------------------------------------*/
void fluid_icons(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar,
		 int itnum           
		);
/*!---------------------------------------------------------------------                                         
\brief initialisation of solution history

<pre>                                                         genk 04/02

in this routine the solution history of the fluid field is initialised.
The time-integration schemes require the solutions at different time-
steps. They are stored in 	   
node->sol_incement: solution history used for calculations	    
      sol_increment[0][i]: solution at (n-1)			    
      sol_increment[1][i]: solution at (n)			    
      sol_increment[2][i]: solution at (n+g)			    
      sol_increment[3][i]: solution at (n+1)			    
In node->sol one findes the converged solutions of the time-steps, which
are written to the output-file and the pss-file.

If the initial fluid fuild comes from the input-file 'fluid-start-data',
these data are copied to sol_increment.
			     
</pre>   
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param *dynvar	   FLUID_DYN_CALC  (i/o)  
\param  itnum      int  	   (i)     actual number of iterations
\return void 
\warning up to now, only fixed-point like iteration checked!!!

------------------------------------------------------------------------*/
void fluid_init(
		FIELD  *actfield,  
                FLUID_DYNAMIC *fdyn	
	       );
/*!---------------------------------------------------------------------                                         
\brief calculating norms for steady state check

<pre>                                                         genk 05/02

in this routine the velocity and pressure norms for the steady state 
check are calculated:
   norm = ||U(n+1) - U(n)|| / ||U(n)||  		
      solution at (n+1): node->sol_increment[3][j]			
      solution at (n)  : node->sol_increment[1][j]			
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)   actual field
\param  numeq_total   int	     (i)   total number of equations
\param *vrat	      double	     (o)   vel.  conv. ratio
\param *prat	      double	     (o)   pres. conv. ratio
\return void 

------------------------------------------------------------------------*/
void fluid_norm(FLUID_DYNAMIC *fdyn, 	     
                FIELD         *actfield,    
		int            numeq_total, 
                double        *vrat,        
		double        *prat         
	       );		       	       				 		     
/*!---------------------------------------------------------------------                                         
\brief copy solution history

<pre>                                                         genk 05/02

in this routine the solution at postion 'from' in the nodal solution 
history is copied to the positon 'to'.	
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)  actual field
\param  from          int	     (i)  pos. in sol(_increment)
\param  to   	      int	     (i)  pos. in sol(_increment)
\param  flag 	      int	     (i)  simply a flag
\return void 
\sa fluid_init()

------------------------------------------------------------------------*/
void fluid_copysol(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,  
                   int            from,     
		   int            to,       
		   int            flag      
		  );
/*!---------------------------------------------------------------------                                         
\brief steady state check

<pre>                                                         genk 05/02

in this routine the convergence ratios for the steady state check are 
calculated and the result is printed to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)  actual field
\param  numeq_total   int	     (i)  total number of equations
\return int steady  

------------------------------------------------------------------------*/
int fluid_steadycheck(FLUID_DYNAMIC *fdyn, 	  
                      FIELD         *actfield,   
		      int            numeq_total 
		     );
/*!---------------------------------------------------------------------                                         
\brief iteration convergence check

<pre>                                                         genk 05/02

in this routine the iteration convergence ratios are compared with
the given tolerance. The result is printed out to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param  vrat          double  	     (i)  vel. conv. ratio
\param  prat          double         (i)  pres. conv. ratio
\param	itnum 	      int	     (i)  actual numb. of iter steps
\param	te 	      double	     (i)  time for element calcul.
\param	ts	      double	     (i)  solver time
\return int converged  

------------------------------------------------------------------------*/
int fluid_convcheck(FLUID_DYNAMIC *fdyn,   
                    double         vrat,  
		    double         prat,  
                    int            itnum, 
		    double         te,    
		    double         ts     
		   );
/*!---------------------------------------------------------------------                                         
\brief print out to the screen

<pre>                                                         genk 05/02

time-integration parameters are printed out to the screen
			     
</pre>   
\param *fdyn 	        FLUID_DYNAMIC   (i)   
\param *dynvar	        FLUID_DYN_CALC  (i/o) 
\return void 

------------------------------------------------------------------------*/
void fluid_algoout(FLUID_DYNAMIC  *fdyn, 
                   FLUID_DYN_CALC *dynvar
		  );		   		     		  
	       
/************************************************************************
 | fluid_stationary.c                                                   |
 ************************************************************************/
/*----------------------------------------------------------------------*
 | routine to control stationary algorithm for fluid for fluid problems |
 | combined with Newton, fixed point iteration and fixed point like     |
 | schemes.                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_stat(void);

/************************************************************************
 | fluid_visual.c                                                       |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief output to pss-file

<pre>                                                         genk 07/02

in this routine the solution (node->sol) is written to the pss-file of
proc 0 used for visualisation with VISUAL2 / VISUAL3
all other data required by VISUAL2 / VISUAL3 are read from the input    
file
			     
</pre>   
\param *actfield	FIELD         (i)  actual field
\param  ntsteps         int           (i)  total num. of time steps
\param *time_a          ARRAY         (i)  time array 
\return void  

------------------------------------------------------------------------*/
void fluid_writevispss(FIELD  *actfield,  
                       int     ntsteps,  
		       ARRAY  *time_a	 
		      );
/*!---------------------------------------------------------------------                                         
\brief input from pss-file

<pre>                                                         genk 07/02

in this routine the solution (node->sol) is read from the pss-file of
proc 0 used for visualisation with VISUAL2 / VISUAL3
all other data required by VISUAL2 / VISUAL3 are read from the input    
file
			     
</pre>   
\param *actfield	FIELD         (i)  actual field
\param  ntsteps         int           (i)  total num. of time steps
\param *time_a          ARRAY         (i)  time array 
\return void  

------------------------------------------------------------------------*/
void fluid_readvispss(FIELD   *actfield, 
                      int     *ntsteps,  
		      ARRAY   *time_a	 
		     );

/************************************************************************
 | inp_fluid_start_data.c                                               |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief input from file 'fluid_start-data'

<pre>                                                         genk 07/02

in this routine the inital fluid data are read form 'fluid_start_data'
and stored in the array 'start'
			     
</pre>   
\return void  
\warning ONLY SEQUENTIEL VERSION !!!!!
------------------------------------------------------------------------*/
void inp_fluid_start_data(void);
   
