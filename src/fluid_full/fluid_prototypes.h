/*----------------------------------------------------------------------*
 |  routine to control fluid dynamic analyis                 genk  03/02|
 *----------------------------------------------------------------------*/
void dyn_fluid();

/*----------------------------------------------------------------------*
 | routine to control implicit and semi-implicit algorithms for fluid   |
 | problems combined with Newton and fixed point iteration schemes.     |                                                  
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_isi(FLUID_DYNAMIC *fdyn);

/*----------------------------------------------------------------------*
 | routine to control stationary algorithm for fluid for fluid problems |
 | combined with Newton, fixed point iteration and fixed point like     |
 | schemes.                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_stat();

/*----------------------------------------------------------------------*
 | routine to control predictor-multicorrector algorithm for fluid      |
 | problems                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_pm();

/*----------------------------------------------------------------------*
 | routine to calculate time integration constants                      |
 | for fluid time algorithms                                genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_tcons(FLUID_DYNAMIC *fdyn,
                FLUID_DYN_CALC *dynvar);

/*----------------------------------------------------------------------*
 | routine to calculate constants for nonlinear iteration   genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_icons(FLUID_DYNAMIC *fdyn,
                FLUID_DYN_CALC *dynvar,
		int itnum);		

/*----------------------------------------------------------------------*
 | input of initial starting field from fluid_start.data                |
 *----------------------------------------------------------------------*/
void inp_fluid_start_data();


/*----------------------------------------------------------------------*
 | routine to initialise solution history                    genk 04/02 |
 *----------------------------------------------------------------------*/
void fluid_init(
		FIELD  *actfield,  
                FLUID_DYNAMIC *fdyn	
	       );    

/*----------------------------------------------------------------------*
 | routine to determine the absolute maximum velocity for use as global |
 | scaling velocity in some stability parameter definitions (Tau_mp)    |
 |                                                          genk  03/02 |
 | at the moment only sequentiell version                               |
 *----------------------------------------------------------------------*/
void fluid_maxvel(
                  FLUID_DYNAMIC *fdyn, 
                  FIELD  *actfield,               		 
		  int actpos,
		  double *res
		 )       ; 	
/*----------------------------------------------------------------------*
 | routine to extract digits from integer number                        |
 |                                                          genk 04/02  |
 *----------------------------------------------------------------------*/
void intextract(
                int num,    /* integer number */
                int *it,    /* integer on position "thousand" */
		int *ih,    /* integer on position "hundred" */
		int *id,    /* integer on position "ten" */
		int *io     /* integer on position "one" */
	       );		 
/*----------------------------------------------------------------------*
 |  routine to initialise the dirichlet boundary conditions             |
 *----------------------------------------------------------------------*/
void fluid_initdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn);		 	  
 /*----------------------------------------------------------------------*
 |  routine to set dirichlet boundary conditions on at time <time>       |
 |                                                           genk 04/02  |
 *----------------------------------------------------------------------*/
void fluid_setdirich(FIELD  *actfield, FLUID_DYNAMIC *fdyn);
/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_increment                                |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Calculates the norms for the iteration check if necessary           |
 |                                                          genk 05/02  |
 *----------------------------------------------------------------------*/
void fluid_result_incre(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ,
			  double *vrat, double *prat, FLUID_DYNAMIC *fdyn );

/*----------------------------------------------------------------------*
 | calculate vel- and pre-norms for steady state check                  |
 | norm = ||U(n+1) - U(n)|| / ||U(n)||                      genk 05/02  |
 |    solution at (n+1): node->sol_increment[3][j]                      |
 |    solution at (n)  : node->sol_increment[1][j]                      |
 *----------------------------------------------------------------------*/
void fluid_norm( FLUID_DYNAMIC *fdyn, FIELD  *actfield, int numeq_total,
                 double *vrat, double *prat );	
/*----------------------------------------------------------------------*
 | routine to copy solution history                         genk 05/02  |
 |    solution at (n+1): node->sol_increment[3][j]                      |
 |    solution at (n)  : node->sol_increment[1][j]                      |
 *----------------------------------------------------------------------*/
void fluid_copysol(FLUID_DYNAMIC *fdyn, FIELD  *actfield, 
                   int from, int to, int flag);	
/*----------------------------------------------------------------------*
 | routine to do the steady state check                     genk 05/02  |
 *----------------------------------------------------------------------*/
int fluid_steadycheck(FLUID_DYNAMIC *fdyn, FIELD  *actfield, int numeq_total);
/*----------------------------------------------------------------------*
 | routine to do the iteration convergence check            genk 05/02  |
 *----------------------------------------------------------------------*/
int fluid_convcheck(FLUID_DYNAMIC *fdyn, double vrat, double prat, 
                    int itnum,double te, double ts);	 
/*----------------------------------------------------------------------*
 | routine to print out time and algo to the screen         genk 05/02  |
 *----------------------------------------------------------------------*/
void fluid_algoout(FLUID_DYNAMIC *fdyn, FLUID_DYN_CALC *dynvar);		      			  
