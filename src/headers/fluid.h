/*!---------------------------------------------------------------------
\file
\brief fluid dynamic data

------------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief integration parameters                                              

<pre>                                                         genk 03/02  

In this structure the coordinates and weights used by gauss integration
are stored.

</pre>

------------------------------------------------------------------------*/
typedef struct _FLUID_DATA
{
double        qxg[MAXQINTP][MAXQINTC];    /*!< coordinates for QUADS and HEX */
double        qwgt[MAXQINTP][MAXQINTC];   /*!< weights for QUADS and HEX */

double        txgr[MAXTINTP][MAXTINTC];   /*!< coordinates in r for TRIS and TETS */
double        txgs[MAXTINTP][MAXTINTC];   /*!< coordinates in s for TRIS and TETS*/
double        txgt[MAXTINTP][MAXTINTC];   /*!< coordinates in t for TRIS and TETS*/
double        twgt[MAXTINTP][MAXTINTC];   /*!< weights for TRIS and TETS*/
} FLUID_DATA;


/*!----------------------------------------------------------------------
\brief calculation parameter                                              

<pre>                                                         genk 03/02  

In this structure all parameters used during the element evaluation are
stored.

</pre>

------------------------------------------------------------------------*/
typedef struct _FLUID_DYN_CALC                
{
double dta;      /*!< actual time increment */
double thsl;     /*!< theta-s,l: const. for "stiffness" terms LHS */
double thsr;     /*!< theta-s,r: const. for "stiffness" terms RHS */
double thpl;     /*!< theta-p,l: const. for "pressure" terms LHS  */
double thpr;     /*!< theta-p,r: const. for "pressure" terms RHS  */
double omt;      /*!< ONE-theta                                   */
double acttime;  /*!< actual time */
double velmax;   /*!< max. velocity, needed for stabilisaton parameter */
double tau[3];   /*!< array for stabilitity parameter */
double tau_tu;   /*!< array for stabilitity parameter for turbulence*/
double tau_tu_dc;/*!< array for DISCONTINUITY CAPTURING for turbulence*/
double sigma;    /*!< const. for nonlinear iteration */   
double washvel;  /*!< wall shear velocity */   
double totarea;  /*!< total area of fluid field */
double coord_scale[2];  /*!<coordinates for scaling the turbulence variables */   
int    itwost;   /*!< control variable for element evaluation */
int    isemim;   /*!< control variable for element evaluation */
int    iprerhs;  /*!< treatment of pressure in time discr. */
int    surftens; /*!< include surface tension effects */
int    fsstnif; 
int    fsstnii;
int    nik; 	 /*!< EVALUATION OF LHS-MATRICES (w/o NONLINEAR TERM)*/
int    nic;	 /*!< EVALUATION OF NONLINEAR LHS N-CONVECTIVE       */
int    nir;	 /*!< EVALUATION OF NONLINEAR LHS N-REACTION	     */
int    nie;	 /*!< EVALUATE ONLY LHS-TERMS FOR EXPLICIT VELOCITY  */
int    nil;	 /*!< EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)  */
int    nif;	 /*!< EVALUATION OF "TIME - RHS"           	     */
int    nii;	 /*!< EVALUATION OF "ITERATION - RHS"		     */
int    nis;	 /*!< STATIONARY CASE (NO TIMEDEPENDENT TERMS)       */
int    niturbu_pro;  /*!< EVALUATION OF "TIME - RHS" for turbulence-model */
int    niturbu_n;    /*!< EVALUATION OF "TIME - RHS" for turbulence-model */
int    kapeps_flag;  /*!< kappa or epsilon equation                       */
int    kapomega_flag;/*!< kappa or omega equation                         */
int    kappan;       /*!< kappan for production-term                      */
int    dis_capt;     /*!< flag for DISCONTINUITY CAPTURING for turbulence model */
int    ishape;   /*!< flag for new element shape                     */
int    ncols;        /*!< number of columns in solution history */
struct  _FLUID_DATA data;
} FLUID_DYN_CALC;

/*!----------------------------------------------------------------------
\brief fluid input parameter                                              

<pre>                                                         genk 03/02  

In this structure all fluid-dynamic variables from the input file are
stored.

</pre>

------------------------------------------------------------------------*/
typedef struct _FLUID_DYNAMIC               
{
char               dyntyp[50];   /*!< dynamictype */
int                numdf;        /*!< number of dofs of the fluid elements */
int                iop;          /* !<time integration method */
int                numcont;      /*!< number of continuation steps */
int                uppss;        /*!< update pss file every n steps */
int                upout;        /*!< store results every n steps */      
int                nstep;        /*!< number of timesteps */
int                step;         /*!< the actual step */
int                stepke;       /*!< the actual step for kappa-epsilon*/
int                ite;          /*!< nonlinear iteration scheme */
int                itemax;       /*!< number of nonlin. iterations */
int                itemax_ke;    /*!< number of nonlin. iterations for kappa-eps */
int                itchk;        /*!< convergence check during nonlin. iteration */
int                itnorm;       /*!< norm for conv. check d. nonlin. iteration */
int                stchk;        /*!< steady state check every n steps */
int                stnorm;       /*!< norm for steady state check */
int                iops;         /*!< starting algorithm */
int                nums;         /*!< number of starting algorithm steps */
int                init;         /*!< initialisation of starting field */
int                iprerhs;      /*!< treatment of pressure in time discr. */
int                viscstr;      /*!< flag for calculation of viscos stresses */
int                freesurf;     /*!< treatment of free surface */
int                surftens;     /*!< include surface tension effects */
int                checkarea;    /*!< check total area of fluid field */
int                turbu;        /*!< the type of turbulence-model */
int                dis_capt;     /*!< flag for DISCONTINUITY CAPTURING for turbulence model */
double             lenght;       /*!< internal lenght of problem */
double             rought;       /*!< roughtness of solid boundaries */
double      coord_scale[2];      /*!< coordinates for scaling the turbulence variables */   
double             maxtime;      /*!< maximal simulation time */
double             time;         /*!< actual time */
double             dt;           /*!< time increment */
double             alpha;        /*!< time integration constant */
double             theta;        /*!< time integration constant */
double             gamma;        /*!< time integration constant */
double             ittol;        /*!< tolerance for iteration convergence check */
double             sttol;        /*!< tolerance for steady state check */
double             thetas;       /*!< constant for starting algorithm) */
struct _ARRAY      start;        /*!< starting field */
struct _FLUID_DYN_CALC dynvar;
} FLUID_DYNAMIC;

