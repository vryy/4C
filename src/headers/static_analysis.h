/*----------------------------------------------------------------------*
 | general static-variables                               m.gee 6/01    |
 | variables used by linear or nonlinear structural static analysis     |
 *----------------------------------------------------------------------*/
typedef struct _STATIC_VAR               
{
int                 geolinear;          /* is linear calculation */
int                 geononlinear;       /* is nonlinear calculation */
enum _NR_CONTROLTYP nr_controltyp;      /* type of control */
int                 nstep;              /* number of steps */
int                 maxiter;            /* max number of iterations in NR */
double              tolresid;           /* tolerance of residual forces */
double              toldisp;            /* tolerance of residual displacements */
double              stepsize;           /* steplenght */
int                 iarc;               /* flag for arscaling in Crisfields Arclenght control */
double              arcscl;             /* arc scaling scaling factor of load part of predictor */
int                 signchcsp;          /* flag for singn changing by CSP */
int                 resevry_disp;       /* write result very updevry step */
int                 resevry_stress;     /* write result very updevry step */
int                 resevery_restart;   /* write restart every res_write_evry step */
				

struct _NODE       *controlnode;        /* ptr to control node */
int                 control_node_global;/* global control node Id (redundant) */
int                 control_dof;        /* dof of control node to be controlled */
} STATIC_VAR;

/*----------------------------------------------------------------------*
 | general static-control-variables                       m.gee 6/01    |
 | variables to perform Newton Raphson                                  |
 *----------------------------------------------------------------------*/
typedef struct _STANLN  
{
double              sp1;                  /* initial stiffness of control node */
double              csp;                  /* current stiffness parameter */
double              rlold;                /* load factor of last step */
double              rlnew;                /* load factor of actual step */
double              rlpre;                /* load factor from predictor */

double              renorm;               /* some norms */
double              rinorm;
double              rrnorm;

double              renergy;

struct _ARRAY       arcfac;               /* vector of load factors of increments */
} STANLN;
