/*----------------------------------------------------------------------*
 | general static-variables                               m.gee 6/01    |
 | variables used by linear or nonlinear structural static analysis     |
 *----------------------------------------------------------------------*/
typedef struct _STATIC_VAR               
{
INT                 linear;             /* is linear calculation */
INT                 nonlinear;          /* is nonlinear calculation */
enum _KINTYP        kintyp;             /* type of kinematic used in nonlinear calculation*/
enum _NR_CONTROLTYP nr_controltyp;      /* type of control */
INT                 nstep;              /* number of steps */
INT                 maxiter;            /* max number of iterations in NR */
DOUBLE              tolresid;           /* tolerance of residual forces */
DOUBLE              toldisp;            /* tolerance of residual displacements */
DOUBLE              stepsize;           /* steplenght */
INT                 iarc;               /* flag for arscaling in Crisfields Arclenght control */
DOUBLE              arcscl;             /* arc scaling scaling factor of load part of predictor */
INT                 signchcsp;          /* flag for singn changing by CSP */
INT                 resevry_disp;       /* write result very updevry step */
INT                 resevry_stress;     /* write result very updevry step */
INT                 resevery_restart;   /* write restart every res_write_evry step */
				

struct _NODE       *controlnode;        /* ptr to control node */
INT                 control_node_global;/* global control node Id (redundant) */
INT                 control_dof;        /* dof of control node to be controlled */
} STATIC_VAR;

/*----------------------------------------------------------------------*
 | general static-control-variables                       m.gee 6/01    |
 | variables to perform Newton Raphson                                  |
 *----------------------------------------------------------------------*/
typedef struct _STANLN  
{
DOUBLE              sp1;                  /* initial stiffness of control node */
DOUBLE              csp;                  /* current stiffness parameter */
DOUBLE              rlold;                /* load factor of last step */
DOUBLE              rlnew;                /* load factor of actual step */
DOUBLE              rlpre;                /* load factor from predictor */

DOUBLE              renorm;               /* some norms */
DOUBLE              rinorm;
DOUBLE              rrnorm;

DOUBLE              renergy;

struct _ARRAY       arcfac;               /* vector of load factors of increments */
} STANLN;
