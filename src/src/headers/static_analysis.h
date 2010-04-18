/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | general static-variables                               m.gee 6/01    |
 | variables used by linear or nonlinear structural static analysis     |
 *----------------------------------------------------------------------*/
typedef struct _STATIC_VAR
{
INT                 linear;             /* is linear calculation */
INT                 nonlinear;          /* is nonlinear calculation */
INT                 multiscale;         /* flag if multiscale model */
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
INT                 graderw;            /* is gradient damage model */
INT                 praedictor;         /* is it praedictor -> 1 first corrector: ->2 else: 3 */


struct _NODE       *controlnode;        /* ptr to control node */
INT                 control_node_global;/* global control node Id (redundant) */
INT                 control_dof;        /* dof of control node to be controlled */

INT                 isrelstepsize;      /* is stepsize constant or variable */
INT                 actstep[20];        /* number of steps for actual variable stepsize */
DOUBLE              actstepsize[20];    /* actual variable stepsize */
INT                 numcurve;           /* number of entries in actstep and actstepsize */

INT                 reldisnode_ID[6];     /* ID's of relative-displacement-output-nodes*/
INT                 reldis_dof[6];        /* dof we want in output for relative displacement */
#ifdef D_MLSTRUCT
DOUBLE              eps_equiv;            /* strain measure, for which multiscale is activated */
#endif
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
