/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
#ifdef D_LS


typedef struct _LS2
{
  INT                  ntyp;         /* flag for element type                     */
  INT                  nGP[2];       /* # of GP in r and s directions             */
  INT                  is_elcut;     /* flag to control element activation        */
  INT                  is_sedge_set; /* flag to control element activation        */
  INT                  is_elsearch;  /* flag to confine element search            */
  INT                  nlayer;       /* corresponding layer number                */
  INT                  nenode;       /* number of enriched nodes                  */
  INT                  enode[4];     /* local node numbers of enriched nodes      */
  INT                  intcnt;       /* index to interface number involved        */
  INT                  prncnt;       /* counter used in printing                  */  
  INT                  rstcnt;       /* counter used in resetting                 */

  struct _LS_INT_DATA  *intdata;     /* struct to store interface data            */
  struct _LS_POLY_DATA *polydata;    /* struct to store polygon data              */
  struct _ELEMENT      *my_fluid;    /* pointer to fluid element associated       */
} LS2;



typedef struct _LS_INT_DATA
{
  DOUBLE               p[2][2];      /* start and end  points of interface        */
  DOUBLE               pd[2];        /* intersection with diagonal                */
  INT                  edge[2];      /* start and end edges                       */
  INT                  polycon[2][3];/* polygon connectivity                      */
  INT                  is_diagcut;   /* is diagonal cut                           */    
  INT                  reconstruct;  /* flag for reconstruction                   */
  struct _ELEMENT     *nbor;         /* element neighbor to eedge                 */
  struct _ELEMENT     *nbor_s;       /* element neighbor to sedge                 */  
} LS_INT_DATA;



typedef struct _LS_POLY_DATA
{
  INT                  ind[2];                                    
  INT                  polygonmat[2][7];   
  DOUBLE               polygonwgt[2][7];
  DOUBLE               polygonGP[2][2][7];
} LS_POLY_DATA;



typedef struct _LS_DYNAMIC
{
  INT                     step;        /* the actual step                                 */
  INT 	                  nstep;       /* (n)umber of time steps                          */
  INT 	                  nfstep;      /* (n)umber of fluid only time steps               */
  INT 	                  nfreinit;    /* frequency of reinitialization                   */
  INT  	                  itemax;      /* (m)aximum num. of itn. at a certain time step   */
  INT                     iop;         /* time (i)ntegration scheme                       */
  INT	                  ite;         /* iteration scheme                                */	
  INT                     itchk;       /* convergence check during nonlin. iteration      */
  INT                     stchk;       /* flag for steady state check                     */
  INT                     itnrm;       /* norm for conv. check d. nonlin. iteration       */
  INT 	                  init;	       /* initialization of level set profile             */
  INT 	                  nif;         /* flag for evaluation of time rhs                 */
  INT 	                  nii;         /* flag for evaluation of iteration rhs            */
  DOUBLE                  time;	       /* actual time                                     */
  DOUBLE                  dt;	       /* time increment                                  */
  DOUBLE                  maxtime;     /* (m)aximum time                                  */
  DOUBLE                  ittol;       /* tolerance to check convergence                  */
  DOUBLE                  sttol;       /* tolerance to check convergence                  */
  DOUBLE                  theta;       /* integration parameter for one step theta scheme */
  struct _LS_GEN_DATA    *lsdata;
} LS_DYNAMIC;



typedef struct _LS_GEN_DATA
{
  INT        numcirc;
  INT        numline;
  INT        is_sharp;
  INT        setvel;
  INT        flag_vel;
  /* circle information */
  DOUBLE     xc1;
  DOUBLE     yc1;
  DOUBLE     rad1;
  DOUBLE     xc2;
  DOUBLE     yc2;
  DOUBLE     rad2;
  DOUBLE     xc3;
  DOUBLE     yc3;
  DOUBLE     rad3;
  /* line information */
  DOUBLE     xs1;
  DOUBLE     ys1;
  DOUBLE     xe1;
  DOUBLE     ye1;
  DOUBLE     xs2;
  DOUBLE     ys2;
  DOUBLE     xe2;
  DOUBLE     ye2;
  DOUBLE     xs3;
  DOUBLE     ys3;
  DOUBLE     xe3;
  DOUBLE     ye3;
  
  INT 	     isstab;           /* flag to control stabilization                   */
  INT 	     localization;     /* flag to control localization                    */
  INT 	     numlayer;         /* flag to control localization region             */
  
  INT 	     reinitialization; /* flag to control reinitialization                */
  INT 	     numreinit;        /* number of reinitialization                      */
  INT 	     reinitflag;       /* initialization flag                             */
  DOUBLE     rdt;              /* dt for reinitialization                         */
  INT        algo;             /* treatment of velocity (imp. or exp.)            */
  DOUBLE     epsilon;          /* smoothing parameter for sign function           */      

  INT        boundary_on_off;
  INT        reconstruct_on_off;
} LS_GEN_DATA;



typedef struct _LS_UPDATE
{
  INT 	      *nael;              /* number of active elements                    */
  INT 	       nand;              /* number of active nodes                       */
  ELEMENT     *first_in_list;     /* first element in the list                    */
  INT          numinterface;      /* number of interfaces                         */
  INT         *typeinterface;     /* type of interface (open or closed)           */
} LS_UPDATE;



typedef struct _LS2_INTG_DATA
{
  /* for quadrilateral element  */
  DOUBLE     wgtq[3][3];
  DOUBLE     xgq[3][3];
  /* for triangular element  */
  DOUBLE     wgtt[3][3];
  DOUBLE     xgtr[3][3];
  DOUBLE     xgts[3][3];
} LS2_INTG_DATA;



typedef enum _LSFLAG
{
  ls_initphase,
  ls_updtphase,
  ls_finaphase
} LSFLAG;



typedef enum _FRONTLSFLAG
{
  front_ls_set,
  front_ls_init,
  front_ls_initmat,    
  front_ls_updt,
  front_ls_updtmat,  
  front_ls_activate,
  front_ls_localize,
  front_ls_polygonize,
  front_ls_write,
  front_ls_finalize
} FRONTLSFLAG;
#endif
