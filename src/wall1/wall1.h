/*----------------------------------------------------------------------*
 | wall1 data                                                al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_DATA
{
double        xgrr[13];
double        wgtr[13];

double        xgss[13];
double        wgts[13];
} W1_DATA;
/*----------------------------------------------------------------------*
 | working array                                             al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_IP_WA
{
double      sig[4]; /* global stresses                              */
double      eps[4]; /* global strains                               */
double         *qn; /* 'backstress' vec                             */
double       epstn; /* equivalent strain                            */
int          yip;   /* stress state: 1=elastic 2=plastic            */
/* concrete */
double       *sigc; /* condensed global stresses                    */
double       *grad; /* actual gradient                              */
double       *dlam; /* uniaxial plastic strain  (- dlam)            */
double       *sigi; /* stresses from last iteration step            */          
double       *epsi; /* strains from last iteration step             */    
double       *di  ; /* components d41,d42,d43,d44 of the            */
                    /* constitutive tensor from last iteration step */
/* rebar */    
double       *rsig;   /* stress */ 
double       *reps;   /* strain */
double       *repstn; /* equivalent strain */
int          *ryip;   /* flag */
/* tension stiffening */
double       stcappae;  /* max. elastic strain of the concrete     */              
double       stcappaut; /* fracture strain of the concrete         */          
double       stalpha;   /* angle of the principal concrete stress  */         
double       stthick ;  /* thickness ot the structure              */        
double       stfbd;     /* tension stiffening stress               */


} W1_IP_WA;
typedef struct _W1_ELE_WA
{
W1_IP_WA      *ipwa;
double         dia;
} W1_ELE_WA;
/*----------------------------------------------------------------------*
 | result struct: stresses                                   al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_ELE_STRESS
{
double      *gprr;
double      *gpss;
double      *gprs;

double      *fps; /* first  principal stress       */                         
double      *sps; /* second principal stress       */                    
double      *aps; /* angle of principal direction  */                    

double      *ferr;
double      *fess;
double      *fers;
} W1_ELE_STRESS;
/*----------------------------------------------------------------------*
 | type of 2D problem                                       al 01/02    |
 *----------------------------------------------------------------------*/
typedef enum _WALL_TYPE
{
                       plane_strain, 
                       plane_stress,
                       rotat_symmet  
} WALL_TYPE;
/*----------------------------------------------------------------------*
 | wall                                                     al 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _WALL1
{
enum _WALL_TYPE    wtype;             /* type of 2D problem, see above */

int           nGP[4];

double        thick;

W1_ELE_WA     *elewa;                         /* element working array */

W1_ELE_STRESS *stress;                        /* element stresses      */

} WALL1;
/*----------------------------------------------------------------------*
 |  prototypes of main and input routines                   al 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | main wall1  control routine                               al 6/01    |
 *----------------------------------------------------------------------*/
void wall1(       PARTITION   *actpart,
                  INTRA       *actintra,
                  ELEMENT     *ele,
                  ARRAY       *estif_global,
                  ARRAY       *emass_global,
                  double      *global_vec,
                  int          global_numeq,
                  CALC_ACTION *action);
/*----------------------------------------------------------------------*
 | read wall element                                         al 9/01    |
 *----------------------------------------------------------------------*/
void w1inp(ELEMENT *ele);
