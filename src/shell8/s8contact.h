/*!---------------------------------------------------------------------
\file
\brief structures for shell contact

---------------------------------------------------------------------*/
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/
#ifdef S8CONTACT
/*!------------------------------------------------------------------------
\brief main structure all kinds of contact stuff is kept in for shell8

m.gee 2/03  

main structure all kinds of contact stuff is kept in   

-------------------------------------------------------------------------*/
typedef struct _SHELLCONTACT
{
int                       numnp;      /*!< number of contact nodes */
struct _SHELLNODE        *cnode;      /*!< vector of contact nodes */
/* searching stuff for nested sector algorithm from HALQUIST 1990 */
double                    maxdiag;    /*!< maximum diagonal of an element */
double                    buckdim[3]; /*!< the dimension of a searching bucket */
double                    min[3];     /*!< minimum coordinates */
int                       nbuck[3];   /*!< number of uckets in each dimension */
struct _CONTACTSLICE     *slice;
} SHELLCONTACT;
/*!------------------------------------------------------------------------
\brief flag indicating contact

m.gee 2/03  

flag indicating contact

-------------------------------------------------------------------------*/
typedef enum _SHELLCONTACTFLAG
{
   s8_c_off,
   s8_c_on
} SHELLCONTACTFLAG;
/*!------------------------------------------------------------------------
\brief flag indicating which side of an element a slave is projecting on

m.gee 2/03  

flag indicating which side of an element a slave is projecting on

-------------------------------------------------------------------------*/
typedef enum _SHELLPROJECTON
{
   s8_c_project_none,
   s8_c_project_ontop,
   s8_c_project_onbot
} SHELLPROJECTON;
/*----------------------------------------------------------------------*/
#define EPSN (1.0E+03)
#define EPST (1.0E+02)
#define NU   (0.5)
/*----------------------------------------------------------------------*/
/*!------------------------------------------------------------------------
\brief one special contact node for shell contact

m.gee 2/03  

one special contact node for shell contact

-------------------------------------------------------------------------*/
typedef struct _SHELLNODE
{
struct _NODE             *node;     /*!< ptr to node */
double                    xr[6];    /*!< reference coodinates of point including director of true length */            
double                    xc[6];    /*!< current   coodinates of point including director of true length */    
/* geometry history */
double                    xc_his[6]; /*! current configuration of the last converged step */

enum _SHELLCONTACTFLAG    topflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      topproj;  /*!< flag to indicate projection onto a top or bottom surface */
ARRAY                     forcetop; /*!< the contact force */
ARRAY                     stifftop; /*!< the contact stiffness */
ARRAY                     lmtop;    /*!< location matrix for stiffness contributions */
ELEMENT                  *topele;   /*!< the element this node is contacting to, otherwise NULL */
double                    xitop[2]; /*!< local coordinates of projection point */
double                    topgap;   /*!< gap function value of the top of this node */
double                    top_tn;
double                    top_ln;
double                    top_lt[2];
double                    top_tT[2];
/* history of the top */
enum _SHELLCONTACTFLAG    histopflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      histopproj;  /*!< flag to indicate projection onto a top or bottom surface */
ELEMENT                  *histopele;   /*!< the element this node is contacting to, otherwise NULL */
double                    hisxitop[2]; /*!< local coordinates of projection point */
double                    histopgap;   /*!< gap function value of the top of this node */
double                    his_top_lt[2];
double                    his_top_tT[2];


enum _SHELLCONTACTFLAG    botflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      botproj;  /*!< flag to indicate projection onto a top or bottom surface */
ARRAY                     forcebot; /*!< the contact force */
ARRAY                     stiffbot; /*!< the contact stiffness */
ARRAY                     lmbot;    /*!< location matrix for stiffness contributions */
ELEMENT                  *botele;   /*!< the element this node is contacting to, otherwise NULL */
double                    xibot[2]; /*!< local coordinates of projection point */
double                    botgap;   /*!< gap function value of the bot of this node */
double                    bot_tn;
double                    bot_ln;
double                    bot_lt[2];
double                    bot_tT[2];
/* history of the bot */
enum _SHELLCONTACTFLAG    hisbotflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      hisbotproj;  /*!< flag to indicate projection onto a top or bottom surface */
ELEMENT                  *hisbotele;   /*!< the element this node is contacting to, otherwise NULL */
double                    hisxibot[2]; /*!< local coordinates of projection point */
double                    hisbotgap;   /*!< gap function value of the bot of this node */
double                    his_bot_lt[2];
double                    his_bot_tT[2];

int                       nneigh;
struct _NODE             *neighbours[30];
struct _CONTACTBUCKET    *buck;
} SHELLNODE;
/*!------------------------------------------------------------------------
\brief restart of one special contact node for shell contact

m.gee 3/03  

one special contact node for shell contact

-------------------------------------------------------------------------*/
typedef struct _SHELLNODERES
{
/* geometry history */
double                    xc_his[6]; /*! current configuration of the last converged step */

enum _SHELLCONTACTFLAG    topflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      topproj;  /*!< flag to indicate projection onto a top or bottom surface */
double                    top_ln;
double                    top_lt[2];
/* history of the top */
enum _SHELLCONTACTFLAG    histopflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      histopproj;  /*!< flag to indicate projection onto a top or bottom surface */

enum _SHELLCONTACTFLAG    botflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      botproj;  /*!< flag to indicate projection onto a top or bottom surface */
double                    bot_ln;
double                    bot_lt[2];
/* history of the bot */
enum _SHELLCONTACTFLAG    hisbotflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      hisbotproj;  /*!< flag to indicate projection onto a top or bottom surface */
} SHELLNODERES;
/*!------------------------------------------------------------------------
\brief structure to restart the contact stuff

m.gee 3/03  

structure to restart the contact stuff

-------------------------------------------------------------------------*/
typedef struct _RESTARTCONTACT
{
int                        numnp;       /*!< number of nodes in field */     
long int                   mainhandle;
long int                  *handles;     /*!< the handles for the cnodes */
} RESTARTCONTACT;
/*!------------------------------------------------------------------------
\brief a slice of the bucket search algorithm

m.gee 3/03  

slice of the bucket search algorithm

-------------------------------------------------------------------------*/
typedef struct _CONTACTSLICE
{
struct _CONTACTSTRIPE    *stripe;
} CONTACTSLICE;
/*!------------------------------------------------------------------------
\brief a stripe of the bucket search algorithm

m.gee 3/03  

slice of the bucket search algorithm

-------------------------------------------------------------------------*/
typedef struct _CONTACTSTRIPE
{
struct _CONTACTBUCKET    *buck;
} CONTACTSTRIPE;
/*!------------------------------------------------------------------------
\brief a bucket of the bucket search algorithm

m.gee 3/03  

slice of the bucket search algorithm

-------------------------------------------------------------------------*/
typedef struct _CONTACTBUCKET
{
int                       ncnode;        /*!< number of cnode in this bucket */    
struct _SHELLNODE       **cnode;         /*!< ptrs to cnodes in this bucket */
int                       ijk[3];        /*!< ijk indizes of this bucket */
} CONTACTBUCKET;





/* prototypes */
/*----------------------------------------------------------------------*
 |  s8_contact2.c                                         m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_make(SHELLNODE *actcnode,
                     ELEMENT   *actele,
                     double    *xi,
                     int        ssurf,
                     int        msurf);
/*----------------------------------------------------------------------*
 |  s8_contact_gap.c                                      m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_gapfunction(SHELLNODE  *actcnode,
                            int        *ssurf,
                            int        *msurf,
                            ELEMENT    *actele,
                            double      xi[],
                            double     *g);
/*----------------------------------------------------------------------*
 |  s8_contact_project.c                                  m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_orthproject(SHELLNODE  *actcnode,
                            int        *ssurf,
                            int        *msurf,
                            ELEMENT    *actele,
                            double      xires[],
                            double     *distance,
                            int        *success);
void s8_contact_functderiv(double     funct[], 
                           double    deriv[][4], 
                           double    deriv2[],
                           double      r, 
                           double      s);
void s8_contact_metrics(double x[][4],
                        double a3[][4],
                        double e3,
                        double gkov[][3],
                        double gkon[][3],
                        double gmkov[][3],
                        double gmkon[][3],
                        double funct[],
                        double deriv[][4],
                        int    iel);
void s8_contact_inv3(double a[][3], double *det);
void s8_contact_trans(double a[][3], int n);
void s8_contact_deta(double gkov[][3], double *deta);
/*----------------------------------------------------------------------*
 |  s8_contact_init.c                                     m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8contact_init(FIELD *actfield, PARTITION* actpart, INTRA *actintra);
/*----------------------------------------------------------------------*
 |  s8_contact_nearestnode.c                              m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_searchupdate(INTRA *actintra, double dt);
void s8_contact_nearestnode(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            int        *msurftop,
                            int        *msurfbot,
                            double     *distop,
                            double     *disbot);
void s8_contact_nearestnode_bruteforce(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            int        *msurftop,
                            int        *msurfbot,
                            double     *distop,
                            double     *disbot);
/*----------------------------------------------------------------------*
 |  s8_contact_history.c                                  m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_history(INTRA *actintra);
void s8_contact_restartwrite(INTRA *actintra,int step);
void s8_contact_restartread(INTRA *actintra,int step);
void s8_contact_setlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra);
void s8_contact_updlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra);
/*! @} (documentation module close)*/
#endif
