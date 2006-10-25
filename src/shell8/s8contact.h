/*!----------------------------------------------------------------------
\file
\brief structures for shell contact

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef S8CONTACT

/*!
\addtogroup CONTACTS8
*//*! @{ (documentation module open)*/


/*!------------------------------------------------------------------------
\brief main structure all kinds of contact stuff is kept in for shell8

m.gee 2/03

main structure all kinds of contact stuff is kept in

-------------------------------------------------------------------------*/
typedef struct _SHELLCONTACT
{
INT                       numnp;      /*!< number of contact nodes */
struct _SHELLNODE        *cnode;      /*!< vector of contact nodes */
/* searching stuff for nested sector algorithm from HALQUIST 1990 */
DOUBLE                    maxdiag;    /*!< maximum diagonal of an element */
DOUBLE                    buckdim[3]; /*!< the dimension of a searching bucket */
DOUBLE                    min[3];     /*!< minimum coordinates */
INT                       nbuck[3];   /*!< number of uckets in each dimension */
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
#define EPSN (1.0E+04)
#define EPST (1.0E+02)
#define NU   (0.4)
/*----------------------------------------------------------------------*/
/*!------------------------------------------------------------------------
\brief one special contact node for shell contact

m.gee 2/03

one special contact node for shell contact

-------------------------------------------------------------------------*/
typedef struct _SHELLNODE
{
struct _NODE             *node;     /*!< ptr to node */
DOUBLE                    xr[6];    /*!< reference coodinates of point including director of true length */
DOUBLE                    xc[6];    /*!< current   coodinates of point including director of true length */
/* geometry history */
DOUBLE                    xc_his[6]; /*! current configuration of the last converged step */

enum _SHELLCONTACTFLAG    topflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      topproj;  /*!< flag to indicate projection onto a top or bottom surface */
ARRAY                     forcetop; /*!< the contact force */
ARRAY                     stifftop; /*!< the contact stiffness */
ARRAY                     lmtop;    /*!< location matrix for stiffness contributions */
ELEMENT                  *topele;   /*!< the element this node is contacting to, otherwise NULL */
DOUBLE                    xitop[2]; /*!< local coordinates of projection point */
DOUBLE                    topgap;   /*!< gap function value of the top of this node */
DOUBLE                    top_tn;   /*!< contact normal force */
DOUBLE                    top_ln;   /*!< augmented langrangian normal contact parameter (== contact force) */
DOUBLE                    top_lt[2];/*!< augmented langrangian normal frictional parameter (== frictional force) */
DOUBLE                    top_tT[2];/*!< contact frictional force */


/* history of the top */
enum _SHELLCONTACTFLAG    histopflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      histopproj;  /*!< flag to indicate projection onto a top or bottom surface */
ELEMENT                  *histopele;   /*!< the element this node is contacting to, otherwise NULL */
DOUBLE                    hisxitop[2]; /*!< local coordinates of projection point */
DOUBLE                    oldhisxitop[2]; /*!< local coordinates of projection point */
DOUBLE                    histopgap;   /*!< gap function value of the top of this node */
DOUBLE                    his_top_ln;
DOUBLE                    his_top_lt[2];


enum _SHELLCONTACTFLAG    botflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      botproj;  /*!< flag to indicate projection onto a top or bottom surface */
ARRAY                     forcebot; /*!< the contact force */
ARRAY                     stiffbot; /*!< the contact stiffness */
ARRAY                     lmbot;    /*!< location matrix for stiffness contributions */
ELEMENT                  *botele;   /*!< the element this node is contacting to, otherwise NULL */
DOUBLE                    xibot[2]; /*!< local coordinates of projection point */
DOUBLE                    botgap;   /*!< gap function value of the bot of this node */
DOUBLE                    bot_tn;
DOUBLE                    bot_ln;
DOUBLE                    bot_lt[2];
DOUBLE                    bot_tT[2];
/* history of the bot */
enum _SHELLCONTACTFLAG    hisbotflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      hisbotproj;  /*!< flag to indicate projection onto a top or bottom surface */
ELEMENT                  *hisbotele;   /*!< the element this node is contacting to, otherwise NULL */
DOUBLE                    hisxibot[2]; /*!< local coordinates of projection point */
DOUBLE                    oldhisxibot[2]; /*!< local coordinates of projection point */
DOUBLE                    hisbotgap;   /*!< gap function value of the bot of this node */
DOUBLE                    his_bot_ln;
DOUBLE                    his_bot_lt[2];

INT                       nneigh;
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
DOUBLE                    xc_his[6]; /*! current configuration of the last converged step */

enum _SHELLCONTACTFLAG    topflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      topproj;  /*!< flag to indicate projection onto a top or bottom surface */
DOUBLE                    top_ln;
DOUBLE                    top_lt[2];
/* history of the top */
enum _SHELLCONTACTFLAG    histopflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      histopproj;  /*!< flag to indicate projection onto a top or bottom surface */

enum _SHELLCONTACTFLAG    botflag;  /*!< flag to indicate active contact */
enum _SHELLPROJECTON      botproj;  /*!< flag to indicate projection onto a top or bottom surface */
DOUBLE                    bot_ln;
DOUBLE                    bot_lt[2];
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
INT                        numnp;       /*!< number of nodes in field */
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
INT                       ncnode;        /*!< number of cnode in this bucket */
struct _SHELLNODE       **cnode;         /*!< ptrs to cnodes in this bucket */
INT                       ijk[3];        /*!< ijk indizes of this bucket */
} CONTACTBUCKET;





/* prototypes */
/*----------------------------------------------------------------------*
 |  s8_contact2.c                                         m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_make(SHELLNODE *actcnode,
                     ELEMENT   *actele,
                     DOUBLE    *xi,
                     INT        ssurf,
                     INT        msurf);
/*----------------------------------------------------------------------*
 |  s8_contact_gap.c                                      m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_gapfunction(SHELLNODE  *actcnode,
                            INT        *ssurf,
                            INT        *msurf,
                            ELEMENT    *actele,
                            DOUBLE      xi[],
                            DOUBLE     *g);
/*----------------------------------------------------------------------*
 |  s8_contact_project.c                                  m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_orthproject(SHELLNODE  *actcnode,
                            INT        *ssurf,
                            INT        *msurf,
                            ELEMENT    *actele,
                            DOUBLE      xires[],
                            DOUBLE     *distance,
                            INT        *success,
                            DOUBLE     *nue);
void s8_contact_functderiv(DOUBLE     funct[],
                           DOUBLE    deriv[][4],
                           DOUBLE    deriv2[],
                           DOUBLE      r,
                           DOUBLE      s);
void s8_contact_metrics(DOUBLE x[][4],
                        DOUBLE a3[][4],
                        DOUBLE e3,
                        DOUBLE gkov[][3],
                        DOUBLE gkon[][3],
                        DOUBLE gmkov[][3],
                        DOUBLE gmkon[][3],
                        DOUBLE funct[],
                        DOUBLE deriv[][4],
                        INT    iel);
void s8_contact_inv3(DOUBLE a[][3], DOUBLE *det);
void s8_contact_trans(DOUBLE a[][3], INT n);
void s8_contact_deta(DOUBLE gkov[][3], DOUBLE *deta);
void s8_contact_timeguess(SHELLNODE *actcnode,
                         ELEMENT   *actele,
                         DOUBLE    *xi,
                         DOUBLE     distance,
                         DOUBLE    *nue,
                         DOUBLE    *dt);
/*----------------------------------------------------------------------*
 |  s8_contact_init.c                                     m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8contact_init(FIELD *actfield, PARTITION* actpart, INTRA *actintra);
/*----------------------------------------------------------------------*
 |  s8_contact_nearestnode.c                              m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_searchupdate(INTRA *actintra, DOUBLE dt);
void s8_contact_nearestnode(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            INT        *msurftop,
                            INT        *msurfbot,
                            DOUBLE     *distop,
                            DOUBLE     *disbot);
void s8_contact_nearestnode_bruteforce(SHELLNODE  *actcnode,
                            SHELLNODE **nearcnodetop,
                            SHELLNODE **nearcnodebot,
                            INT        *msurftop,
                            INT        *msurfbot,
                            DOUBLE     *distop,
                            DOUBLE     *disbot);
/*----------------------------------------------------------------------*
 |  s8_contact_history.c                                  m.gee 3/03    |
 *----------------------------------------------------------------------*/
void s8_contact_history(INTRA *actintra);
void s8_contact_restartwrite(INTRA *actintra,INT step);
void s8_contact_restartread(INTRA *actintra,INT step);
void s8_contact_setlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra);
void s8_contact_updlagr(FIELD *actfield, PARTITION *actpart, INTRA *actintra);
void s8_contact_historyback(INTRA *actintra);

/*! @} (documentation module close)*/
#endif
