#ifdef WALLCONTACT
/*!------------------------------------------------------------------------
\file
\brief wallcontact.h headerfile for 2-D bilinear contact interfaces, containing
structures and prototypes 

------------------------------------------------------------------------*/
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief main structure for 2-D contact (bilinear discretization)

<pre>                                                           m.gee 10/02
defined in wall_contact_detection.c                             
</pre>

-------------------------------------------------------------------------*/

typedef enum _CONTACTFLAG
{
      contact_off,  /*!< Contact flag indicating no contact */
      contact_on    /*!< Contact flag indicating contact */
} CONTACTFLAG;


typedef struct _WALL_CONTACT
{
int                 ndline;        /*!< Number of desing lines labelled as master or slave */
struct _DLINE     **dline;         /*!< Master and slave desing lines are stored in these structures*/
int                 ng_masterline; /*!< Number of master g_lines*/
int                 ng_slaveline;  /*!< Number of slave g_lines*/
struct _GLINE     **g_masterline;  /*!< Master g_lines are stored in these structures*/
struct _GLINE     **g_slaveline;   /*!< Slave g_lines are stored in these structures*/
int                 ng_masternode; /*!< Number of master g_nodes*/
int                 ng_slavenode;  /*!< Number of slave g_nodes*/
struct _GNODE     **g_masternode;  /*!< Master g_nodes are stored in these structures*/
struct _GNODE     **g_slavenode;   /*!< Slave g_nodes are stored in these structures*/
struct _GNODE     **contact_set;   /*!< Active contact g_nodes are stored in this structure*/
int                 set_size;      /*!< Size of the active set (g_nodes which are in contact)*/
int                 CET_flag;      /*!< Flag indicating Cons. Enforcement Tech. used (0: Penalty Met. 1: Aug. Lagr. Meth.)*/
int                 FR_flag;       /*!< Flag indicating Contact Pr. Type (0: Frictionless 1: Frictional)*/  
double              n_pen_par;     /*!< Normal Penalty Parameter*/
double              t_pen_par;     /*!< Tangential Penalty Parameter*/
double              fr_coef;       /*!< Coefficient of friction*/
#ifdef GEMM
double              dt;            /*!<Time step size used in case of E-M int. scheme with contact*/
#endif
enum _CONTACTFLAG   contactflag;   /*!< Contact flag*/ 

} WALL_CONTACT;



typedef enum _CONTACTTYPE
{
      contact_none,     /*!< not used */
      contact_master,   /*!< the object is of master type */
      contact_slave,    /*!< the object is of slave type */
      contact_self      /*!< not used */
} CONTACTTYPE;



typedef struct _HISTORY  /*!< History structure used to keep track of some variables*/
{
double           cr_local_coord;     /*!< Current closest point projection (local coordinate)*/
double           pr_local_coord;     /*!< Previous(time step) closest point projection (local coordinate)*/
struct _GNODE   *pr_masters[2];      /*!< G_nodes of the previous master segment */
struct _NODE    *pr_closest;         /*!< Previous(time step) closest master node to the slave node*/
double           pr_multipliers[2];  /*!< Lagrange Multipliers */
CONTACTFLAG      pr_flag;            /*!< Contact flag of the slave node indicating the state of contact at the previous step*/
double           cr_g;               /*!< Current value of gap function */
double           R_Metric;           /*!< Reference Metric of the closest point */
double           cr_force;           /*!< Current normal traction component*/
double           cr_tan;             /*!< Current tangential traction component*/
double           pr_t_tan;           /*!< Previous tangential traction component*/
#ifdef GEMM
double           tau_n[3];           /*!< Previous tau_n*/
double           g_n;                /*!< Previous value of gap function*/
double           mid_projection;     /*!< Closest point projection at mid-configuration*/
double           mid_metric;         /*!< Metric of the closest point at the mid configuration*/
double           g_mid;              /*!< Value of gap function at the mid-configuration*/
double           g_tilda;            /*!< Algorithmic gpa rate*/
double           mid_velocity[6];    /*!< Relative velocity of slave node - master segment obtained at the mid-configuration*/
double           mid_normal[3];      /*!< Outward normal at the closest point obtained at mid-configuration*/
double           mid_tangent[3];     /*!< Tangential vector at the closest point obtained at mid-configuration*/
#endif
} HISTORY;

/*! @} (documentation module close)*/
#endif
