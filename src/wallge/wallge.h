/*!----------------------------------------------------------------------
\file
\brief headerfile for gradient enhance wall element,
       containing structures

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALLGE

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief gradient enhanced wall element

<pre>                                                            ah 08/03
This structure contains the working array at the Gaussian points for the
gradient enhanced wall element
</pre>

*----------------------------------------------------------------------*/

typedef struct _WGE_IP_WA
{
DOUBLE       sig[4];      /* global stresses  */
INT          yip;         /* flag  */
DOUBLE       kappa;       /* max. equivalent strain, ever reached  */
DOUBLE       damage;      /* damage variable -> for output  */
DOUBLE       aequistrain;      /* local equivalent strain  */
DOUBLE       aequistrain_nl;      /* nonlocal equivalent strain  */
} WGE_IP_WA;

/*!----------------------------------------------------------------------
\brief gradient enhanced wall element

<pre>                                                            ah 08/03
This structure contains the working arrays for the
gradient enhanced wall element
</pre>

*----------------------------------------------------------------------*/
typedef struct _WGE_ELE_WA
{
WGE_IP_WA      *iptwa;
DOUBLE         *matdata; /* element material data */
} WGE_ELE_WA;

/*!----------------------------------------------------------------------
\brief gradient enhanced wall element

<pre>                                                            ah 08/03
This structure contains information about the wall type: plane strain...
</pre>

*----------------------------------------------------------------------*/
typedef enum _WALLGE_TYPE
{
                       pl_strain,
                       pl_stress,
                       rot_symmet
} WALLGE_TYPE;
/*!----------------------------------------------------------------------
\brief gradient enhanced wall element

<pre>                                                            ah 08/03
This structure contains all specific information for the
gradient enhanced wall element
</pre>

*----------------------------------------------------------------------*/
typedef struct _WALLGE
{
enum _WALLGE_TYPE    wgetype;       /*!< type of 2D problem, see above */
INT           nGP[4];
DOUBLE        thick;              /*!< thickness perpend. to wallplane */
WGE_ELE_WA   *elwa;
enum
    {
    wge_xy,                         /*!< stresses in global orientaion */
    wge_rs                          /*!< stresses in local orientaion */
    } stresstyp;
#define WALLGE_STRESSTYPE { "wge_xy", "wge_rs", NULL }
struct _ARRAY4D  stress_GP;
struct _ARRAY4D  stress_ND;
} WALLGE;

/*!----------------------------------------------------------------------
\brief wallge element data

<pre>                                                            ah 08/03
This structure contains the coordinates and weights for numerical
integration of the gradient enhanced wall element
</pre>

*----------------------------------------------------------------------*/

typedef struct _WALLGE_DATA
{
DOUBLE        xgrr[13];         /*!< natural coordinates r of gaussian points */
DOUBLE        wgtr[13];         /*!< weights at natural coordinates r */

DOUBLE        xgss[13];         /*!< natural coordinates s of gaussian points */
DOUBLE        wgts[13];         /*!< weights at natural coordinates s */
} WALLGE_DATA;


/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                          mn 05/03    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | main wallge control routine                           mn 05/03    |
 *----------------------------------------------------------------------*/
void wallge(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container);
/*----------------------------------------------------------------------*
 | read input of wallge element                         mn 05/03    |
 *----------------------------------------------------------------------*/
void wge_inp(ELEMENT *ele);
/*----------------------------------------------------------------------*/

#endif /*D_WALLGE*/
/*! @} (documentation module close)*/


