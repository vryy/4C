/*!----------------------------------------------------------------------
\file
\brief headerfile for 1D interface element, 
       containing structures

*----------------------------------------------------------------------*/
#ifdef D_INTERF

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief 1D interface element

<pre>                                                            ah 08/03
This structure contains the working array at the Gaussian points for the 
1D interface element
</pre>

*----------------------------------------------------------------------*/

typedef struct _IF_IP_WA
{
INT          yip;         /* displacement jump normal  */              
DOUBLE       jump_ut_pl;  /* plastic displacement jump tangential  */              
DOUBLE       Q[2][2];     /* material tangent  */              
DOUBLE       dn;          /* damage normal  */              
DOUBLE       dt;          /* damage tangential  */              
DOUBLE       Tn;          /* stress normal  */              
DOUBLE       Tt;          /* stress tangential  */              
} IF_IP_WA;

/*!----------------------------------------------------------------------
\brief 1D interface element

<pre>                                                            ah 08/03
This structure contains the working arrays for the 
1D interface element
</pre>

*----------------------------------------------------------------------*/
typedef struct _IF_ELE_WA
{
IF_IP_WA      *ipwa;
/*DOUBLE       *matdata; *//* element material data, actual density ... */
} IF_ELE_WA;

/*!----------------------------------------------------------------------
\brief 1D interface element

<pre>                                                            ah 08/03
This structure contains all specific information for the 
1D interface element
</pre>

*----------------------------------------------------------------------*/
typedef struct _INTERF
{
INT           nGP;   
DOUBLE        thick;              /*!< thickness perpend. to wallplane */
IF_ELE_WA     *elewa;         
enum
    {
    if_xy,                         /*!< stresses in global orientaion */
    if_tn                          /*!< stresses in local orientaion */
    } stresstyp; 
struct _ARRAY4D  stress_GP;
struct _ARRAY4D  stress_ND;
} INTERF;

/*!----------------------------------------------------------------------
\brief interface element data

<pre>                                                            ah 08/03
This structure contains the coordinates and weights for numerical
integration of the 1D interface element
</pre>

*----------------------------------------------------------------------*/
typedef struct _INTERF_DATA
{
DOUBLE        xgr[3];  /*!< natural coordinates xsi of gaussian points */
DOUBLE        wgtr[3]; /*!< weights at natural coordinates xsi         */
} INTERF_DATA;


/*----------------------------------------------------------------------*
 |  PROTOTYPES OF ALL ALE ROUTINES                          mn 05/03    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | main interface control routine                           mn 05/03    |
 *----------------------------------------------------------------------*/
void interf(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY       *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container); 
         
/*----------------------------------------------------------------------*
 | read input of interface element                         mn 05/03    |
 *----------------------------------------------------------------------*/
void interf_inp(ELEMENT *ele);
/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/


