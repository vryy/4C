/*!----------------------------------------------------------------------
\file
\brief headerfile for 3D hex element, containing structures and prototypes

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief 3D hex element data

<pre>                                                             al 06/01
This structure contains the coordinates and weights for numerical
integration of a 3D hex element
</pre>

*----------------------------------------------------------------------*/
typedef struct _C1_DATA
{
double        xgrr[13]; /*!< natural coordinates r of gaussian points */
double        wgtr[13]; /*!< weights at natural coordinates r */

double        xgss[13]; /*!< natural coordinates s of gaussian points */
double        wgts[13]; /*!< weights at natural coordinates s */

double        xgtt[13]; /*!< natural coordinates t of gaussian points */
double        wgtt[13]; /*!< weights at natural coordinates t */
} C1_DATA;
/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains eas data 
</pre>

*----------------------------------------------------------------------*/
typedef struct _C1_EASDAT
{
/* eas ...*/
double     *disl;  /*!< element displacements (last step) */
double     *ehdis; /*!< element displacements             */ 
double     *hil;   /*!< element parameters                */ 
double     *hih;   /*!< element parameters                */ 
} C1_EASDAT;

/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains the working array for g.p. values 
</pre>

*----------------------------------------------------------------------*/
typedef struct _C1_IP_WA
{
/* mises ...*/
double      sig[6]; /*!< global stresses                              */
double      eps[9]; /*!< global strains                               */
double       epstn; /*!< equivalent strain                            */
int          yip;   /*!< stress state: 1=elastic 2=plastic            */
/* hashin delamination (plasticity) ...*/
double      kappa;  /*!< damage threshold value                       */ 
int         imod;   /*!< flag, indicates if sigy has been reached     */
} C1_IP_WA;


/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains the working array for element values 
</pre>

*----------------------------------------------------------------------*/
typedef struct _C1_ELE_WA
{
/* eas ...*/
C1_EASDAT     *eas;     /*!< eas specific data                          */
C1_IP_WA      *ipwa;    /*!< working array for integration points       */
double        *matdata; /*!< element material data, actual density ...  */
int           *optdata; /*!< optimization variable number ...           */
} C1_ELE_WA;
/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains gp and nodal coordinates and stress values 
</pre>

*----------------------------------------------------------------------*/
typedef struct _C1_ELE_STRESS
{
/*
gpstrs[gp][ 0.. 7] stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz
gpstrs[gp][ 6..11] stress-rr stress-ss stress-tt stress-rs stress-st stress-tr
gpstrs[gp][12..14] stress-11 stress-22 stress-33
gpstrs[gp][15..23] ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3
*/
double      **gpstrs;  /*!< stress values of gaussian points   */
double      **gpcoor;  /*!< coordinates xyz of gaussian points */
/*
npstrs[nn] [ 0.. 2]: stress-rr   stress-ss   stress-tt
npstrs[nn] [ 3.. 5]: stress-rs   stress-st   stress-tr  
npstrs[nn] [ 6.. 8]: stress-xx   stress-yy   stress-zz
npstrs[nn] [ 9..11]: stress-xy   stress-yz   stress-xz 
npstrs[nn] [12..14]: stress-11   stress-22   stress-33
npstrs[nn] [15..17]: ang-r1  ang-s1  ang-t1
npstrs[nn] [18..20]: ang-r2  ang-s2  ang-t2
npstrs[nn] [21..23]: ang-r3  ang-s3  ang-t3
*/
double      **npstrs;  /*!< stress values of nodal points      */
} C1_ELE_STRESS;
/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains all specific information for a 3D hex element
</pre>

*----------------------------------------------------------------------*/
typedef struct _BRICK1
{
int           nGP[3];    /*!< number of gaussian points in rst direction*/
int           nhyb;      /*!< flag whether to use the eas formulation   */
int           form;      /*!< ==2: T.L.  */

C1_ELE_WA     *elewa;                        /*!< element working array */

C1_ELE_STRESS *stress;                       /*!< element stresses      */
/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              fh 6/02
This structure contains arrays for stress values for a 3D hex element
</pre>

*----------------------------------------------------------------------*/
struct _ARRAY4D  stress_GP;  /*!< array for stress values - unused      */
struct _ARRAY4D  stress_ND;  /*!< array for stress values - unused      */
} BRICK1;
/*----------------------------------------------------------------------*
 |  prototypes of main and input routines                   al 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | c1_main.c                                                 al 6/01    |
 | main brick1  control routine                                         |
 *----------------------------------------------------------------------*/
void brick1(PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,
            ARRAY       *estif_global,
            ARRAY       *emass_global,
            ARRAY      *intforce_global,
            CALC_ACTION *action,
            CONTAINER   *container);   /*!< contains variables defined in container.h */
/*----------------------------------------------------------------------*
 | c1_inpele.c                                               al 6/01    |
 | read brick1 element                                                  |
 *----------------------------------------------------------------------*/
void c1inp(ELEMENT *ele);
#endif
/*! @} (documentation module close)*/
