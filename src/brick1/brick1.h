/*!----------------------------------------------------------------------
\file
\brief headerfile for 3D hex element, containing structures and prototypes

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1
#include "brick1_doxygen.h"
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
DOUBLE        xgrr[13]; /*!< natural coordinates r of gaussian points */
DOUBLE        wgtr[13]; /*!< weights at natural coordinates r */

DOUBLE        xgss[13]; /*!< natural coordinates s of gaussian points */
DOUBLE        wgts[13]; /*!< weights at natural coordinates s */

DOUBLE        xgtt[13]; /*!< natural coordinates t of gaussian points */
DOUBLE        wgtt[13]; /*!< weights at natural coordinates t */
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
DOUBLE     *disl;  /*!< element displacements (last step) */
DOUBLE     *ehdis; /*!< element displacements             */
DOUBLE     *hil;   /*!< element parameters                */
DOUBLE     *hih;   /*!< element parameters                */
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
DOUBLE      sig[6]; /*!< global stresses                              */
DOUBLE      eps[9]; /*!< global strains                               */
DOUBLE       epstn; /*!< equivalent strain                            */
INT          yip;   /*!< stress state: 1=elastic 2=plastic            */
/* hashin delamination (plasticity) ...*/
DOUBLE      kappa;  /*!< damage threshold value                       */
INT         imod;   /*!< flag, indicates if sigy has been reached     */
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
DOUBLE        *matdata; /*!< element material data, actual density ...  */
INT           *optdata; /*!< optimization variable number ...           */
} C1_ELE_WA;
/*!----------------------------------------------------------------------
\brief 3D hex element

<pre>                                                              al 06/01
This structure contains all specific information for a 3D hex element
</pre>

*----------------------------------------------------------------------*/
typedef struct _BRICK1
{
INT           nGP[3];    /*!< number of gaussian points in rst direction*/
INT           nhyb;      /*!< flag whether to use the eas formulation   */
INT           form;      /*!< ==2: T.L.  */

C1_ELE_WA     *elewa;                        /*!< element working array */
/*----------------------------- stresses in local and global systems ---*/
/*
stress_GP[place][ 0.. 5][gp] stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz
stress_GP[place][ 6..11][gp] stress-rr stress-ss stress-tt stress-rs stress-st stress-tr
stress_GP[place][12..14][gp] stress-11 stress-22 stress-33
stress_GP[place][15..23][gp] ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3
stress_GP[place][24..26][gp] x-coord.    y-coord.    z-coord.

stress_ND[place][ 0.. 2][nn] stress-rr   stress-ss   stress-tt
stress_ND[place][ 3.. 5][nn] stress-rs   stress-st   stress-tr
stress_ND[place][ 6.. 8][nn] stress-xx   stress-yy   stress-zz
stress_ND[place][ 9..11][nn] stress-xy   stress-yz   stress-xz
stress_ND[place][12..14][nn] stress-11   stress-22   stress-33
stress_ND[place][15..17][nn] ang-r1  ang-s1  ang-t1
stress_ND[place][18..20][nn] ang-r2  ang-s2  ang-t2
stress_ND[place][21..23][nn] ang-r3  ang-s3  ang-t3
*/
enum
    {
    c1_nostr,  /*!<  default                           */
    c1_gpxyz,  /*!<  gp: global system on element level    */
    c1_gprst,  /*!<  gp: local  system on element level    */
    c1_gp123,  /*!<  gp: principal-stresses                */
    c1_npxyz,  /*!<  extrapolation to nodes: global system on element level    */
    c1_nprst,  /*!<  extrapolation to nodes: local  system on element level    */
    c1_np123,  /*!<  extrapolation to nodes: principal-stresses                */
    c1_npeqs   /*!<  extrapolation to nodes: equivalent stress                 */
    }            stresstyp;
struct _ARRAY4D  stress_GP;  /*!< array for stress values */
struct _ARRAY4D  stress_ND;  /*!< array for stress values */
} BRICK1;
/*!----------------------------------------------------------------------
\brief rst-flag for edge, surface - load evaluation

<pre>                                                              al 06/01
</pre>
*----------------------------------------------------------------------*/
typedef enum _RSTF
{ /* r..r-dir., s..s-dir., t..t-dir., n..-1.0, p..+1.0 */
  rpp, rpn, rnp, rnn, psp, psn, nsp, nsn, ppt, pnt, npt, nnt,/*line*/
  rsn, rsp, rnt, rpt, nst, pst,/*surface*/
  rst/*volume*/
} RSTF;
/*----------------------------------------------------------------------*
 |  prototypes of main and input/output routines            al 11/01    |
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
            CONTAINER   *container);
/*----------------------------------------------------------------------*
 | c1_inpele.c                                               al 6/01    |
 | read brick1 element                                                  |
 *----------------------------------------------------------------------*/
void c1inp(ELEMENT *ele);
/*----------------------------------------------------------------------*
 |  routine to calcualate and to write gp stresses of a brick1          |
 |  element to visualize in gid -> Hexahedra elements        al 1/03    |
 *----------------------------------------------------------------------*/
void c1_out_gid_sol_str(
                        FILE       *out, /* File pointer to flavia.res */
                        FIELD *actfield, /* active field               */
                        INT       place, /* current solution           */
                        INT         init /* allocate/free memory       */
                        );
#endif
/*! @} (documentation module close)*/
