/*!----------------------------------------------------------------------
\file
\brief headerfile for 3D beam element, contains structures and prototypes

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief beam3 data set for numerical integration

<pre>                                                              fh 09/02
This structure contains all data set for numerical integration
</pre>

*----------------------------------------------------------------------*/
typedef struct _B3_DATA
{
DOUBLE	        xgrr[13]; /*!< gauss point coordinate r */
DOUBLE          wgtr[13]; /*!< weight of gauss point coordinate r */
DOUBLE          xlst[13]; /*!< lobatto point coordinate s,t */
DOUBLE          wlst[13]; /*!< weight of lobatto point coordinate s,t */
DOUBLE          xr[3];	  /*!< local r-coordinate of element nodes */
} B3_DATA;


/*!----------------------------------------------------------------------
\brief beam3 working array for element input data set

<pre>                                                              fh 09/02
This structure contains all input data set (number of GP, Hinge Code,
material, cross section, geometry) of a specific beam element
</pre>

*----------------------------------------------------------------------*/
typedef struct _BEAM3
{

INT           nGP[1];  /*!< number of gauss points for the element */
INT           hc[25];  /*!< hinge code for the element */

DOUBLE        nref[3]; /*!< reference node to determine local z-axis */
DOUBLE        area;    /*!< area of element cross section */
DOUBLE        gs;      /*!< 1/shear correction factor of element cross section */
DOUBLE        width;   /*!< width of element cross section */
DOUBLE        height;  /*!< height of element cross section */
DOUBLE        iyy;     /*!< Iyy of element cross section */
DOUBLE        izz;     /*!< Izz of element cross section */
DOUBLE        iyz;     /*!< Iyz of element cross section */
DOUBLE        it;      /*!< It of element cross section */
DOUBLE        iuu;     /*!< first principle moment of inertia */
DOUBLE        ivv;     /*!< second principle moment of inertia */
INT           ike;     /*!< linear elastic stiffness matrix (1) or not (2) */
DOUBLE        alpha;   /*!< angle of principle axes of cross section */
DOUBLE        length;  /*!< length of the element */

struct _ARRAY	elewa; /*!< element working array for plasticity*/

/*----------------------------------------------------------------------*
 | result struct: stresses                                  fh 09/02    |
 *----------------------------------------------------------------------*/
struct _ARRAY4D  force_GP;
struct _ARRAY4D  force_ND;
} BEAM3;


/*!----------------------------------------------------------------------
\brief controls beam3 element

<pre>                                                              fh 09/02
This routine is the main beam3 control routine

</pre>
\param *actfield        FIELD        (i)  actual field
\param *actpart         PARTITION    (i)  actual partition
\param *actintra        INTRA        (i)  actual intra
\param *actele          ELEMENT     (i/o) actual element
\param *estif_global    ARRAY       (i/o) global element stiffness vector
\param *emass_global    ARRAY       (i/o) global element mass vector
\param *intforce_global ARRAY       (i/o) global internal force vector
\param *action          CALC_ACTION (i/o)  actual action
\param *container       CONTAINER   (i/o) container

\warning There is nothing special in this routine
\return void
\sa calling:   b3_init() , b3_static_ke() , b3_bop, b3_trans_stf() ,
               b3_load() , b3_loadlin() , b3_cal_stress() ,
	       b3_cal_stresslin() , b3_setdirich()
    called by: calinit() , calelm()

*----------------------------------------------------------------------*/
void beam3(FIELD       *actfield,
           PARTITION   *actpart,
           INTRA       *actintra,
           ELEMENT     *ele,
           ARRAY       *estif_global,
           ARRAY       *emass_global,
           ARRAY       *intforce_global,
           CALC_ACTION *action,
	   CONTAINER   *container);

/*!----------------------------------------------------------------------
\brief reads all datas of the beam element from input file

<pre>                                                              fh 09/02
This routine reads all datas of the actual beam element from input file

</pre>
\param *ele      ELEMENT  (i/o)  actual element


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/
void b3inp(ELEMENT *ele);
#endif
/*! @} (documentation module close)*/
