/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
/*!---------------------------------------------------------------------
\brief fluid2

<pre>                                                         genk 03/02

In this structure all variables used for element evaluation by the 2D
fluid element fluid2 are stored.
The stabilisation flags and input parameters have been stored within the
union 'stabi' which is a condition assigned to the volume or surface which
contains the fluid elements. (chfoe 01/04)

</pre>

--------------------------------------------------------------------------*/
#ifdef D_FLUID2
typedef struct _FLUID2
{
INT                nGP[2];   /*!< number of gaussian points in rs direct. */
INT                is_ale;   /*!< flag whether there is ale to me or not  */

/*---------------------------------------------------- stabilisation ---*/
enum _STABILISATION_TYP	stab_type;	/* enum of stabilisation	*/
union
{
  struct _STAB_PAR_GLS  *gls;	/*! pointer to stabilisation parameters	*/
/*  struct _STAB_PRES_PRO *pp; */
} stabi;

/*--------------------------------- element sizes for stability parameter */
DOUBLE             hk[3];    /*!< vel/pres/cont                           */
struct _ARRAY      tau_old;
/*-------------------------- flag for turbulence  1=algebraic, 2=ke-model, 3=ko-model */
INT                turbu;

/*------------------------------------------------ free surface parameter */
INT                fs_on;   /*! element belongs to free surface           */

/*------------------------------------------------- structure for submesh */
#ifdef FLUID2_ML
INT                smisal;        /* flag for element submesh creation    */
DOUBLE             smcml;         /* charact. mesh length for submesh     */
struct _ARRAY      xyzsm;         /* coordinates of submesh nodes         */
struct _ARRAY      solsm;         /* sol. of current timestep             */
struct _ARRAY      solsmn;        /* sol. of last timestep                */

/*--------------------------------------------- structure for sub-submesh */
struct _ARRAY      xyzssm;        /* coordinates of sub-submeshnodes      */
#endif
/*-------------------------------------------------------- stress results */
struct _ARRAY      stress_ND;    /*!< nodal stresses                         */
/*------------------------------------------------------- nodal curvature */
struct _ARRAY      kappa_ND;
} FLUID2;

#endif
/*! @} (documentation module close)*/
