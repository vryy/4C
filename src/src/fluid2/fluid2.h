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

#ifndef FLUID2_H
#define FLUID2_H

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
INT                create_ale;

struct _ELEMENT   *ale_ele;  /*!< pointer to the corresponding ale element */


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

/*----- flag, if there is a lif&drag or fsi coupling line to this element */
INT                force_on;
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

#ifdef QUASI_NEWTON
  struct _ARRAY   estif;
#endif

#ifdef D_FLUID2_TDS
  struct _ARRAY    sub_pres;
    /*
     * USAGE:
     *    
     *            sub_pres.a.dv[number of gausspoint]
     *
     *            whereas
     *
     *                                 number of gausspoint in r-direction
     *                                 |
     *                                 |        number of gausspoint
     *                                 |        in s-direction
     *                                 |        |
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */
    
  struct _ARRAY    sub_pres_acc;
    /* only required for incremental acceleration genalpha
     * 
     * USAGE:
     *    
     *            sub_pres.a.dv[number of gausspoint]
     *
     *            whereas
     *
     *                                 number of gausspoint in r-direction
     *                                 |
     *                                 |        number of gausspoint
     *                                 |        in s-direction
     *                                 |        |
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */

  struct _ARRAY    sub_pres_trial;
    /*
     * USAGE:
     *    
     *            sub_pres.a.dv[number of gausspoint]
     *
     *            whereas
     *
     *                                 number of gausspoint in r-direction
     *                                 |
     *                                 |        number of gausspoint
     *                                 |        in s-direction
     *                                 |        |
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */
    
  struct _ARRAY    sub_pres_acc_trial;
    /* only required for incremental acceleration genalpha
     * 
     * USAGE:
     *    
     *            sub_pres.a.dv[number of gausspoint]
     *
     *            whereas
     *
     *                                 number of gausspoint in r-direction
     *                                 |
     *                                 |        number of gausspoint
     *                                 |        in s-direction
     *                                 |        |
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */

    
  struct _ARRAY  sub_vel;
    /*
     * USAGE:
     *    
     *            sub_vel.a.da[dim][number of gausspoint]
     *
     *            whereas
     *
     *                        xdim=0, ydim=1
     *
     *            and
     *
     *            
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */

  struct _ARRAY  sub_vel_acc;
    /* only required for incremental acceleration genalpha
     * 
     * USAGE:
     *    
     *            sub_vel.a.da[dim][number of gausspoint]
     *
     *            whereas
     *
     *                        xdim=0, ydim=1
     *
     *            and
     *
     *            
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */
  
     struct _ARRAY  sub_vel_trial;
    /*
     * USAGE:
     *    
     *            sub_vel.a.da[dim][number of gausspoint]
     *
     *            whereas
     *
     *                        xdim=0, ydim=1
     *
     *            and
     *
     *            
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */

  struct _ARRAY  sub_vel_acc_trial;
    /* only required for incremental acceleration genalpha
     * 
     * USAGE:
     *    
     *            sub_vel.a.da[dim][number of gausspoint]
     *
     *            whereas
     *
     *                        xdim=0, ydim=1
     *
     *            and
     *
     *            
     *            number of gausspoint=lr*(nis)+ls
     *                                      |
     *                                     total number of gausspoint
     *                                     in s-direction
     * */ 
#endif

    
} FLUID2;


#endif
#endif
/*! @} (documentation module close)*/
