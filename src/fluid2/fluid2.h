/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
/*!---------------------------------------------------------------------
\brief fluid2                                                 

<pre>                                                         genk 03/02 

In this structure all variables used for element evaluation by the 2D
fluid element fluid2 are stored.

</pre>

--------------------------------------------------------------------------*/
#ifdef D_FLUID2
typedef struct _FLUID2
{
INT                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
INT                nGP[2];   /*!< number of gaussian points in rs direct. */
INT                is_ale;   /*!< flag whether there is ale to me or not  */
struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */

/*--------------------------------------------------- stabilisation flags */
INT                istabi;   /*!< stabilasation: 0=no; 1=yes              */
INT                iadvec;   /*!< adevction stab.: 0=no; 1=yes            */
INT                ipres;    /*!< pressure stab.: 0=no; 1=yes             */
INT                ivisc;    /*!< diffusion stab.: 0=no; 1=GLS-; 2=GLS+   */
INT                icont;    /*!< continuity stab.: 0=no; 1=yes           */
INT                istapa;   /*!< version of stab. parameter              */
INT                istapc;   /*!< flag for stab parameter calculation     */
INT                mk;       /*!< 0=mk fixed 1=min(1/3,2*C); -1 mk=1/3    */ 
INT                ihele[3]; /*!< x/y/z length-def. for vel/pres/cont stab*/
INT                ninths;   /*!< number of integ. points for streamlength*/

/*---------------------------------------------------- stabilisation norm */
INT                norm_p;   /*!< p-norm: p+1<=infinity; 0=Max.norm       */

/*----------------------------------------------- stabilisation constants */
DOUBLE             clamb;

/*------------------------------------ statiblisation control information */
INT                istrle;   /*!< has streamlength to be computed         */
INT                iarea;    /*!< calculation of area length              */
INT                iduring;  /*!< calculation during INT.-pt.loop         */
INT                idiaxy;   /*!< has diagonals etc. to be computed       */
INT                itau[3];  /*!< flags for tau_? calculation 
                                  (-1: before; 1: during                  */

/*--------------------------------- element sizes for stability parameter */
DOUBLE             hk[3];    /*!< vel/pres/cont                           */

/*-------------------------- flag for turbulence  1=algebraic, 2=ke-model, 3=ko-model */
INT                turbu;

/*------------------------------------------------ free surface parameter */
INT                fs_on;   /*! element belongs to free surface           */

/*------------------------------------------------- structure for submesh */
INT                smisal;        /* flag for element submesh creation    */
DOUBLE             smcml;         /* charact. mesh length for submesh     */
struct _ARRAY      xyzsm;         /* coordinates of submesh nodes         */
struct _ARRAY      solsm;         /* sol. of current timestep             */
struct _ARRAY      solsmn;        /* sol. of last timestep                */

/*--------------------------------------------- structure for sub-submesh */
struct _ARRAY      xyzssm;        /* coordinates of sub-submeshnodes      */

/*-------------------------------------------------------- stress results */
struct _ARRAY      stress_ND; /*!< nodal stresses                         */
/*------------------------------------------------------- nodal curvature */
struct _ARRAY      kappa_ND;
} FLUID2;

#endif
/*! @} (documentation module close)*/
