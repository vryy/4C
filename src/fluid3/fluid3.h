/*!---------------------------------------------------------------------
\file
\brief The 3D fluid element

------------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief fluid2                                                 

<pre>                                                         genk 03/02 

In this structure all variables used for element evaluation by the 3D
fluid element fluid2 are stored.

</pre>

--------------------------------------------------------------------------*/
typedef struct _FLUID3
{
int                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
int                nGP[3];   /*!< number of gaussian points in rs direct. */
int                is_ale;   /*!< flag whether there is ale to me or not  */
struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */

/*--------------------------------------------------- stabilisation flags */
int                istabi;   /*!< stabilasation: 0=no; 1=yes		  */
int                iadvec;   /*!< adevction stab.: 0=no; 1=yes  	  */
int                ipres;    /*!< pressure stab.: 0=no; 1=yes		  */
int                ivisc;    /*!< diffusion stab.: 0=no; 1=GLS-; 2=GLS+   */
int                icont;    /*!< continuity stab.: 0=no; 1=yes 	  */
int                istapa;   /*!< version of stab. parameter		  */
int                istapc;   /*!< flag for stab parameter calculation	  */
int                mk;       /*!< 0=mk fixed 1=min(1/3,2*C); -1 mk=1/3    */
int                ihele[3]; /*!< x/y/z length-def. for vel/pres/cont stab*/
int                ninths;   /*!< number of integ. points for streamlength*/

/*---------------------------------------------------- stabilisation norm */
int                norm_p;   /*!< p-norm: p+1<=infinity; 0=Max.norm       */

/*----------------------------------------------- stabilisation constants */
double             clamb;

/*------------------------------------ statiblisation control information */
int                istrle;   /*!< has streamlength to be computed	  */
int                ivol ;    /*!< calculation of area length		  */
int                iduring;  /*!< calculation during int.-pt.loop	  */
int                itau[3];  /*!< has diagonals etc. to be computed	  */
int                idiaxy;   /*!< flags for tau_? calculation 
			     	  (-1: before; 1: during		  */

/*--------------------------------- element sizes for stability parameter */
double             hk[3];    /*!< vel/pres/cont                           */
} FLUID3;

