/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | general fsi variables                                  genk 09/02    |
 *----------------------------------------------------------------------*/
typedef struct _FSI_DYNAMIC
{
FSI_COUPLING       ifsi;            /*!< coupling algorithm */
INT                ipre;            /*!< type of predictor */
INT                inrmfsi;         /*!< convergence criterion */
INT                ichecke;         /*!< energy check */
INT                isdmax;          /*!< Max. No. of steepest descent iterations */
INT                nstep;           /*!< number of steps */
INT                itemax;          /*!< max. number of iterations over fields */
INT                uppss;           /*!<  */
INT                upres;           /*!< update .flavia.res every step */
INT                uprestart;       /*!< write restart every step */
INT                step;            /*!<  */
INT                iale;            /*!<  */
DOUBLE             time;            /*!<  */
DOUBLE             dt;              /*!< time increment */
DOUBLE             maxtime;         /*!< total time */
DOUBLE             entol;           /*!< tolerance for energy check over fields */
DOUBLE             relax;           /*!< actual relaxation parameter */
DOUBLE             convtol;         /*!< tolerance for iteration over fields */
DOUBLE             deltaeint;       /*!< energy production at the interface */
ARRAY              sid;             /*!< structural interface dofs */
INT                numsid;          /*!< number of structural interface dofs */
INT                actpos;          /*!<  */
INT                coupmethod;      /*!< flag, 0=mortar , 1=conforming */
enum {
      cf_none,       /*! No evaluation of coupling force                */
      cf_stress,     /*! Evaluation using stress values (derivatives)   */
      cf_nodeforce   /*! Evaluation via consistent nodal forces         */
     } coupforce;    /*!< how to calculate fsi coupling force */
} FSI_DYNAMIC;


/*----------------------------------------------------------------------*
 | general ale dynamic variables                            ck 12/02    |
 *----------------------------------------------------------------------*/
typedef struct _ALE_DYNAMIC
{
enum
   {
                  classic_lin,   /*!< classic linear calculation */
                  incr_lin,      /*!< classic linear calculation (same as above, but solving for displ. increments) */
                  min_Je_stiff,  /*!< incremental calculation
                                      stiffened with min J_element^2 */
                  two_step,      /*!< calculation in 2 steps per timestep */
                  springs,       /*!< springs rather than continous pseudo material */
                  springs_fixed_ref,    /*!< springs using the initial configuration for stiffness computation and reference */
                  laplace,       /*!< Laplace smoothing algorithm */
                  LAS            /*!< large amplitude sloshing */
   } typ;                        /*!< switch dynamic algorithm */

enum
   {
             no_quality,         /*!< no element quality monitoring */
             aspect_ratio,       /*!< aspect ratio element quality criterion */
             corner_angle,       /*!< corner angle element quality criterion */
	     min_detF            /*!< minimal elemental Jacobian determinant */
   } measure_quality;            /*!< switch, which quality measure to watch */

INT                nstep;        /*!< number of steps */
INT                step;         /*!< actual step */
INT                updevry_disp; /*!< write result very updevry step */
INT                num_initstep; /*!< number of initial steps with prestress */
DOUBLE             dt;           /*!< stepsize */
DOUBLE             maxtime;      /*!< maximum total time */
DOUBLE             time;         /*!< actual time */
INT                coupmethod;   /*!< flag, 0=mortar , 1=conforming, 2=nonconforming */
} ALE_DYNAMIC;


