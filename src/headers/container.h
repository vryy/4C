/*!----------------------------------------------------------------------
\brief file pointers

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

<pre>                                                         m.gee 10/02
This Container is used to transport variables from the 'steuerroutinen' 
to the 'Elementroutinen'
</pre>

*----------------------------------------------------------------------*/
typedef struct _CONTAINER
{
/*------------------------------------------------------------- file I/O */
enum _FIELDTYP fieldtyp;     /*!< typ of field */
INT            handsize;     /*!< has to do with restart */
long int      *handles;      /*!< has to do with restart */                         
DOUBLE        *dvec;         /*!< global redundant vector passed to elements */
DOUBLE        *dirich;       /*!< ? */
INT            global_numeq; /*!< size of dvec */
DOUBLE        *dirichfacs;   /*!< factors for rhs-entries due to prescribed displacements */
INT            isdyn;        /*!< flag for dynamic or static calculation */
                             /*!< isdyn = 0 for static calculation */
                             /*!< isdyn = 1 for dynamic calculation */
DOUBLE         ekin;         /*!< kinetic energy, calculated at element level */
INT            kintyp;       /*!< kinematic to be used */
                             /*!< kintyp = 0: linear kinematic */
                             /*!< kintyp = 1: updated lagrange */
                             /*!< kintyp = 2: total lagrange */
INT            kstep;        /*!< time in increment step we are in */
INT            actndis;      /*!< which discretisation we have */
INT            inherit;
INT            point_neum;

INT            quality;      /*!< element quality measure */

DOUBLE         min, max;     /*<! scaling parameters for ale two_step */
DOUBLE         min_stiff;
DOUBLE         max_stiff; 
                                                   
INT            pos;          /*<! sol_increment[pos] contains dbc in ale */

#ifdef D_FLUID               /*!< ab hier fuer fluid */ 
DOUBLE        *ftimerhs;       
DOUBLE        *fiterhs;
DOUBLE        *liftdrag;
DOUBLE        *fidrichrhs;   /*!< for storing the pressure rhs values */
DOUBLE        *ftimerhs_pro;
INT            nii;
INT            nif;
INT            nim;
struct _DBCSR *gradmatrix;   /*!< gradient matrix Projection Method */
struct _DBCSR *lumpedmass;   /*!< gradient matrix Projection Method */
INT            turbu;
INT            niturbu_pro;
INT            niturbu_n;
enum _FLUID_STRESS str;         
INT            is_relax;      /*!< flag, if calculation is for relaxation parameter */
                              /*!< is_relax = 0 -> fluid results are read from sol_increment[3][i] */
			      /*!< is_relax = 1 -> fluid results are read from sol_increment[7][i] */
INT		gen_alpha;    /*!< generalised alpha time integration algorithm flag */
#endif

#ifdef D_OPTIM                /* include optimization code to ccarat        */
DOUBLE         getvalue ;     /*!< optimization */   
DOUBLE        *getvector;     /*!< optimization */   
#endif                        /* stop including optimization code to ccarat */
} CONTAINER;
