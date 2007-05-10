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

  enum _FIELDTYP fieldtyp;     /*!< typ of field */
  INT            handsize;     /*!< has to do with restart */
  long int      *handles;      /*!< has to do with restart */
  DOUBLE        *dvec;         /*!< global redundant vector passed to elements */
  DOUBLE        *dirich;       /*!< ? */
  INT            global_numeq; /*!< size of dvec */
  DOUBLE        *dirichfacs;   /*!< factors for rhs-entries due to prescribed displacements */
  INT            isdyn;        /*!< flag for dynamic or static calculation
                                 isdyn = 0 for static calculation
                                 isdyn = 1 for dynamic calculation */
  DOUBLE         ekin;         /*!< kinetic energy, calculated at element level */
  INT            kintyp;       /*!< kinematic to be used
                                 kintyp = 0: linear kinematic
                                 kintyp = 1: updated lagrange
                                 kintyp = 2: total lagrange */
  INT            kstep;        /*!< time in increment step we are in */
  INT            disnum;       /*!< which discretisation we have */
  INT            inherit;
  INT            point_neum;

  INT            quality;      /*!< element quality measure */

#ifdef D_ALE
  DOUBLE         min, max;     /*<! scaling parameters for ale two_step */
  DOUBLE         min_stiff;
  DOUBLE         max_stiff;
#endif

  INT            pos;          /*<! sol_increment[pos] contains dbc in ale */
  INT            coupl_typ;    /*!< conforming or non-conf. discretization */
  DOUBLE         relax_param;  /*!< the relaxation parameter omega */


#ifdef D_FLUID               /* ab hier fuer fluid */
  DOUBLE        *frhs;
  DOUBLE        *liftdrag;
  DOUBLE        *fidrichrhs;   /*!< for storing the pressure rhs values */
  DOUBLE        *ftimerhs_pro;
  INT            nii;


  INT            error_norm;   /*!< error norm for error calculation
                                 = 0 ... infinity norm
                                 = 1 ... L1 norm
                                 = 2 ... L2 norm */
  DOUBLE         vel_error;    /*!< vel error in the above given norm */
  DOUBLE         pre_error;    /*!< pre error in the above given norm */

  DOUBLE         vel_norm;     /*!< above given norm of the anal. sol */
  DOUBLE         pre_norm;     /*!< above given norm of the anal. sol */


#ifdef D_FLUID_PM
  DOUBLE        *fgradprhs;
#endif

  INT            turbu;

#ifdef D_FLUID2TU
  INT            niturbu_pro;
  INT            niturbu_n;
  DOUBLE        *ftimerhs;     /* used within turbulence (only??) */
#endif

  enum _FLUID_STRESS str;
  INT           *iedgnod;
  INT            ngnode;
  INT            is_relax;      /*!< flag, if calculation is for relaxation parameter */
  /*!< is_relax = 0 -> fluid results are read from sol_increment[3][i] */
  /*!< is_relax = 1 -> fluid results are read from sol_increment[7][i] */

#endif /* ifdef D_FLUID */

#ifdef D_OPTIM                /* include optimization code to ccarat        */
  DOUBLE         getvalue ;     /*!< optimization */
  DOUBLE        *getvector;     /*!< optimization */
#endif                        /* stop including optimization code to ccarat */

#ifdef D_TSI
  INT            disnum_s;   /* structure discretisation index ( ==0 ) */
  INT            disnum_t;   /* thermo-discretisation index ( ==0 ) */
  INT            isoltemdn;  /* position in sol array of presc. temperatures */
#endif
} CONTAINER;
