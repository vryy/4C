/*!----------------------------------------------------------------------
\brief file pointers

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

DOUBLE         min, max, min_stiff, max_stiff; /*! scaling parameters for 
                                                   ale two_step */
INT            pos;          /*! sol_increment[pos] contains dbc in ale */

#ifdef D_FLUID
DOUBLE        *ftimerhs;     /*!< ab hier fuer fluid */   
DOUBLE        *fiterhs;
DOUBLE        *ftimerhs_pro;
INT            nii;
INT            nif;
INT            turbu;
INT            niturbu_pro;
INT            niturbu_n;
enum _FLUID_STRESS str;         
#endif

#ifdef D_OPTIM                /* include optimization code to ccarat        */
DOUBLE         getvalue ;     /*!< optimization */   
DOUBLE        *getvector;     /*!< optimization */   
#endif                        /* stop including optimization code to ccarat */
} CONTAINER;
