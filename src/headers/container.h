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
int            handsize;     /*!< has to do with restart */
long int      *handles;      /*!< has to do with restart */                         
double        *dvec;         /*!< global redundant vector passed to elements */
double        *dirich;       /*!< ? */
int            global_numeq; /*!< size of dvec */
double        *dirichfacs;   /*!< factors for rhs-entries due to prescribed displacements */
int            isdyn;        /*!< flag for dynamic or static calculation */
                            /*!< isdyn = 0 for static calculation */
                            /*!< isdyn = 1 for dynamic calculation */
double         ekin;         /*!< kinetic energy, calculated at element level */
int            kintyp;      /*!< kinematic to be used */
                             /*!< kintyp = 0: linear kinematic */
                             /*!< kintyp = 1: updated lagrange */
                             /*!< kintyp = 2: total lagrange */
int            kstep;        /*!< time in increment step we are in */
int            actndis;      /*!< which discretisation we have */
int            inherit;
int            point_neum;

#ifdef D_FLUID
double        *ftimerhs;     /*!< ab hier fuer fluid */   
double        *fiterhs;
double        *ftimerhs_pro;
int            nii;
int            nif;
int            turbu;
int            niturbu_pro;
int            niturbu_n;
enum _FLUID_STRESS str;         
#endif

#ifdef D_OPTIM                /* include optimization code to ccarat        */
double         getvalue ;     /*!< optimization */   
double        *getvector;     /*!< optimization */   
#endif                        /* stop including optimization code to ccarat */
} CONTAINER;
