/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         lch 10/02
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
int            kstep;        /*!< time in increment step we are in */
} CONTAINER;
