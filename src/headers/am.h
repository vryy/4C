/*!------------------------------------------------------------------------
\file
\brief am.h header to the AM-System

------------------------------------------------------------------------*/
/*! 
\addtogroup AMSYSTEM 
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief main structure all kinds of fields are kept

m.gee 6/01  

main structure all kinds of fields are kept with         

-------------------------------------------------------------------------*/
typedef struct _ARRAY
{
char                name[9];           /*!< name of the field (just for fun) */
int                 fdim;              /*!< first dimension of field         */
int                 sdim;              /*!< scnd dimension of field          */
enum         
   {
    cca_XX,                            /*!< not defined    */
    cca_DA,                            /*!< double array   */
    cca_DV,                            /*!< double vector  */
    cca_IA,                            /*!< integer array  */
    cca_IV                             /*!< integer vector */
   }                Typ;               /*!< enum type of field */
union
   {
    int     *iv;                       /*!< integer vector */
    double  *dv;                       /*!< double vector  */
    int    **ia;                       /*!< integer array  */
    double **da;                       /*!< double array   */
   }                a;                 /*!< ptr used for calculations        */
#ifdef DEBUG 
struct _TRACEARRAY  *mytracer;         /*!< bugtracing information */
#endif
} ARRAY;

/*!-------------------------------------------------------------------------
\brief main structure all kinds of 3D & 4D fields are kept with   

m.gee 12/01 

main structure all kinds of 3D & 4D fields are kept with 

-------------------------------------------------------------------------*/
typedef struct _ARRAY4D
{
char                name[9];           /*!< name of the field (just for fun) */
int                 fdim;              /*!< first dimension of field         */
int                 sdim;              /*!< scnd dimension of field          */
int                 tdim;              /*!< third dimension of field         */
int                 fodim;             /*!< fourth dimension of field        */
enum         
   {
    cca_XX4D,                              /*!< not defined    */
    cca_D3,                                /*!< double 3D-array   */
    cca_D4,                                /*!< double 4D-array  */
    cca_I3,                                /*!< integer 3D-array  */
    cca_I4                                 /*!< integer 4D-array  */
   }                Typ;               /*!< enum type of field */
union
   {
     double   ***d3;                   /*!< 3D - double array */
     double  ****d4;                   /*!< 4D - double array */
     int      ***i3;                   /*!< 3D - integer array */
     int     ****i4;                   /*!< 4D - integer array */
   }                a;                 /*!< name of union */
#ifdef DEBUG                           
struct _TRACEARRAY  *mytracer;         /*!< bugtracing information */
#endif
} ARRAY4D;
/*! @} (documentation module close)*/
