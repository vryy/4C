/*------------------------------------------------------------------------*
 | structures belonging to the AM-System                    m.gee 6/01    |
 *------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*
 | main structure all kinds of fields are kept with         m.gee 6/01    |
 *------------------------------------------------------------------------*/
typedef struct _ARRAY
{
char                name[9];           /* name of the field (just for fun) */
int                 fdim;              /* first dimension of field         */
int                 sdim;              /* scnd dimension of field          */
enum         
   {
    XX,                                /* not defined    */
    DA,                                /* double array   */
    DV,                                /* double vector  */
    IA,                                /* integer array  */
    IV                                 /* integer vector */
   }                Typ;               /* enum type of field */
union
   {
    int     *iv;                       /* integer vector */
    double  *dv;                       /* double vector  */
    int    **ia;                       /* integer array  */
    double **da;                       /* double array   */
   }                a;                 /* ptr used for calculations        */
#ifdef DEBUG 
int                 place_in_trace;    /* place in bugtracing system       */
#endif
} ARRAY;

/*------------------------------------------------------------------------*
 | main structure all kinds of 3D & 4D fields are kept with   m.gee 12/01 |
 *------------------------------------------------------------------------*/
typedef struct _ARRAY4D
{
char                name[9];           /* name of the field (just for fun) */
int                 fdim;              /* first dimension of field         */
int                 sdim;              /* scnd dimension of field          */
int                 tdim;              /* third dimension of field         */
int                 fodim;             /* fourth dimension of field        */
enum         
   {
    XX4D,                              /* not defined    */
    D3,                                /* double 3D-array   */
    D4,                                /* double 4D-array  */
    I3,                                /* integer 3D-array  */
    I4                                 /* integer 4D-array  */
   }                Typ;               /* enum type of field */
union
   {
     double   ***d3;
     double  ****d4;
     int      ***i3;
     int     ****i4;
   }                a;                 /* ptr used for calculations        */
#ifdef DEBUG 
int                 place_in_trace;    /* place in bugtracing system       */
#endif
} ARRAY4D;
