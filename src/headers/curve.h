/*----------------------------------------------------------------------*
 | curve types                                            m.gee 2/02    |
 | there can be different types of load curves                          |
 *----------------------------------------------------------------------*/
typedef enum _CURVTYP
{
                       curve_none,
                       curve_polygonal,  /* polygonal, piecewise linear curve */    
                       curve_explicit    /* excplicit function */
} CURVTYP;
/*----------------------------------------------------------------------*
 | general dynamic-curves                                 m.gee 4/01    |
 | the arrays time and value are dependent on type of load curve and    |
 | size of the curve                                                    |
 *----------------------------------------------------------------------*/
typedef struct _CURVE
{
INT                    Id;            /* Id of the load curve */
enum _CURVTYP          curvetyp;      /* type of load curve */
INT                    bystep;        /* flag whether curve operated by number of steps or in absolut time */
INT                    numex;         /* number of explicit function */
ARRAY                  time;          /* array for time steps */
ARRAY                  value;         /* array for values at time steps */
DOUBLE                 T;             /* forgot about it..... */
DOUBLE                 c1;            /* constant for explicit functions */
DOUBLE                 c2;            /* constant for explicit functions */
} CURVE;
