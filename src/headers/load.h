/*----------------------------------------------------------------------*
 | general dynamic-curves                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _CURVE
{
int                    Id;            /* Id of the load curve */
char                   curvetyp[50];  /* type of load curve, will be made an enum some day */
int                    bystep;        /* flag whether curve operated by number of steps or in absolut time */
ARRAY                  time;          /* array for time steps */
ARRAY                  value;         /* array for values at time steps */
double                 T;             /* forgot about it..... */
} CURVE;
