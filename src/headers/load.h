/*----------------------------------------------------------------------*
 | general dynamic-curves                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _CURVE
{
int                    Id;
char                   curvetyp[50];
int                    bystep;
ARRAY                  time;
ARRAY                  value;
double                 T;
} CURVE;
