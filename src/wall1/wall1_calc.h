/*----------------------------------------------------------------------*
 | wall1 data                                                al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_DATA
{
double        xgrr[13];
double        wgtr[13];

double        xgss[13];
double        wgts[13];
} W1_DATA;
/*----------------------------------------------------------------------*
 | working array                                             al 6/01    |
 *----------------------------------------------------------------------*/
typedef struct _W1_IP_WA
{
double      sig[4];
double      eps[4];
double       qn[4];

double       epstn;
int          yip;
} W1_IP_WA;
typedef struct _W1_ELE_WA
{
W1_IP_WA      *ipwa;
} W1_ELE_WA;

