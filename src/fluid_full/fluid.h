/*----------------------------------------------------------------------*
 | fluid2                                                  genk 3/02    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID2
{
int                ntyp;     /* flag for element type: 1=rectangle; 2=triangle */
int                nGP[2];   /* number of gaussian points in rs directions */
int                is_ale;   /* flag whether there is ale to me or not */
struct _ELEMENT   *my_ale;   /* pointer to my ale element, otherwise NULL */
/*------------------------------------------------- stabilisation flags */
int                istabi;   /* stabilasation: 0=no; 1=yes */
int                iadvec;   /* adevction stab.: 0=no; 1=yes */
int                ipres;    /* pressure stab.: 0=no; 1=yes */
int                ivisc;    /* diffusion stab.: 0=no; 1=GLS-; 2=GLS+ */
int                icont;    /* continuity stab.: 0=no; 1=yes */
int                istapa;   /* version of stab. parameter */
int                istapc;   /* flag for stab parameter calculation */
int                mk;       /* 0=mk fixed 1=min(1/3,2*C); -1 mk=1/3 */ 
int                ihele[3]; /* x/y/z length-def. for vel/pres/cont stab */
int                ninths;   /* number of integration points for streamlength */
/*-------------------------------------------------- stabilisation norm */
int                norm_p;   /* p-norm: p+1<=infinity; 0=Max.norm */
/*--------------------------------------------- stabilisation constants */
double             clamb;
double             c1;
double             c2;
double             c3;
/*---------------------------------- statiblisation control information */
int                istrle;      /* has streamlength to be computed */
int                iarea;       /* calculation of area length */
int                iduring;     /* calculation during int.-pt.loop */
int                itau[3];     /* flags for tau_? calculation (-1: before; 1: during */
int                idiaxy;      /* has diagonals etc. to be computed */
/*---------------------------- element sizes for stability parameter */
double             hk[3];      /* vel/pres/cont */
} FLUID2;

/*----------------------------------------------------------------------*
 | fluid2 data                                              genk 03/02  |
 *----------------------------------------------------------------------*/
typedef struct _F2_DATA
{
double        qxg[MAXQINTP][MAXQINTC];    /* coordinates for QUADS */
double        qwgt[MAXQINTP][MAXQINTC];   /* weights for QUADS */

double        txgr[MAXTINTP][MAXTINTC];   /* coordinates in r for TRIS */
double        txgs[MAXTINTP][MAXTINTC];   /* coordinates in s for TRIS */
double        twgt[MAXTINTP][MAXTINTC];   /* weights for TRIS */
} F2_DATA;

typedef union _FLUID_DATA
{
   struct _F2_DATA  f2data;
} FLUID_DATA;
