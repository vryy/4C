#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
/*----------------------------------------------------------------------*
 | integration of element loads                              ah 07/02   |
 | in here, line and surface loads are integrated                       |
 *----------------------------------------------------------------------*/
void w1_eleload(ELEMENT  *ele,
                W1_DATA  *data,
                MATERIAL *mat,
                double	 *loadvec,
                int	  init)
{
int          lr,ls;
int          i,j,k;
int          inode,idof;
int          nir,nis;
int          iel;
int          nd;
const int    numdf  = 2;   /* number dof per node   */
int          foundsurface;

double       thickness;
double       e1,e2;
double       facr,facs,wgt;
double       det;
double       xi,yi;

static ARRAY eload_a; static double **eload;
static ARRAY funct_a; static double  *funct;
static ARRAY deriv_a; static double **deriv;
static ARRAY xjm_a;   static double **xjm;

W1_DATA      actdata;
MATERIAL    *actmat;
NODE        *actnode;
/*--------------------- variables needed for integration of line loads */
int             ngline;
int             ngnode;
int             gnode[3];
int             line;
int             foundline;
GLINE          *gline[4];
NEUM_CONDITION *lineneum[4];
int             ngp;
int             gp;
double          xgp[3];
double          xgp_n[3];
double          wgp[3];
int             dir;
double          ds;
double          ap[3],ar[3];

/*----------------------------------------------------------------------*/
/* init phase        (init=1)                                           */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("w1_eleload");
#endif

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_eleload */





/*----------------------------------------------------------------------*
 | integration of element loads                             ah 07/02    |
 *----------------------------------------------------------------------*/
void w1_fextsurf(ELEMENT    *ele,
                double    **eload,
                double     *funct,
                int         fac,
                int         iel)
{
double    force[2];
int i,j;/*---------- ------------i = loaddirection x or y---------------*/
#ifdef DEBUG 
dstrc_enter("w1_fextsurf");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
switch(ele->g.gsurf->neum->neum_type)
{
case pres_domain_load:
/*----------------------------------------------------------------------*/
/*                                      uniform prescribed surface load */
/*----------------------------------------------------------------------*/
   for (i=0; i<1; i++)
   {
    force[i] = ele->g.gsurf->neum->neum_val.a.dv[i];
   }
/*-------------------- add load vector component to element load vector */
   for (j=0; j<iel; j++)
   {
      for (i=0; i<1; i++)
      {
         eload[i][j] += funct[j] * force[i] * fac;
      }
   }
break;
/*----------------------------------------------------------------------*/
case neum_live:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
case neum_consthydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
case neum_increhydro_z:
dserror("load case unknown");
break;
/*----------------------------------------------------------------------*/
default:
break;

/*----------------------------------------------------------------------*/
}/* end of switch(ele->g.gsurf->neum->neum_type)*/
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of  w1_fextsurf*/
