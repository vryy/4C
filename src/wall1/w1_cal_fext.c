#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
typedef enum _RSF
{ /* r..r-dir., s..s-dir., n..-1.0, p..+1.0 */
  rp, rn, ps, ns 
} RSF;

/*----------------------------------------------------------------------*
/*----------------------------------------------------------------------*
 | integration of element loads                              ah 07/02   |
 | in here, line and surface loads are integrated                       |
 *----------------------------------------------------------------------*/
 /*----------------------------------------------------------------------*
 |                 s                            s                        |
 |                 |                            |                        |
 |         1-------4--------0          1------line 1-----0               |
 |         |       |        |          |        |        |               |
 |         |       |        |          |        |        |               |
 |         |quad8/9|        |          l quad4  |        l               |
 |         |       |        |          i        |        i               |
 |      r--5------(8)-------7--r    r--n--------|--------n--r            |
 |         |       |        |          e        |        e               |
 |         |       |        |          2        |        4               |
 |         |       |        |          |        |        |               |
 |         |       |        |          |        |        |               |
 |         2-------6--------3          2------line 3-----3               |
 |                 |                            |                        |
 |                 s                            s                        |
 *----------------------------------------------------------------------*/

void w1_eleload(ELEMENT  *ele,     /* actual element                  */
                W1_DATA  *data,    /* wall1- Data                     */
                double	 *loadvec, /* global element load vector fext */
                int	  init)    /* flag if init or calculation     */
{
int          lr,ls;              /* integration directions          */
int          i,j,k;              /* some loopers                    */
int          inode,idof;         /* some loopers                    */
int          nir,nis;            /* number of GP's in r-s direction */
int          iel;                /* number of element nodes         */
int          nd;                 /* element DOF                     */
const int    numdf  = 2;         /* dof per node                    */
int          foundsurface;       /* flag for surfaceload present    */

/* double       thickness;            */
double       e1,e2;              /* GP-koordinates in r-s-system   */
double       fac,facr,facs,wgt;  /* integration factor  GP-info    */
double       det;                /* det of jacobian matrix         */

static ARRAY eload_a; static double **eload;  /* static element load vector */
static ARRAY funct_a; static double  *funct;  /* ansatz-functions           */
static ARRAY deriv_a; static double **deriv;  /* derivatives of ansatz-funct*/
static ARRAY xjm_a;   static double **xjm;    /* jacobian matrix            */

/*--------------------- variables needed for integration of line loads */
int             foundline;   /* flag for lineload present or not       */
int             ngline;      /* number of geometrylines to the element */
GLINE          *gline[4];    /* geometrylines of the element           */
NEUM_CONDITION *lineneum[4]; /* short grep on line-neum. conditions    */
int             line;        /* looper over lines                      */
int             ngnode;      /* number of geometry-nodes on g-line     */
int             ngr,ngs;     /* number of GP'e for line-integration    */
double          xgp[3];      /* Coordinates of GP'es                   */
double          ygp[3];      /* Coordinates of GP'es                   */
double          wgx[3];      /* weights at GP'es                       */
double          wgy[3];      /* weights at GP'es                       */
double          ds;       /* dx/dr line increment for line integration */
double          facline;     /*integration factor for line integration */
double          forceline[2];/* lineload value in x and y direction(inp)*/
RSF rsgeo;                   /* integration direction on line          */

/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("w1_eleload");
#endif

if (init==1)
{
  eload   = amdef("eload"  ,&eload_a,MAXDOFPERNODE,MAXNOD_WALL1,"DA");
  funct   = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1            ,"DV");       
  deriv   = amdef("deriv"  ,&deriv_a,2            ,MAXNOD_WALL1,"DA");       
  xjm     = amdef("xjm_a"  ,&xjm_a  ,numdf        ,numdf       ,"DA");
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/

else if (init==-1)
{
amdel(&eload_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*-------------------------------------------- init the gaussian points */
w1intg(ele,data,1);

nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
iel     = ele->numnp;
nd      = iel*numdf;

/*--------------------------------- check for presence of surface loads */
foundsurface=0;
if (!(ele->g.gsurf->neum)) goto endsurface;
if (ele->g.gsurf->neum->neum_type!=pres_domain_load)
{
  dserror("load case not implemented");
  goto endsurface;
}
foundsurface=1;

/*thickness  = ele->e.w1->thick;*/
/*------------------------ integration loop over GP's of actual surface */

for (lr=0; lr<nir; lr++)/*---------------------------- loop r-direction */
{
   /*=============================== gaussian point and weight at it ===*/
   e1   = data->xgrr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)/*------------------------- loop s direction */
   {
      /*============================ gaussian point and weight at it ===*/
      e2   = data->xgss[ls];
      facs = data->wgts[ls];
      /*----------- shape functions (and (not needed)their derivatives) */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*----------------------------------------------- jacobian matrix */       
      w1_jaco (funct,deriv,xjm,&det,ele,iel); 
      /*------------------------------------------- integration factor  */ 
      fac = facr * facs * det; 
      /*------------------------------------------------- surface-load  */
      w1_fextsurf(ele,eload,funct,fac,iel); 
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */
/*                      the surface load of this element has been done, */
/*                             -> switch off the surface load condition */
ele->g.gsurf->neum=NULL;
/*----------------------------------------------------------------------*/
endsurface:
/*----------------------------------------------------------------------*/
/*--------- integration of line loads on lines adjacent to this element */
/*------------------------------------ check for presence of line loads */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   lineneum[i] = gline[i]->neum;
   if (lineneum[i]) foundline=1;
}
if (foundline==0) goto endline;
/*-------------------------- loop over lines (with neumann conditions) */
for (line=0; line<ngline; line++)
{

   if (lineneum[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*-------------------- coordinates and weights of integration points */
   /*--------- original GP-coordinates and weights for area-integration */
   ngr = nir; 
   ngs = nis; 
    for (i=0; i<ngr; i++) { xgp[i] = data->xgrr[i];
                            wgx[i] = data->wgtr[i]; }
    for (i=0; i<ngs; i++) { ygp[i] = data->xgss[i];
                            wgy[i] = data->wgts[i]; }
   /*--------- degeneration to line-integration-info for r-s directions */
    switch (line)
    {
    case 0:  /* line1 (s=const=+1) - first,last,(quadratic:middle) node */                     
       ngs    =  1; /*s=const -> only one integration point in s-direct */ 
       ygp[0] =  1.;/* line 1 -> s=const=+1                             */
       wgy[0] =  1.;
       rsgeo = rp; /*  r = {+ -> 0 -> -} |s = +1  */
    break;
    case 1:
       ngr    =  1; 
       xgp[0] = -1.;
       wgx[0] =  1.;
       rsgeo = ns; /*  r = -1 | s={+ -> 0 -> -}  */
    break;
    case 2:
       ngs    = 1; 
       ygp[0] = -1.;
       wgy[0] =  1.;
       rsgeo = rn; /*  r = {- -> 0 -> +} |s = -1  */
    break;
    case 3:
       ngr    =  1; 
       xgp[0] =  1.;
       wgx[0] =  1.;
       rsgeo = ps; /*  r =  1 | s={- -> 0 -> +}  */
    break;
    }
    /*----------------------------------- integration loop on actual line */

    for (lr=0; lr<ngr; lr++)/*-------------------------- loop r-direction */
    {
       /*============================= gaussian point and weight at it ===*/
       e1   = xgp[lr];
       facr = wgx[lr];
       for (ls=0; ls<ngs; ls++)/*----------------------- loop s direction */
       {
          /*========================== gaussian point and weight at it ===*/
          e2   = ygp[ls];
          facs = wgy[ls];
          /*--------- shape functions (and (not needed)their derivatives) */
          w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
          /*--------------------------------------------- jacobian matrix */       
          w1_jaco (funct,deriv,xjm,&det,ele,iel); 
          /*---------------------------------------------- line increment */
          ds = 0.0;
          switch (rsgeo)
          {
          case rp: case rn: 
            ds = DSQR(xjm[0][0])+DSQR(xjm[0][1]);
            ds = sqrt(ds);
          break;
          case ps: case ns: 
            ds = DSQR(xjm[1][0])+DSQR(xjm[1][1]);
            ds = sqrt(ds);
          break;
          }
          /*----------------------------------------- integration factor  */
          facline = ds * facr * facs ;
          /*-------------------------------------------------------------*/
          /*                                uniform prescribed line load */
          /*-------------------------------------------------------------*/
             for (i=0; i<numdf; i++)
             {
                forceline[i] = gline[line]->neum->neum_val.a.dv[i];
             }
             /*-------- add load vector component to element load vector */
             for (j=0; j<iel; j++)
             {
                 for (i=0; i<numdf; i++)
                 {
                    eload[i][j] += funct[j] * forceline[i] * facline;
                 }
             }
       }/*========================================== end of loop over ls */ 
    }/*============================================= end of loop over lr */
    /* line number lie has been done,switch of the neumann pointer of it */
   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*-----------------------------------------------------------------------*/
endline:
/*-----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/* add static array eload to global external element load vector        */
/*----------------------------------------------------------------------*/
if (foundsurface+foundline != 0)
for (inode=0; inode<iel; inode++)
{
   for (idof=0; idof<numdf; idof++)
   {
         loadvec[inode*numdf+idof] += eload[idof][inode];
   }
}
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
                double      fac,
                int         iel)
{
double    force[2];
int i,j;        /*-loopers(i = loaddirection x or y)(j=node)------------*/
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
   for (i=0; i<2; i++)
   {
    force[i] = ele->g.gsurf->neum->neum_val.a.dv[i];
   }
/*-------------------- add load vector component to element load vector */
   for (j=0; j<iel; j++)
   {
      for (i=0; i<2; i++)
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
