/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_eleload' which calculates element load
       vector due to external element (line and surface) loads for gradient
       enhanced wall element

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"
/*----------------------------------------------------------------------*/
typedef enum _RSF
{ /* r..r-dir., s..s-dir., n..-1.0, p..+1.0 */
  rp, rn, ps, ns 
} RSF;

static ARRAY eload_a;  static DOUBLE **eload;  /* static element load vector */
static ARRAY functd_a; static DOUBLE  *functd; /* ansatz-func. for displacements*/
static ARRAY derivd_a; static DOUBLE **derivd;  /* derivatives of ansatz-funct*/
static ARRAY xjm_a;    static DOUBLE **xjm;    /* jacobian matrix            */

/*! 
\addtogroup WALLGE 
*//*! @{ (documentation module open)*/

#ifdef D_WALLGE

/*!----------------------------------------------------------------------
\brief  calculates element load
       vector due to external element (line and surface) loads for gradient
       enhanced wall element

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
void wge_eleload(ELEMENT     *ele,       /* actual element              */
                 WALLGE_DATA *data,      /* element integration data    */
                 DOUBLE      *loadvec,   /* external element forces     */
                 INT          init) 
{

INT          i,j;                /* element Diplacement-DOF         */
INT          lr,ls;              /* integration directions          */
INT          nodei,nodestarti;              /* integration directions          */
const INT    numdfd  = 2;        /* displacement-dof per node       */
INT          nd;                 /* element Diplacement-DOF         */
INT          nir,nis;            /* number of GP's in r-s direction */
INT          ield;               /* number of element nodes         */
INT          foundsurface;       /* flag for surfaceload present    */

DOUBLE       e1,e2;              /* GP-koordinates in r-s-system   */
DOUBLE       fac,facr,facs,wgt;  /* integration factor  GP-info    */
DOUBLE       det;                /* det of jacobian matrix         */

/*--------------------- variables needed for integration of line loads */
INT             foundline;   /* flag for lineload present or not       */
INT             ngline;      /* number of geometrylines to the element */
GLINE          *gline[4];    /* geometrylines of the element           */
NEUM_CONDITION *lineneum[4]; /* short grep on line-neum. conditions    */
INT             line;        /* looper over lines                      */
INT             ngnode;      /* number of geometry-nodes on g-line     */
INT             ngr,ngs;     /* number of GP'e for line-integration    */
DOUBLE          xgp[3];      /* Coordinates of GP'es                   */
DOUBLE          ygp[3];      /* Coordinates of GP'es                   */
DOUBLE          wgx[3];      /* weights at GP'es                       */
DOUBLE          wgy[3];      /* weights at GP'es                       */
DOUBLE          ds;       /* dx/dr line increment for line integration */
DOUBLE          facline;     /*integration factor for line integration */
DOUBLE          forceline[2];/* lineload value in x and y direction(inp)*/
RSF rsgeo;                   /* integration direction on line          */

#ifdef DEBUG 
dstrc_enter("wge_eleload");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* init phase        (init=1)                                          */
/*---------------------------------------------------------------------*/
if (init==1)
{
  eload   = amdef("eload"  ,&eload_a,MAXDOFPERNODE,MAXNOD_WALL1,"DA");
  functd  = amdef("functd"  ,&functd_a,MAXNOD_WALL1,1          ,"DV");       
  derivd  = amdef("derivd"  ,&derivd_a,2          ,MAXNOD_WALL1,"DA");       
  xjm     = amdef("xjm_a"  ,&xjm_a  ,numdfd        ,numdfd     ,"DA");
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
amdel(&eload_a);
amdel(&functd_a);
amdel(&derivd_a);
amdel(&xjm_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*------------------------------------------- get integration parameter */
wgeintg(ele,data,1);
nir     = ele->e.wallge->nGP[0];
nis     = ele->e.wallge->nGP[1];
ield    = ele->numnp;
nd      = ield * numdfd;
/*--------------------------------- check for presence of surface loads */
foundsurface=0;
if (!(ele->g.gsurf->neum)) goto endsurface;
if (ele->g.gsurf->neum->neum_type!=pres_domain_load)
{
  dserror("load case not implemented");
  goto endsurface;
}
foundsurface=1;
for (lr=0; lr<nir; lr++)/*------------------------- loop r-direction */
{
  /*============================ gaussian point and weight at it ===*/
  e1   = data->xgrr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)/*---------------------- loop s direction */
  {
     /*========================= gaussian point and weight at it ===*/
     e2   = data->xgss[ls];
     facs = data->wgts[ls];
     /*----------------------------------------- shape functions ---*/
     w1_funct_deriv(functd,derivd,e1,e2,ele->distyp,1);
     /*------------------------------------------ jacobian matrix ---*/       
     w1_jaco(functd,derivd,xjm,&det,ele,ield); 
     /*--------------------------------------- integration factor ---*/ 
     fac = facr * facs * det; 
     /*--------------------------------------------------------------*/ 
     w1_fextsurf(ele,eload,functd,fac,ield); 
  }
}
/*----------------------------------------------------------------------*/
endsurface:
/*----------------------------------------------------------------------*/

/*--------- integration of line loads on lines adjacent to this element */
foundline=0;
/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
/*----------------------------------------------------------------------*/
/*-------------- loop over lines, check for neumann conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   lineneum[i] = gline[i]->neum;
   if (lineneum[i]==NULL) continue;
   if (lineneum[i]) foundline=1; 
}
if (foundline==0) goto endline;
/*-------------------------- loop over lines (with neumann conditions) */
for (line=0; line<ngline; line++)
{
   if (lineneum[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*--------- original GP-coordinates and weights for area-integration */
   ngr = nir; 
   ngs = nis; 
    for (i=0; i<ngr; i++) 
    { 
      xgp[i] = data->xgrr[i];
      wgx[i] = data->wgtr[i]; 
     }
    for (i=0; i<ngs; i++) 
    { 
       ygp[i] = data->xgss[i];
       wgy[i] = data->wgts[i]; 
     }
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
         /*------------------------------------------- shape functions ---*/
          w1_funct_deriv(functd,derivd,e1,e2,ele->distyp,1);
          /*--------------------------------------------- jacobian matrix */       
          w1_jaco(functd,derivd,xjm,&det,ele,ield); 
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
          for (i=0; i<numdfd; i++)
          {
             forceline[i] = gline[line]->neum->neum_val.a.dv[i];
          }
          for (j=0; j<ield; j++)
          {
              for (i=0; i<numdfd; i++)
              {
                 eload[i][j] += functd[j] * forceline[i] * facline;
              }
          }
          
          
        }/*========================================== end of loop over ls */
      }/*============================================= end of loop over lr */
   /* line number lie has been done,switch of the neumann pointer of it */
   ele->g.gsurf->gline[line]->neum=NULL;
}/* end loop line over lines */
/*----------------------------------------------------------------------*/
endline:
/*----------------------------------------------------------------------*/
/* add static array eload to global external element load vector        */
/*----------------------------------------------------------------------*/
if (foundline != 0 || foundsurface != 0)
{
  for (nodei=0; nodei<4; nodei++)
  {
    nodestarti = 3* nodei;
    loadvec[nodestarti]   = eload[0][nodei];
    loadvec[nodestarti+1] = eload[1][nodei];
    loadvec[nodestarti+2] = 0.0;
  }
  if(ield > 4)
  {
    for (nodei=4; nodei<ield; nodei++)
    {
      nodestarti = 2* nodei + 4;
      loadvec[nodestarti]   = eload[0][nodei];
      loadvec[nodestarti+1] = eload[1][nodei];
   }
  }
}
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of wge_eleload */

#endif /*D_WALLGE*/

/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
