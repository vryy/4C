/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | integration of mass matrix                                           |
 | for shell8 element                                     m.gee 4/03    |
 *----------------------------------------------------------------------*/
void s8static_mass(ELEMENT   *ele,                         /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   DOUBLE    *force,                      /* for internal forces (initialized!) */
                   INT        kstep,                      /* actual step in nonlinear analysis */
                   INT        init)                       /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*/
/* if force==NULL no internal forces are calculated                     */
/*----------------------------------------------------------------------*/
{
INT                 k;                                      /* some loopers */
INT                 ngauss;
INT                 imass;                                  /* flag for calculating mass matrix */
DOUBLE              density;                                /* density of the material */
DOUBLE              facv,facw,facvw;                        /* variables for mass integration */
DOUBLE             *thick;

INT                 nir,nis,nit;                            /* num GP in r/s/t direction */
INT                 lr, ls, lt;                             /* loopers over GP */
INT                 iel;                                    /* numnp to this element */
INT                 nd;                                     /* ndofs to this element (=numdf*numnp) */
const INT           numdf=NUMDOF_SHELL8;                    /* ndofs per node to this element */

DOUBLE              e1,e2,e3;                               /*GP-coords*/
DOUBLE              fac,facr,facs,fact;                     /* weights at GP */
DOUBLE              xnu;                                    /* value of shell shifter */

DOUBLE              condfac;                                /* sdc conditioning factor */
DOUBLE              h2;                                     /* half nodal height */

DOUBLE              detr;                                   /* jacobian determinant in ref conf. at gp */
DOUBLE              detc;                                   /* jacobian determinant in cur conf. at gp */

DOUBLE              h[3];                                   /* working array */
DOUBLE              da;                                     /* area on mid surface */

static ARRAY        a3r_a;       static DOUBLE **a3r;       /* a3 in reference config (lenght h/2) */
static ARRAY        a3c_a;       static DOUBLE **a3c;       /* a3 in current   config (lenght h/2 + disp) */

static ARRAY        a3kvpr_a;    static DOUBLE **a3kvpr;    /* partiel derivatives of normal vector ref.config. */
static ARRAY        a3kvpc_a;    static DOUBLE **a3kvpc;    /* partiel derivatives of normal vector cur.config. */

static ARRAY        xrefe_a;     static DOUBLE **xrefe;     /* coords of midsurface in ref config */
static ARRAY        xcure_a;     static DOUBLE **xcure;     /* coords of midsurface in cur condig */
                                        DOUBLE **a3ref;     /* elements directors (lenght 1) */

static ARRAY        funct_a;     static DOUBLE  *funct;     /* shape functions */
static ARRAY        deriv_a;     static DOUBLE **deriv;     /* derivatives of shape functions */

/* mid surface basis vectors and metric tensors */
static ARRAY        akovr_a;     static DOUBLE **akovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        akonr_a;     static DOUBLE **akonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        amkovr_a;    static DOUBLE **amkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        amkonr_a;    static DOUBLE **amkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        akovc_a;     static DOUBLE **akovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        akonc_a;     static DOUBLE **akonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        amkovc_a;    static DOUBLE **amkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        amkonc_a;    static DOUBLE **amkonc;    /* kontravar.--------------"------------ current.config. */

/* shell body basis vectors and metric tensors */
static ARRAY        gkovr_a;     static DOUBLE **gkovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        gkonr_a;     static DOUBLE **gkonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        gmkovr_a;    static DOUBLE **gmkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        gmkonr_a;    static DOUBLE **gmkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        gkovc_a;     static DOUBLE **gkovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        gkonc_a;     static DOUBLE **gkonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        gmkovc_a;    static DOUBLE **gmkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        gmkonc_a;    static DOUBLE **gmkonc;    /* kontravar.--------------"------------ current.config. */

                                 static DOUBLE **emass;     /* element mass matrix */



#ifdef DEBUG
dstrc_enter("s8static_mass");
#endif
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
if (init==1)
{
iel = MAXNOD_SHELL8; /* maximum number of nodes for this type of shell */

xrefe     = amdef("xrefe"  ,&xrefe_a,3,MAXNOD_SHELL8,"DA");
xcure     = amdef("xcure"  ,&xcure_a,3,MAXNOD_SHELL8,"DA");
a3r       = amdef("a3r"    ,&a3r_a,3,MAXNOD_SHELL8,"DA");
a3c       = amdef("a3c"    ,&a3c_a,3,MAXNOD_SHELL8,"DA");

a3kvpr    = amdef("a3kvpr" ,&a3kvpr_a,3,2,"DA");
a3kvpc    = amdef("a3kvpc" ,&a3kvpc_a,3,2,"DA");

funct     = amdef("funct"  ,&funct_a,MAXNOD_SHELL8,1,"DV");
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL8,"DA");

akovr     = amdef("akovr"  ,&akovr_a,3,3,"DA");
akonr     = amdef("akonr"  ,&akonr_a,3,3,"DA");
amkovr    = amdef("amkovr" ,&amkovr_a,3,3,"DA");
amkonr    = amdef("amkonr" ,&amkonr_a,3,3,"DA");

akovc     = amdef("akovc"  ,&akovc_a,3,3,"DA");
akonc     = amdef("akonc"  ,&akonc_a,3,3,"DA");
amkovc    = amdef("amkovc" ,&amkovc_a,3,3,"DA");
amkonc    = amdef("amkonc" ,&amkonc_a,3,3,"DA");

gkovr     = amdef("gkovr"  ,&gkovr_a,3,3,"DA");
gkonr     = amdef("gkonr"  ,&gkonr_a,3,3,"DA");
gmkovr    = amdef("gmkovr" ,&gmkovr_a,3,3,"DA");
gmkonr    = amdef("gmkonr" ,&gmkonr_a,3,3,"DA");

gkovc     = amdef("gkovc"  ,&gkovc_a,3,3,"DA");
gkonc     = amdef("gkonc"  ,&gkonc_a,3,3,"DA");
gmkovc    = amdef("gmkovc" ,&gmkovc_a,3,3,"DA");
gmkonc    = amdef("gmkonc" ,&gmkonc_a,3,3,"DA");

goto end;
}
/*----------------------------------------------------------------------*/
/* calculation phase                                                    */
/*----------------------------------------------------------------------*/
/*-------------------------------------------- init the gaussian points */
s8intg(ele,data,0);
/*------------------------------------ check calculation of mass matrix */
imass = 1;
amzero(emass_global);
emass = emass_global->a.da;
s8_getdensity(mat,&density);
thick = ele->e.s8->thick_node.a.dv;
/*----------------------------------------------- integrationsparameter */
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
nit     = ele->e.s8->nGP[2];
iel     = ele->numnp;
nd      = iel*NUMDOF_SHELL8;
condfac = ele->e.s8->sdc;
a3ref   = ele->e.s8->a3ref.a.da;
/*----------------------------------------------------- geometry update */
for (k=0; k<iel; k++)
{
   h2 = ele->e.s8->thick_node.a.dv[k];
   h2 /= 2.0;
   h2 *= condfac;

   a3r[0][k] = a3ref[0][k] * h2;
   a3r[1][k] = a3ref[1][k] * h2;
   a3r[2][k] = a3ref[2][k] * h2;

   xrefe[0][k] = ele->node[k]->x[0];
   xrefe[1][k] = ele->node[k]->x[1];
   xrefe[2][k] = ele->node[k]->x[2];

   xcure[0][k] = xrefe[0][k] + ele->node[k]->sol.a.da[0][0];
   xcure[1][k] = xrefe[1][k] + ele->node[k]->sol.a.da[0][1];
   xcure[2][k] = xrefe[2][k] + ele->node[k]->sol.a.da[0][2];

   a3c[0][k] = a3r[0][k]     + ele->node[k]->sol.a.da[0][3];
   a3c[1][k] = a3r[1][k]     + ele->node[k]->sol.a.da[0][4];
   a3c[2][k] = a3r[2][k]     + ele->node[k]->sol.a.da[0][5];
}
/*=================================================== integration loops */
ngauss=0;
for (lr=0; lr<nir; lr++)
{
   /*================================== gaussian point and weight at it */
   e1   = data->xgpr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      /*=============================== gaussian point and weight at it */
      e2   = data->xgps[ls];
      facs = data->wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      s8_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------------------------------------ init mass matrix variables */
         facv  = 0.0;
         facw  = 0.0;
         facvw = 0.0;
      /*------------------------------------- metrics at gaussian point */
      s8_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,&detr,
                 funct,deriv,iel,a3kvpr,0);
      s8_tvmr(xcure,a3c,akovc,akonc,amkovc,amkonc,&detc,
                 funct,deriv,iel,a3kvpc,0);
      /*------------------------- make h as cross product in ref config.
                                       to get area da on shell mid surf */
      h[0] = akovr[1][0]*akovr[2][1] - akovr[2][0]*akovr[1][1];
      h[1] = akovr[2][0]*akovr[0][1] - akovr[0][0]*akovr[2][1];
      h[2] = akovr[0][0]*akovr[1][1] - akovr[1][0]*akovr[0][1];
      /*------------------------------------- make director unit lenght
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);
      /*============================== loop GP in thickness direction t */
      for (lt=0; lt<nit; lt++)
      {
         /*---------------------------- gaussian point and weight at it */
         e3   = data->xgpt[lt];
         fact = data->wgtt[lt];
         /*-------------------- basis vectors and metrics at shell body */
         s8_tmtr(xrefe,a3r,e3,gkovr,gkonr,gmkovr,gmkonr,&detr,
                    funct,deriv,iel,condfac,0);

         s8_tmtr(xcure,a3c,e3,gkovc,gkonc,gmkovc,gmkonc,&detc,
                    funct,deriv,iel,condfac,0);
         /*--------------------------------- metric at gp in shell body */
         s8_tvhe(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,&detr,&detc,
                 amkovc,amkovr,akovc,akovr,a3kvpc,a3kvpr,e3,condfac);
         /*----------- calc shell shifter and put it in the weight fact */
         xnu   = (1.0/condfac)*(detr/da);
         fact *= xnu;
         /*-------------------------- mass matrix thickness integration */
         facv  += data->wgtt[lt] * detr;
         facw  += data->wgtt[lt] * detr * DSQR(e3);
         facvw += data->wgtt[lt] * detr * e3;
      }/*========================================== end of loop over lt */
      /*------------- mass matrix : gaussian point on shell mid surface */
      if (imass)
      {
         fac    = facr * facs * density;
         facv  *= fac;
         facw  *= fac;
         facvw *= fac;
         s8_tmas(funct,thick,emass,iel,numdf,facv,facw,facvw);
      }
   ngauss++;
   }/*============================================= end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------- end */
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8static_mass */
#endif
