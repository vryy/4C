/*!----------------------------------------------------------------------
\file
\brief contains the routine

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../interf/interf.h"
#include "../interf/interf_prototypes.h"
#include "s2ml.h"
#include "s2ml_prototypes.h"


/*!
\addtogroup MLSTRUCT
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine for calculation of

\param *actpart         PARTITION   (i)   my partition
\param *container       CONTAINER   (i/o) contains variables defined in container.h

\return void
\sa calling:   ifinit, ifstatic_ke, if_cal_stress, if_eleload;
    called by: global_calelm();

*----------------------------------------------------------------------*/
void s2ml_stiff_interf(MATERIAL  *actsmmat,    /* actual submesh material*/
                       ELEMENT   *actmaele,    /* actual macro element   */
                       ELEMENT   *actsmele,    /* actual submesh element */
                       ARRAY     *estif_ma_ma, /* element stiffness macro-macro  */
                       ARRAY     *estif_ma_mi, /* element stiffness macro-micro  */
                       ARRAY     *estif_mi_ma, /* element stiffness micro-macro  */
                       ARRAY     *estif_mi_mi, /* element stiffness micro-micro  */
                       ARRAY     *intforce_ma, /* element internal force macro   */
                       ARRAY     *intforce_mi, /* element internal force micro   */
                       INT        istore,      /* update?                */
                       INT        init)        /* allocate arrays?       */
{
const DOUBLE    q12 = ONE/TWO;

INTERF_DATA     actsmdata;

INT       i,lr,k;    /* loopers */
INT       ielmi;     /* numnp to this submesh-element     */
INT       numdfmi;   /* numdof of this submesh-element     */
INT       flag=0;    /* flag for distinction between differnet element orientation cases  */
INT       nir;       /* number of gaussian points            */
INT       ip;        /* GP-counter */
INT       smallscale=1;/* flag for mat law->it's call by smelement in multiscale */
DOUBLE    c_parabel=0.0; /* y = a + b*x + c*x^2  */
DOUBLE    b_parabel=0.0; /* y = a + b*x + c*x^2  */
DOUBLE    help;      /* Zwischenwert */
DOUBLE    width;     /* element thickness in wall plane                    */
DOUBLE    Thick;     /* element thickness perpendicular to wall plane    */
DOUBLE    e1;        /* xi-coordinate of gaussian point            */
DOUBLE    facr;      /* weight at gaussian point                         */
DOUBLE    fac;       /* integration factor    */
DOUBLE    co,si;     /* local orientation */
DOUBLE    det;       /* determinants of jacobian matrix  */

static ARRAY    xrefe_a;     /* coordinates of element nodes */
static DOUBLE **xrefe;
static ARRAY    functmi_a;     /* shape functions for [u]*/
static DOUBLE  *functmi;
static ARRAY    bopmi_a;     /* Micro-"B-operator"=A N'L for nt-direction*/
static DOUBLE **bopmi;
static ARRAY    bopma_a;     /* Macro-"B-operator"=A N'L V for nt-direction*/
static DOUBLE **bopma;
static ARRAY    D_a;         /* material tangente */
static DOUBLE **D;

static DOUBLE **stiff_ma_ma;          /* element stiffness macro-macro */
static DOUBLE **stiff_ma_mi;          /* element stiffness macro-micro */
static DOUBLE **stiff_mi_ma;          /* element stiffness micro-macro */
static DOUBLE **stiff_mi_mi;          /* element stiffness micro-micro */
static DOUBLE  *fint_ma;              /* element internal force macro  */
static DOUBLE  *fint_mi;              /* element internal force micro  */

DOUBLE L[4];                 /* distance between corner nodes of the element */
DOUBLE x_mid[3];             /* x-coordinates on xi-axis */
DOUBLE y_mid[3];             /* y-coordinates on xi-axis */
DOUBLE jumpuma[2];           /* displacement jumps (tn) due to macro-displ. */
DOUBLE jumpumi[2];           /* displacement jumps (tn) due to micro-displ. */
DOUBLE DELTAjumpuma[2];      /* incremental displacement jumps (tn) due to macro-displ. */
DOUBLE DELTAjumpumi[2];      /* incremental displacement jumps (tn) due to micro-displ. */
DOUBLE jumpu_tot[2];         /* total displacement jumps (tn)  */
DOUBLE DELTAjumpu_tot[2];    /* total incremental displacement jumps (tn)  */
DOUBLE T[2];                 /* stress */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_stiff_interf");
#endif
/*----------------------------------------------------------------------*/
if (init==1)
{
xrefe     = amdef("xrefe"  ,&xrefe_a,  2,8, "DA");
functmi   = amdef("functmi",&functmi_a,3,1, "DV");
bopmi     = amdef("bopmi"  ,&bopmi_a  ,2,16,"DA");
bopma     = amdef("bopma"  ,&bopma_a  ,2,8, "DA");
D         = amdef("D"      ,&D_a      ,2,2, "DA");

goto end;
}
ielmi     = actsmele->numnp;
numdfmi   = 2* ielmi;
/*----------- check orientation of element (which is my xi direction)---*/
for (k=0; k<ielmi; k++)
{
 xrefe[0][k] = actsmele->node[k]->x[0];   /* coordinates in x-direction */
 xrefe[1][k] = actsmele->node[k]->x[1];   /* coordinates in y-direction */
}
L[0] = sqrt( (xrefe[0][1] - xrefe[0][0]) * (xrefe[0][1] - xrefe[0][0])
     +       (xrefe[1][1] - xrefe[1][0]) * (xrefe[1][1] - xrefe[1][0]));
L[1] = sqrt( (xrefe[0][2] - xrefe[0][1]) * (xrefe[0][2] - xrefe[0][1])
     +       (xrefe[1][2] - xrefe[1][1]) * (xrefe[1][2] - xrefe[1][1]));

/*-------------------------------------------- integration parameters ---*/
ifintg(actsmele,&actsmdata);

/*--------------------- coordinates of "nonexisting nodes" on xi-axes ---*/
switch(actsmele->distyp)
{
case quad4:

     if(L[0]>L[1])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][3]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][3]);
      x_mid[1]   = q12*(xrefe[0][1] + xrefe[0][2]);
      y_mid[1]   = q12*(xrefe[1][1] + xrefe[1][2]);
      flag = 1;
      width=L[1];
     }
     else if (L[1]>L[0])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][1]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][1]);
      x_mid[1]   = q12*(xrefe[0][2] + xrefe[0][3]);
      y_mid[1]   = q12*(xrefe[1][2] + xrefe[1][3]);
      flag = 2;
      width=L[0];
     }

break;
/*-----------------------------------------------------------------------*/
case quad8:

     if(L[0]>L[1])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][3]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][3]);
      x_mid[1]   = q12*(xrefe[0][1] + xrefe[0][2]);
      y_mid[1]   = q12*(xrefe[1][1] + xrefe[1][2]);
      x_mid[2]   = q12*(xrefe[0][4] + xrefe[0][6]);
      y_mid[2]   = q12*(xrefe[1][4] + xrefe[1][6]);
      flag = 1;
      width=L[1];
     }
     else if (L[1]>L[0])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][1]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][1]);
      x_mid[1]   = q12*(xrefe[0][2] + xrefe[0][3]);
      y_mid[1]   = q12*(xrefe[1][2] + xrefe[1][3]);
      x_mid[2]   = q12*(xrefe[0][5] + xrefe[0][7]);
      y_mid[2]   = q12*(xrefe[1][5] + xrefe[1][7]);
      flag = 2;
      width=L[0];
     }
     help      = (x_mid[0]-x_mid[1])/(x_mid[0]-x_mid[2]);
     c_parabel = (y_mid[0]-y_mid[1]-(y_mid[0]-y_mid[2])*help)/
                 (x_mid[0]*x_mid[0]-x_mid[1]*x_mid[1]-
                  (x_mid[0]*x_mid[0]-x_mid[2]*x_mid[2])*help);
     b_parabel = (y_mid[0]-y_mid[1]-c_parabel*(x_mid[0]*x_mid[0]-x_mid[1]*x_mid[1]))/
                 (x_mid[0]-x_mid[1]);

break;
   default:
     dserror("unkonwn discretisation for interface");
   break;
}
Thick = actsmele->e.interf->thick;
nir   = actsmele->e.interf->nGP;
/*---- stifness and internal forces have to be reinitialized to zero ---*/
amzero(estif_ma_ma);
amzero(estif_ma_mi);
amzero(estif_mi_ma);
amzero(estif_mi_mi);
amzero(intforce_ma);
amzero(intforce_mi);
/*------------------------------------------------------- short cuts ---*/
stiff_ma_ma     = estif_ma_ma->a.da;
stiff_ma_mi     = estif_ma_mi->a.da;
stiff_mi_ma     = estif_mi_ma->a.da;
stiff_mi_mi     = estif_mi_mi->a.da;
fint_ma         = intforce_ma->a.dv;
fint_mi         = intforce_mi->a.dv;
/*=======================================================================*/
ip = -1;
/*================================================== integration loop ===*/
for (lr=0; lr<nir; lr++)
{
  ip++;
   /*======================== submesh gaussian point and weight at it ===*/
   e1   = actsmdata.xgr[lr];
   facr = actsmdata.wgtr[lr];
   /*----------------------- submesh shape functions, angle and detJ ---*/
   if_funcderiv(e1,actsmele->distyp,x_mid,y_mid,b_parabel,c_parabel,functmi,&co,&si,&det);
   /*------------------------------------ submesh integration factor ---*/
   fac  = facr * det * Thick;
   /*--------------------- calculate submesh micro-operator B'=A N'L ---*/
   amzero(&bopmi_a);
   if_bop(actsmele->distyp,bopmi,functmi,co,si,flag);
   /*----------------- calculate submesh macro-operator Bbar=A N'L V ---*/
   amzero(&bopma_a);
   s2ml_ifbopmajumpu(actsmele,actmaele,bopmi,bopma,jumpuma,DELTAjumpuma,jumpumi,DELTAjumpumi);
   /*-------------------------------------- total displacement jumps ---*/
   for (i=0; i<2; i++)
   {
     jumpu_tot[i]      = jumpuma[i] + jumpumi[i];
     DELTAjumpu_tot[i] = DELTAjumpuma[i] + DELTAjumpumi[i];
   }
  /*------------------------ call material law -> stress and tangent ---*/
   switch(actsmmat->mattyp)
   {
   case m_ifmat:
         if_mat(actmaele,actsmmat,NULL,D,T,ip,istore,0,smallscale,actsmele,jumpu_tot,DELTAjumpu_tot);
   break;
   case m_interf_therm:
       if_mat_thermodyn(actmaele,actsmmat,NULL,D,T,ip,istore,0,smallscale,actsmele,jumpu_tot,DELTAjumpu_tot);
   break;
   default:
     dserror(" unknown type of material law");
   break;
   }
   /*-------------------------------------------------------------------*/
   /*-------------------------- element stiffness matrix macro macro ---*/
   if_ke(4,flag,stiff_ma_ma,bopma,D,fac);
   /*-------------------------- element stiffness matrix macro micro ---*/
   s2ml_cal_stiff(stiff_ma_mi,bopma,bopmi,D,fac,8,numdfmi,2);
   /*-------------------------- element stiffness matrix micro macro ---*/
   s2ml_cal_stiff(stiff_mi_ma,bopmi,bopma,D,fac,numdfmi,8,2);
   /*-------------------------- element stiffness matrix micro micro ---*/
   if_ke(ielmi,flag,stiff_mi_mi,bopmi,D,fac);
   /*---------------------------------- element internal force macro ---*/
   if_fint(4,T,fac,bopma,fint_ma);
   /*---------------------------------- element internal force micro ---*/
   if_fint(ielmi,T,fac,bopmi,fint_mi);
}/*============================================= end of loop over lr ===*/

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of s2ml_stiff_interf */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
#endif
