/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1oint' which calclate displacement
       derivatives for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/
#ifdef D_BRICK1
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../brick1/brick1.h"
#include "../brick1/brick1_prototypes.h"
#include "../headers/optimization.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;
/*!----------------------------------------------------------------------
\brief integration routine for BRICK1 element

<pre>                                                              al 06/02
This routine performs integration of an 3D-hex-element.

</pre>
\param           *ele ELEMENT  (i)   element data
\param          *data C1_DATA  (i)   hex element data
\param           *mat MATERIAL (i)   material data
\param  *estif_global ARRAY    (o)   element stiffness matrix
\param         *force DOUBLE   (o)   vector for internal forces
\param         *init  INT      (i)   flag for initialization (alloc mem...)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_oint()

*----------------------------------------------------------------------*/
void c1_oint(
             ELEMENT   *ele,
             C1_DATA   *data,
             MATERIAL  *mat,
             DOUBLE    *retval,  /* return value */
             INT        init     /* ==2 calc.strain energy */
             )
{
INT                 i,j;              /* some loopers */
INT                 nir,nis,nit;      /* num GP in r/s/t direction */
INT                 lr, ls, lt;       /* loopers over GP */
INT                 ip;
INT                 iel;              /* numnp to this element */
INT                 nd;
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 0;/* controls evaluation of new stresses    */
const INT           numdf =3;
const INT           numeps=6;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE              rs;

DOUBLE dfie[80];
DOUBLE disd[9];
DOUBLE F[6];  /* element stress vector   (stress-resultants) */
DOUBLE DF[6]; /* derivative of element stress vector         */
DOUBLE strain[6];
DOUBLE xyze[60];
DOUBLE edis[60];
DOUBLE grdis[60];/* displacement derivatives                  */
DOUBLE  g[6][6]; /* transformation matrix s(glob)= g*s(loc)   */
DOUBLE gi[6][6]; /* inverse of g          s(loc) = gi*s(glob) */
/*-------------------------------------------  for eas elements only ---*/
INT    cc;

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    DD_a;      /* material tensor */
static DOUBLE **DD;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    bnop_a;   /* BN-operator */
static DOUBLE **bn;

static DOUBLE **estiflo;
static ARRAY    estiflo_a; /* local element stiffness matrix ke for eas */

DOUBLE det;

INT    iform;             /* index for nonlinear formulation of element */
INT    calstr;            /* flag for stress calculation                */

DOUBLE  ste;  /* strain energy               */
DOUBLE dste;  /* derivative of strain energy */
DOUBLE  stm;  /* mass                        */
DOUBLE  stv;  /* volume                      */
DOUBLE  stf;  /* frequency                   */

INT    ieig;     /* counter for eigen values   */
DOUBLE kro, kru; /* Kreiselmeier Steinhaeuser  */
DOUBLE reig[10][60];/* vec. with eigenforms    */
DOUBLE teta_s1[10]; /* helping vector for freq */
DOUBLE teta1[10];   /* helping vector for freq */
DOUBLE hvar1[60];   /* helping vector for freq */
DOUBLE w2s1[10];    /* helping vector for freq */
DOUBLE krexpo;      /* exponent in Kreiselmeier*/
DOUBLE dens, density;  /* density                     */
DOUBLE lmvec[60];   /* lumped mass vector     */
DOUBLE facm, totmas, emasdg;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_oint");
#endif
/*----------------------------------------------------------------------*/
istore = 0;
calstr = 0;
newval = 0;
/*------------------------------------------------- some working arrays */
if (init==1)
{
  funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
  deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
  D         = amdef("D"      ,&D_a   ,6,6              ,"DA");
  DD        = amdef("DD"     ,&DD_a  ,6,6              ,"DA");
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");
  bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
  bn        = amdef("bnop" ,&bnop_a,3     ,       MAXNOD_BRICK1 ,"DA");
  estiflo   = amdef("estiflo"  ,&estiflo_a ,60,60,"DA");
goto end;
}
/*------------------------------------------- integration parameters ---*/
c1intg(ele,data);
/*------------------------------------------- integration parameters ---*/
  nir     = ele->e.c1->nGP[0];
  nis     = ele->e.c1->nGP[1];
  nit     = ele->e.c1->nGP[2];
  iel     = ele->numnp;
  nd      = numdf * iel;
/*---------------------------------- setup individual element arrays ---*/
  if(iel==8)
  {
    cc=0;
    for (i=0;i<iel;i++) for (j=0;j<3;j++) xyze[cc++] = ele->node[i]->x[j];
    cc=0;
    for (i=0;i<iel;i++) for (j=0;j<3;j++) edis[cc++] = ele->node[i]->sol.a.da[0][j];
    /*------------------------------------------------------------------*/
    if (init==6)
    {
    cc=0;
    for (i=0;i<iel;i++) for (j=0;j<3;j++) grdis[cc++] = ele->node[i]->sol.a.da[1][j];
    }
    /*------------------------------------------------------------------*/
    if (init==7)/* eigen frequency optimization */
    {
      for (ieig=0;ieig<opt->oeig->numeigv;ieig++)
      {/* loop eigenvalues */
        cc=0;
        for (i=0;i<iel;i++) for (j=0;j<3;j++)
                    reig[ieig][cc++] = ele->node[i]->sol.a.da[ieig][j];
      }/* loop eigenvalues */
    }
    /*------------------------------------------------------------------*/
  }
  else
  {/*iel==20*/
   cc=0;
   xyze[cc++] = ele->node[0]->x[0];
   xyze[cc++] = ele->node[0]->x[1];
   xyze[cc++] = ele->node[0]->x[2];
   xyze[cc++] = ele->node[1]->x[0];
   xyze[cc++] = ele->node[1]->x[1];
   xyze[cc++] = ele->node[1]->x[2];
   xyze[cc++] = ele->node[2]->x[0];
   xyze[cc++] = ele->node[2]->x[1];
   xyze[cc++] = ele->node[2]->x[2];
   xyze[cc++] = ele->node[3]->x[0];
   xyze[cc++] = ele->node[3]->x[1];
   xyze[cc++] = ele->node[3]->x[2];
   xyze[cc++] = ele->node[4]->x[0];
   xyze[cc++] = ele->node[4]->x[1];
   xyze[cc++] = ele->node[4]->x[2];
   xyze[cc++] = ele->node[5]->x[0];
   xyze[cc++] = ele->node[5]->x[1];
   xyze[cc++] = ele->node[5]->x[2];
   xyze[cc++] = ele->node[6]->x[0];
   xyze[cc++] = ele->node[6]->x[1];
   xyze[cc++] = ele->node[6]->x[2];
   xyze[cc++] = ele->node[7]->x[0];
   xyze[cc++] = ele->node[7]->x[1];
   xyze[cc++] = ele->node[7]->x[2];
   xyze[cc++] = ele->node[8]->x[0];
   xyze[cc++] = ele->node[8]->x[1];
   xyze[cc++] = ele->node[8]->x[2];
   xyze[cc++] = ele->node[9]->x[0];
   xyze[cc++] = ele->node[9]->x[1];
   xyze[cc++] = ele->node[9]->x[2];
   xyze[cc++] = ele->node[10]->x[0];
   xyze[cc++] = ele->node[10]->x[1];
   xyze[cc++] = ele->node[10]->x[2];
   xyze[cc++] = ele->node[11]->x[0];
   xyze[cc++] = ele->node[11]->x[1];
   xyze[cc++] = ele->node[11]->x[2];
   xyze[cc++] = ele->node[16]->x[0];
   xyze[cc++] = ele->node[16]->x[1];
   xyze[cc++] = ele->node[16]->x[2];
   xyze[cc++] = ele->node[17]->x[0];
   xyze[cc++] = ele->node[17]->x[1];
   xyze[cc++] = ele->node[17]->x[2];
   xyze[cc++] = ele->node[18]->x[0];
   xyze[cc++] = ele->node[18]->x[1];
   xyze[cc++] = ele->node[18]->x[2];
   xyze[cc++] = ele->node[19]->x[0];
   xyze[cc++] = ele->node[19]->x[1];
   xyze[cc++] = ele->node[19]->x[2];
   xyze[cc++] = ele->node[12]->x[0];
   xyze[cc++] = ele->node[12]->x[1];
   xyze[cc++] = ele->node[12]->x[2];
   xyze[cc++] = ele->node[13]->x[0];
   xyze[cc++] = ele->node[13]->x[1];
   xyze[cc++] = ele->node[13]->x[2];
   xyze[cc++] = ele->node[14]->x[0];
   xyze[cc++] = ele->node[14]->x[1];
   xyze[cc++] = ele->node[14]->x[2];
   xyze[cc++] = ele->node[15]->x[0];
   xyze[cc++] = ele->node[15]->x[1];
   xyze[cc++] = ele->node[15]->x[2];
   cc=0;
   edis[cc++] = ele->node[0]->sol.a.da[0][0];
   edis[cc++] = ele->node[0]->sol.a.da[0][1];
   edis[cc++] = ele->node[0]->sol.a.da[0][2];
   edis[cc++] = ele->node[1]->sol.a.da[0][0];
   edis[cc++] = ele->node[1]->sol.a.da[0][1];
   edis[cc++] = ele->node[1]->sol.a.da[0][2];
   edis[cc++] = ele->node[2]->sol.a.da[0][0];
   edis[cc++] = ele->node[2]->sol.a.da[0][1];
   edis[cc++] = ele->node[2]->sol.a.da[0][2];
   edis[cc++] = ele->node[3]->sol.a.da[0][0];
   edis[cc++] = ele->node[3]->sol.a.da[0][1];
   edis[cc++] = ele->node[3]->sol.a.da[0][2];
   edis[cc++] = ele->node[4]->sol.a.da[0][0];
   edis[cc++] = ele->node[4]->sol.a.da[0][1];
   edis[cc++] = ele->node[4]->sol.a.da[0][2];
   edis[cc++] = ele->node[5]->sol.a.da[0][0];
   edis[cc++] = ele->node[5]->sol.a.da[0][1];
   edis[cc++] = ele->node[5]->sol.a.da[0][2];
   edis[cc++] = ele->node[6]->sol.a.da[0][0];
   edis[cc++] = ele->node[6]->sol.a.da[0][1];
   edis[cc++] = ele->node[6]->sol.a.da[0][2];
   edis[cc++] = ele->node[7]->sol.a.da[0][0];
   edis[cc++] = ele->node[7]->sol.a.da[0][1];
   edis[cc++] = ele->node[7]->sol.a.da[0][2];
   edis[cc++] = ele->node[8]->sol.a.da[0][0];
   edis[cc++] = ele->node[8]->sol.a.da[0][1];
   edis[cc++] = ele->node[8]->sol.a.da[0][2];
   edis[cc++] = ele->node[9]->sol.a.da[0][0];
   edis[cc++] = ele->node[9]->sol.a.da[0][1];
   edis[cc++] = ele->node[9]->sol.a.da[0][2];
   edis[cc++] = ele->node[10]->sol.a.da[0][0];
   edis[cc++] = ele->node[10]->sol.a.da[0][1];
   edis[cc++] = ele->node[10]->sol.a.da[0][2];
   edis[cc++] = ele->node[11]->sol.a.da[0][0];
   edis[cc++] = ele->node[11]->sol.a.da[0][1];
   edis[cc++] = ele->node[11]->sol.a.da[0][2];
   edis[cc++] = ele->node[16]->sol.a.da[0][0];
   edis[cc++] = ele->node[16]->sol.a.da[0][1];
   edis[cc++] = ele->node[16]->sol.a.da[0][2];
   edis[cc++] = ele->node[17]->sol.a.da[0][0];
   edis[cc++] = ele->node[17]->sol.a.da[0][1];
   edis[cc++] = ele->node[17]->sol.a.da[0][2];
   edis[cc++] = ele->node[18]->sol.a.da[0][0];
   edis[cc++] = ele->node[18]->sol.a.da[0][1];
   edis[cc++] = ele->node[18]->sol.a.da[0][2];
   edis[cc++] = ele->node[19]->sol.a.da[0][0];
   edis[cc++] = ele->node[19]->sol.a.da[0][1];
   edis[cc++] = ele->node[19]->sol.a.da[0][2];
   edis[cc++] = ele->node[12]->sol.a.da[0][0];
   edis[cc++] = ele->node[12]->sol.a.da[0][1];
   edis[cc++] = ele->node[12]->sol.a.da[0][2];
   edis[cc++] = ele->node[13]->sol.a.da[0][0];
   edis[cc++] = ele->node[13]->sol.a.da[0][1];
   edis[cc++] = ele->node[13]->sol.a.da[0][2];
   edis[cc++] = ele->node[14]->sol.a.da[0][0];
   edis[cc++] = ele->node[14]->sol.a.da[0][1];
   edis[cc++] = ele->node[14]->sol.a.da[0][2];
   edis[cc++] = ele->node[15]->sol.a.da[0][0];
   edis[cc++] = ele->node[15]->sol.a.da[0][1];
   edis[cc++] = ele->node[15]->sol.a.da[0][2];
    if (init==6)
    {
   cc=0;
   grdis[cc++] = ele->node[0]->sol.a.da[1][0];
   grdis[cc++] = ele->node[0]->sol.a.da[1][1];
   grdis[cc++] = ele->node[0]->sol.a.da[1][2];
   grdis[cc++] = ele->node[1]->sol.a.da[1][0];
   grdis[cc++] = ele->node[1]->sol.a.da[1][1];
   grdis[cc++] = ele->node[1]->sol.a.da[1][2];
   grdis[cc++] = ele->node[2]->sol.a.da[1][0];
   grdis[cc++] = ele->node[2]->sol.a.da[1][1];
   grdis[cc++] = ele->node[2]->sol.a.da[1][2];
   grdis[cc++] = ele->node[3]->sol.a.da[1][0];
   grdis[cc++] = ele->node[3]->sol.a.da[1][1];
   grdis[cc++] = ele->node[3]->sol.a.da[1][2];
   grdis[cc++] = ele->node[4]->sol.a.da[1][0];
   grdis[cc++] = ele->node[4]->sol.a.da[1][1];
   grdis[cc++] = ele->node[4]->sol.a.da[1][2];
   grdis[cc++] = ele->node[5]->sol.a.da[1][0];
   grdis[cc++] = ele->node[5]->sol.a.da[1][1];
   grdis[cc++] = ele->node[5]->sol.a.da[1][2];
   grdis[cc++] = ele->node[6]->sol.a.da[1][0];
   grdis[cc++] = ele->node[6]->sol.a.da[1][1];
   grdis[cc++] = ele->node[6]->sol.a.da[1][2];
   grdis[cc++] = ele->node[7]->sol.a.da[1][0];
   grdis[cc++] = ele->node[7]->sol.a.da[1][1];
   grdis[cc++] = ele->node[7]->sol.a.da[1][2];
   grdis[cc++] = ele->node[8]->sol.a.da[1][0];
   grdis[cc++] = ele->node[8]->sol.a.da[1][1];
   grdis[cc++] = ele->node[8]->sol.a.da[1][2];
   grdis[cc++] = ele->node[9]->sol.a.da[1][0];
   grdis[cc++] = ele->node[9]->sol.a.da[1][1];
   grdis[cc++] = ele->node[9]->sol.a.da[1][2];
   grdis[cc++] = ele->node[10]->sol.a.da[1][0];
   grdis[cc++] = ele->node[10]->sol.a.da[1][1];
   grdis[cc++] = ele->node[10]->sol.a.da[1][2];
   grdis[cc++] = ele->node[11]->sol.a.da[1][0];
   grdis[cc++] = ele->node[11]->sol.a.da[1][1];
   grdis[cc++] = ele->node[11]->sol.a.da[1][2];
   grdis[cc++] = ele->node[16]->sol.a.da[1][0];
   grdis[cc++] = ele->node[16]->sol.a.da[1][1];
   grdis[cc++] = ele->node[16]->sol.a.da[1][2];
   grdis[cc++] = ele->node[17]->sol.a.da[1][0];
   grdis[cc++] = ele->node[17]->sol.a.da[1][1];
   grdis[cc++] = ele->node[17]->sol.a.da[1][2];
   grdis[cc++] = ele->node[18]->sol.a.da[1][0];
   grdis[cc++] = ele->node[18]->sol.a.da[1][1];
   grdis[cc++] = ele->node[18]->sol.a.da[1][2];
   grdis[cc++] = ele->node[19]->sol.a.da[1][0];
   grdis[cc++] = ele->node[19]->sol.a.da[1][1];
   grdis[cc++] = ele->node[19]->sol.a.da[1][2];
   grdis[cc++] = ele->node[12]->sol.a.da[1][0];
   grdis[cc++] = ele->node[12]->sol.a.da[1][1];
   grdis[cc++] = ele->node[12]->sol.a.da[1][2];
   grdis[cc++] = ele->node[13]->sol.a.da[1][0];
   grdis[cc++] = ele->node[13]->sol.a.da[1][1];
   grdis[cc++] = ele->node[13]->sol.a.da[1][2];
   grdis[cc++] = ele->node[14]->sol.a.da[1][0];
   grdis[cc++] = ele->node[14]->sol.a.da[1][1];
   grdis[cc++] = ele->node[14]->sol.a.da[1][2];
   grdis[cc++] = ele->node[15]->sol.a.da[1][0];
   grdis[cc++] = ele->node[15]->sol.a.da[1][1];
   grdis[cc++] = ele->node[15]->sol.a.da[1][2];
   for (i=0; i<80; i++) dfie[i] = 0.0;
   }
   /*------------------------------------------------*/
    if (init==7)
    {/*init==7*/
      for (ieig=0;ieig<opt->oeig->numeigv;ieig++)
      {/* loop eigenvalues */
         cc=0;
         reig[ieig][cc++] = ele->node[0]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[0]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[0]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[1]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[1]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[1]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[2]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[2]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[2]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[3]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[3]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[3]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[4]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[4]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[4]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[5]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[5]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[5]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[6]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[6]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[6]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[7]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[7]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[7]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[8]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[8]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[8]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[9]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[9]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[9]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[10]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[10]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[10]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[11]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[11]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[11]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[16]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[16]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[16]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[17]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[17]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[17]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[18]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[18]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[18]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[19]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[19]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[19]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[12]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[12]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[12]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[13]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[13]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[13]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[14]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[14]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[14]->sol.a.da[ieig][2];
         reig[ieig][cc++] = ele->node[15]->sol.a.da[ieig][0];
         reig[ieig][cc++] = ele->node[15]->sol.a.da[ieig][1];
         reig[ieig][cc++] = ele->node[15]->sol.a.da[ieig][2];
      }/* loop eigenvalues */
    }/*init==7*/
   /*------------------------------------------------*/
  }
/*-------------------------------------  type of element formulation ---*/
  iform   = ele->e.c1->form;/*=1:linear:=2 total lagrangian formulation */

  ste  = 0.;
  stv  = 0.;
  dste = 0.;
/*------------------------------------ check calculation of mass matrix */
if (init==7)
{
  krexpo = opt->oeig->rhoks;                 /* exponent in Kreiselmeier*/
  /*---------------------------------------------------- get density ---*/
  #ifdef D_OPTIM                   /* include optimization code to ccarat */
   if(ele->e.c1->elewa->matdata==NULL) c1_getdensity(mat, &density);
   else density = ele->e.c1[0].elewa[0].matdata[0];
  #else
  c1_getdensity(mat, &density);
  #endif /* stop including optimization code to ccarat :*/
  for (i=0; i<60; i++) lmvec[i] = 0.0;
  amzero(&estiflo_a);
}
/*-------------------------------------------  initialize total mass ---*/
  totmas = 0.0;
  emasdg = 0.0;
/*================================================ integration loops ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgrr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
    /*============================== gaussian point and weight at it ===*/
    e2   = data->xgss[ls];
    facs = data->wgts[ls];
    for (lt=0; lt<nit; lt++)
    {
      ip++;
      /*============================ gaussian point and weight at it ===*/
      e3   = data->xgtt[lt];
      fact = data->wgtt[lt];
      /*------------------------- shape functions and their derivatives */
      c1_funct_deriv(funct,deriv,e1,e2,e3,iel,1);
      /*------------------------------------ compute jacobian matrix ---*/
      c1_jaco (deriv,xjm,&det,xyze,iel);
      fac = facr * facs *  fact * det;
      /*------------------------------------------- calculate volume ---*/
      if (init==5||init==3)
      {
          stv += fac;
          continue;
      }
      amzero(&bop_a);
      /*-- local element coordinate system for anisotropic materials ---*/
      c1tram (xjm,g,gi);
      /*--------------------------------------- calculate operator B ---*/
      c1_bop(bop,bn,deriv,xjm,det,iel);
      /*--------------------------- compute displacement derivatives ---*/
      c1_disd (bop,edis,disd,iel) ;
      /*---------------- include initial displacements to b-operator ---*/
      if(iform==2 && mat->mattyp!=m_pl_mises_ls)
      {
        c1_bdis (bop,disd,iel) ;
      }
      /*------------------------------- get actual strains -> strain ---*/
      c1_eps (disd,strain,iform);
      /*------------------------------------------ call material law ---*/
      c1_call_mat(ele, mat,ip,F,strain,D,disd,g,gi,istore,newval);
      /*------------------------------------ calculate strain energy ---*/
      if (init==2)
      {
          if(mat->mattyp==m_nhmfcc)
          {
            ste += strain[0]*fac;
            continue;
          }


          ste += (F[0]*strain[0] + F[1]*strain[1] + F[2]*strain[2] +
                  F[3]*strain[3] + F[4]*strain[4] + F[5]*strain[5])
                  *fac*0.5;
          continue;

      }
      /*----------- calculate derivatives material matrix and stress ---*/
      if(init==4)
      {
        for (i=0;i<6;i++) DF[i] = 0.;
        if(mat->mattyp==m_nhmfcc)
        {
          c1_call_matd(ele, mat,DF,strain,DD,g);
          dste -= strain[0]*fac;


        }
        else
        {
          c1_call_matd(ele, mat,DF,strain,DD,g);
          dste -= (DF[0]*strain[0] + DF[1]*strain[1] + DF[2]*strain[2] +
                   DF[3]*strain[3] + DF[4]*strain[4] + DF[5]*strain[5])
                   *fac*0.5;
        }
      }
      /*------------- calculate derivatives for selfadjoint problems ---*/
      if(init==6)
      {
        for (i=0;i<6;i++) DF[i] = 0.;
        c1_call_matd(ele, mat,DF,strain,DD,g);
        /* derivative of internal force vector */
        c1dfi (DF,fac,bop,disd,nd,dfie);
      }
      /*--------------- calculate derivatives for frequency problems ---*/
      if(init==7)
      {
       /*- eval. derivatives of stiffness and mass matrix(with dens=1.0)*/
        for (i=0;i<6;i++) DF[i] = 0.;
        c1_call_matd(ele, mat,DF,strain,DD,g);
        c1_keku(estiflo,bop,DD,fac,nd,numeps);

        density = 1.0;
        facm = fac * density;
        totmas += facm;
        c1cptp (funct,lmvec,NULL,&emasdg,iel,1,facm);
      }
  }/*============================================== end of loop over lt */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
  /*---------------------------------------------------- mass matrix ---*/
  if (init==7)
  {
  /*--------------------------------------------------------------------*/
    fac=3.0*totmas/emasdg;
    for (i=0; i<nd; i++) lmvec[i] *= fac;

    for (ieig=0;ieig<opt->oeig->numeigv;ieig++)
    {/* loop eigenvalues */
      /*-------------------------------------------- k,s - w|2 * m,s ---*/
      for (i=0; i<(iel*3); i++)
      {
        estiflo[i][i] -= opt->oeig->eigv[ieig] * lmvec[i];
      }
      /*------------------------- w|2,s = rt * (k,s - w|2 * m,s) * r ---*/
      w2s1[ieig]=0.0;
      for (i=0; i<(iel*3); i++) hvar1[i]=0.0;
      for (i=0; i<(iel*3); i++)
      {
        for (j=0; j<(iel*3); j++) hvar1[i] += reig[ieig][j]*estiflo[j][i];
	w2s1[ieig] += hvar1[i]*reig[ieig][i];
      }
      /*---------------------------------------------------- scaling ---*/
      w2s1[ieig]=w2s1[ieig]/opt->oeig->eigs[ieig];
      /*-------------------------------- teta_s = w|2,s / (8*PI*fie) ---*/
      teta_s1[ieig] =  w2s1[ieig] / (4.0*PI*sqrt(opt->oeig->eigv[ieig]));
      teta1[ieig]   =  sqrt(opt->oeig->eigv[ieig]) / (2*PI);
    }
      /*------------------------- Kreisselmeier-Steinhaeuser (MAUTE) ---*/
      stf = 0.0;
      if(opt->oeig->numeigv==1)
      {
        stf = - teta_s1[0];
      }
      else
      {
        kro = 0.0;
        kru = 0.0;
        for (j=0;j<opt->oeig->numeigv;j++) kro -= exp(-krexpo*teta1[j])*teta_s1[j];
        for (j=0;j<opt->oeig->numeigv;j++) kru += exp(-krexpo*teta1[j])           ;
        stf = kro / kru;
      }
  /*--------------------------------------------------------------------*/
  }
/*----------------------------------------------------------------------*/
  switch(init)
  {
  case 5:
          (*retval) = stv;
  break;
  case 2:
          (*retval) = ste;
  break;
  case 3:
          if(ele->e.c1->elewa->matdata==NULL) c1_getdensity(mat,&dens);
          else
          {
            dens = ele->e.c1->elewa->matdata[0];
          }
          stm       = stv*dens;
          (*retval) = stm;
  break;
  case 4:
          (*retval) = dste;
  break;
  case 6:
          rs = 0.;
          for (i=0;i<nd;i++) rs += grdis[i] * dfie[i];
          (*retval) = -1. * rs;
  break;
  case 7:
          (*retval) = stf;
  break;
  }

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_oint */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | derivative of internal force vector                   al    9/01     |
 *----------------------------------------------------------------------*/
void c1dfi(DOUBLE  *F,  /*  force vector integral (stress-resultants)  */
           DOUBLE   fac, /*  multiplier for numerical integration       */
           DOUBLE **bop, /*  b-operator matrix                          */
           DOUBLE *disd, /* displacement derivatives               */
           INT      nd,  /*  total number degrees of freedom of element */
           DOUBLE  *fie) /*  internal force vector                      */
{
/*----------------------------------------------------------------------*/
INT i,j,k;
DOUBLE n11,n22,n33,n12,n23,n31;
DOUBLE dd11,dd22,dd33,dd12,dd21,dd13,dd23,dd31,dd32;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1dfi");
#endif
/*---------------------------- set values of force vector components ---*/
  n11 = F[0]*fac;
  n22 = F[1]*fac;
  n33 = F[2]*fac;
  n12 = F[3]*fac;
  n23 = F[4]*fac;
  n31 = F[5]*fac;

  dd11 = disd[0] + 1.;
  dd22 = disd[1] + 1.;
  dd33 = disd[2] + 1.;
  dd12 = disd[3];
  dd21 = disd[4];
  dd13 = disd[5];
  dd31 = disd[6];
  dd23 = disd[7];
  dd32 = disd[8];
/*----------------------------- updated lagrange or geometric linear ---*/
  for (j=2; j<nd; j+=3)
  {
    k=j-1;
    i=j-2;
    /*
    fie[i]+=  bop[0][i]*dd11*n11 + bop[1][i]*dd12*n22 + bop[2][i]*dd13*n33 +
              bop[3][i]*(dd12+dd11)*n12 +
              bop[4][i]*(dd13+dd11)*n23 +
              bop[5][i]*(dd13+dd12)*n31;
    fie[k]+=  bop[0][k]*dd21*n11 + bop[1][k]*dd22*n22 + bop[2][k]*dd23*n33 +
              bop[3][k]*(dd22+dd21)*n12 +
              bop[4][k]*(dd23+dd21)*n23 +
              bop[5][k]*(dd23+dd22)*n31;
    fie[j]+=  bop[0][j]*dd31*n11 + bop[1][j]*dd32*n22 + bop[2][j]*dd33*n33 +
              bop[3][j]*(dd32+dd31)*n12 +
              bop[4][j]*(dd33+dd31)*n23 +
              bop[5][j]*(dd33+dd32)*n31;
    /**/
    fie[i]+=  bop[0][i]*n11 + bop[1][i]*n22 + bop[2][i]*n33 +
              bop[3][i]*n12 +
              bop[4][i]*n23 +
              bop[5][i]*n31;
    fie[k]+=  bop[0][k]*n11 + bop[1][k]*n22 + bop[2][k]*n33 +
              bop[3][k]*n12 +
              bop[4][k]*n23 +
              bop[5][k]*n31;
    fie[j]+=  bop[0][j]*n11 + bop[1][j]*n22 + bop[2][j]*n33 +
              bop[3][j]*n12 +
              bop[4][j]*n23 +
              bop[5][j]*n31;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1dfi */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
