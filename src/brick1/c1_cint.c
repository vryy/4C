/*!----------------------------------------------------------------------
\file
\brief contains the control program 'c1_cint' for integration over
       the element volume for a 3D hex element

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*!
\addtogroup BRICK1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief integration routine for BRICK1 element

<pre>                                                              al 06/02
This routine performs integration of an 3D-hex-element.

</pre>
\param           ele ELEMENT*  (i)   element data
\param          data C1_DATA*  (i)   hex element data
\param          mat MATERIAL*  (i)   material data
\param    estif_global ARRAY*  (o)   element stiffness matrix
\param          force DOUBLE*  (o)   vector for internal forces
\param             init  INT*  (i)   flag for initialization (alloc mem...)


 *----------------------------------------------------------------------*
 |                                                                      |
 |    6-------18-------2           6----------------2                   |
 |    |\               |\          |\               |\                  |
 |    | \              | \         | \        S     | \                 |
 |    |  13            |  9        |  \       |     |  \                |
 |    |   \            |   \       |   \      |     |   \               |
 |   14    \           10   \      |    \  \  |     10   \              |
 |    |     5-------17-------1     |     5----------------1             |
 |    |     |          |     |     |     |   \|     |     |             |
 |    |     |          |     |     | T---|----o--------   |             |
 |    |     |          |     |     |     |    |\    |     |             |
 |    7-----|-19-------3     |     7-----|----|-\---3     |             |
 |     \    12          \    8      \    |    |  \   \    |             |
 |      \   |            \   |       \   |    |   R   \   |             |
 |       15 |             11 |        \  |    |        \  |             |
 |        \ |              \ |         \ |              \ |             |
 |         \|               \|          \|               \|             |
 |          4-------16-------0           4----------------0             |
 |                                                                      |
 *----------------------------------------------------------------------*


\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_main()

*----------------------------------------------------------------------*/
void c1_cint(
             ELEMENT   *ele,
             C1_DATA   *data,
             MATERIAL  *mat,
             ARRAY     *estif_global,
             ARRAY     *emass_global,
             DOUBLE    *force,
             INT        init
             )
{
INT                 i,j,k;            /* some loopers */
INT                 nir,nis,nit;      /* num GP in r/s/t direction */
INT                 lr, ls, lt;       /* loopers over GP */
INT                 ip;
INT                 iel;              /* numnp to this element */
INT                 nd, nd1;
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 0;/* controls evaluation of new stresses    */
const INT           numdf =3;
const INT           numeps=6;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE disd[9];
DOUBLE F[6]; /* element stress vector   (stress-resultants) */
DOUBLE fielo[81];
DOUBLE strain[6];
DOUBLE xyze[60];
DOUBLE edis[60];
DOUBLE  g[6][6]; /* transformation matrix s(glob)= g*s(loc)   */
DOUBLE gi[6][6]; /* inverse of g          s(loc) = gi*s(glob) */

/*-------------------------  for postprocessing - stress calculation ---*/
DOUBLE gpstrs[27][26]; /* [number of gp   ][stresses] */
DOUBLE nostrs[20][26]; /* [number of nodes][stresses] */
DOUBLE srst[6];
DOUBLE s123[12];
DOUBLE gpcod[3]; /* natural coordinates of g.p.*/
/*-------------------------------------------  for eas elements only ---*/
INT    l1, l3, ihyb, cc;
DOUBLE ehdis[3][10],fi[6][6],ff[6][6];
DOUBLE det0, det1;
DOUBLE bn1[3][10];
DOUBLE disd1[9];
DOUBLE fieh[30];
DOUBLE epsh[6];

static DOUBLE **estiflo;
static ARRAY    estiflo_a; /* local element stiffness matrix ke for eas */

static DOUBLE **estif9;
static ARRAY    estif9_a;   /* element stiffness matrix ke for eas */

/*----------------------------------------------------------------------*/
static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
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
static DOUBLE **estif;    /* element stiffness matrix ke */
static DOUBLE **emass;     /* element mass matrix */
DOUBLE          lmvec[ 60]; /* lumped     mass vector */
DOUBLE          consm[400]; /* consistent mass matrix */
DOUBLE emasdg;  /* factor for lumped mass matrix */
DOUBLE              density, facm, totmas;
DOUBLE det;
/*----------------------------------------------------------------------*/
INT    iform;             /* index for nonlinear formulation of element */
INT    calstr;            /* flag for stress calculation                */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1_cint");
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
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");

  bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
  bn        = amdef("bnop" ,&bnop_a,3     ,       MAXNOD_BRICK1 ,"DA");
  estif9    = amdef("estif9"  ,&estif9_a ,54,54,"DA");

  estiflo   = amdef("estiflo"  ,&estiflo_a ,60,60,"DA");
goto end;
}
else if(init==2)
{
  istore = 1;
}
else if(init==3)
{
  if(ele->e.c1->stresstyp != c1_nostr) calstr = 1;
  newval = 0;
}
else if(init==0)
{
  newval = 1;
}
/*------------------------------------------- integration parameters ---*/
c1intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
amzero(&estiflo_a);


for (i=0; i<81; i++) fielo[i] = 0.0;
/*------------------------------------ check calculation of mass matrix */
if (init==4 || init==5)
{
  /*---------------------------------------------------- get density ---*/
  #ifdef D_OPTIM                   /* include optimization code to ccarat */
   if(ele->e.c1->elewa->matdata==NULL) c1_getdensity(mat, &density);
   else density = ele->e.c1[0].elewa[0].matdata[0];
  #else
  c1_getdensity(mat, &density);
  #endif /* stop including optimization code to ccarat :*/
  for (i=0; i<60; i++)  lmvec[i] = 0.0;
  for (i=0; i<400; i++) consm[i] = 0.0;
  if(emass_global!=NULL) amzero(emass_global);
  emass     = emass_global->a.da;
}
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
  }
/*-------------------------------------  type of element formulation ---*/
  iform   = ele->e.c1->form;/*=1:linear:=2 total lagrangian formulation */
/*-------------------------------------------  for eas elements only ---*/
  ihyb = ele->e.c1->nhyb;
  if(ihyb>0)
  {
    l1=3*ihyb;
    l3=ihyb;
    for (i=0; i<30; i++) fieh[i] =0.;
    amzero(&estif9_a);
    /*--------------------------------- update of strain paramenters ---*/
    c1upenh(ele, edis, ehdis, l1, l3);
    /*----- determine jacobi matrix at center of element (for eas) ---*/
    c1_funct_deriv(funct,deriv,0.,0.,0.,iel,1);
    c1_jaco (deriv,xjm,&det0,xyze,iel);
    c1t0 (fi,ff,xjm);
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
      amzero(&bop_a);
      /*------------------------------------------------ mass matrix ---*/
      if (init==4)
      {
        facm = fac * density;
        totmas += facm;
        c1cptp (funct,lmvec,consm,&emasdg,iel,1,facm);
      }
      /*-- local element coordinate system for anisotropic materials ---*/
      c1tram (xjm,g,gi);
      /*----------------------------------  eas element (small def.) ---*/
      if(ihyb>0)
      {
        det1 = det;
        c1bop9 (bop, bn1,fi, disd1,ehdis, det0,det1, e1,e2,e3, iel, l1,l3);
      }
      /*--------------------------------------- calculate operator B ---*/
      c1_bop(bop,bn,deriv,xjm,det,iel);
      /*--------------------------- compute displacement derivatives ---*/
      c1_disd (bop,edis,disd,iel) ;
      /*---------------- include initial displacements to b-operator ---*/
      if(iform==2 && mat->mattyp!=m_pl_mises_ls)
      {
        c1_bdis (bop,disd,iel) ;
       if(ihyb >0) c1bdish (bop,bn1,disd1,iel,l3) ;
      }
      /*------------------------------- get actual strains -> strain ---*/
      c1_eps (disd,strain,iform);
      /*--------------------------------------------------- eas part ---*/
      if(ihyb>0)
      {
        c1_eps (disd1,epsh,iform);
        cc=0;
        for (i=0; i<ihyb; i++) {
          for (j=0; j<6; j++) {
            for (k=0; k<3; k++) {
              strain[j] +=  bop[j][k + (i+iel)*3] * ehdis[k][i];}}}
      }
      /*------------------------------------------ call material law ---*/
      c1_call_mat(ele, mat,ip,F,strain,D,disd,g,gi,istore,newval);
      /*----------- calculate element stresses at integration points ---*/
      if(calstr==1)
      {
        for (i=0; i<6; i++) srst[i]=F[i]; /* stresses at gauss point    */
        for (i=0; i<6; i++) s123[i]=F[i]; /* princ. stress,  directions */
        c1trss2local(srst, gi);

        c1_cstr (srst ,s123, gpstrs[ip]); /* gpstrs[ip][0..24] */
        /*------------- global coordinates of actual  gaussian point ---*/
        c1gcor  (funct,xyze, iel  ,gpcod);
        /*--------------------------------------------- store values ---*/
        for (i=0; i<25; i++)
        {
          ele->e.c1->stress_GP.a.d3[0][i][ip]= gpstrs[ip][i];
        }
        for (i=0; i<3; i++)
        {
          ele->e.c1->stress_GP.a.d3[0][i+24][ip]= gpcod[i];
        }
        /*--------------------------------------------------------------*/
      }
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*-------------------------------- elastic stiffness matrix ke ---*/
        if(ihyb>0)
        {
          nd1 = 3*iel + l1;
          c1_keku(estif9,bop,D,fac,nd1,numeps);
        }
        else
        {
          c1_keku(estiflo,bop,D,fac,nd,numeps);
        }
      /*------------------------------ geometric stiffness-matrix kg ---*/
        /* besser eigener kenner fuer geom. n.l.! hier T.L.Abfrage!*/
        if(iform==2 && mat->mattyp!=m_pl_mises_ls) c1vkg (estiflo,bn,F,fac,iel);
      /*--------------- nodal forces fi from integration of stresses ---*/
        if (force)
        {
           c1fi (F,fac,bop,nd,fielo);
      /*------------------------ compute residuals----------------------*/
          if(ihyb>0)
          {
            c1res (F,fac,bop,iel,fieh,l3);
          }
        }
      }
  }/*============================================== end of loop over lt */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
  if(calstr==1)
  {
    /*---- extrapolation of stress from gauss points to nodal points ---*/
    c1_sext(nostrs, funct, deriv,xjm,xyze,
            gpstrs,
            data->xgrr, data->xgss, data->xgtt, nir,nis,nit,iel);
    /*------------------------------------------------- store values ---*/
    for (i=0; i<25; i++) /* number of stress components */
    {
      for (j=0; j<iel; j++) /* number of nodes */
      {
        ele->e.c1->stress_ND.a.d3[0][i][j] = nostrs[j][i];
      }
    }
    goto end;
  }
/*----------------------------------------------------------------------*/
  if(ihyb>0)
  {
    c1rkefi(ele, estif9, estiflo, fieh, fielo, l1);
  }
/*----------------------------------------------------------------------*/
  /*---------------------------------------------------- mass matrix ---*/
  if (init==4)
  {
    fac=3.0*totmas/emasdg;
    for (i=0; i<nd; i++) fielo[i] = lmvec[i]*fac;
  }
  if (init==5)
  {
    cc=0;
    for (i=0; i<nd; i++) for (j=0; j<nd; j++) emass[i][j] = consm[cc++];
  }
/*----------------------------------------------------------------------*/
/* reorder stiffness and element forces for 'gid' element-node topo...  */
    if(iel==20)
    {
      c1kgid(estiflo, estif);
      if (force) c1fgid(fielo,force);
    }
    else
    {
      for (i=0; i<24; i++){
      for (j=0; j<24; j++){ estif[i][j] = estiflo[i][j];}}
      if (force) for (i=0; i<24; i++) force[i] = fielo[i];
    }
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_cint */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief evaluates element forces

<pre>                                                              al 06/02
This routine evaluates element forces of an 3D-hex-element.

</pre>
\param    F   DOUBLE*   (i)   force vector integral (stress-resultants)
\param  fac   DOUBLE    (i)   multiplier for numerical integration
\param  bop   DOUBLE**  (i)   b-operator matrix
\param   nd   INT       (i)   total number degrees of freedom of element
\param  fie   DOUBLE*   (o)   internal force vector

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1fi( DOUBLE  *F,   /*  force vector integral (stress-resultants)  */
           DOUBLE   fac, /*  multiplier for numerical integration       */
           DOUBLE **bop, /*  b-operator matrix                          */
           INT      nd,  /*  total number degrees of freedom of element */
           DOUBLE  *fie) /*  internal force vector                      */
{
/*----------------------------------------------------------------------*/
INT i,j,k;
DOUBLE n11,n22,n33,n12,n23,n31;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1fi");
#endif
/*---------------------------- set values of force vector components ---*/
  n11 = F[0]*fac;
  n22 = F[1]*fac;
  n33 = F[2]*fac;
  n12 = F[3]*fac;
  n23 = F[4]*fac;
  n31 = F[5]*fac;
/*----------------------------- updated lagrange or geometric linear ---*/
  for (j=2; j<nd; j+=3)
  {
    k=j-1;
    i=j-2;
    fie[i]+=  bop[0][i]*n11 + bop[1][i]*n22 + bop[2][i]*n33 +
              bop[3][i]*n12 + bop[4][i]*n23 + bop[5][i]*n31;
    fie[k]+=  bop[0][k]*n11 + bop[1][k]*n22 + bop[2][k]*n33 +
              bop[3][k]*n12 + bop[4][k]*n23 + bop[5][k]*n31;
    fie[j]+=  bop[0][j]*n11 + bop[1][j]*n22 + bop[2][j]*n33 +
              bop[3][j]*n12 + bop[4][j]*n23 + bop[5][j]*n31;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1fi */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief evaluates material transformation matricies

<pre>                                                              al 06/02
This routine evaluates material transformation matricies for a 3D-hex-element.

</pre>
\param       xjm   DOUBLE**  (i)  jacobian matrix r,s,t-direction
\param   g[6][6]   DOUBLE    (o)  transformation matrix s(glob)=g*s(loc)
\param  gi[6][6]   DOUBLE    (o)  inverse of g          s(loc) =gi*s(glob)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1tram( DOUBLE **xjm,   /* jacobian matrix r,s,t-direction         */
             DOUBLE g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
             DOUBLE gi[6][6])/* inverse of g          s(loc) =gi*s(glob)*/
{
/*----------------------------------------------------------------------*/
DOUBLE dm;
DOUBLE x1r, x2r, x3r, x1s, x2s, x3s;
DOUBLE x1n, x2n, x3n, x1c, x2c, x3c;
DOUBLE a11, a21, a31, a12, a22, a32, a13, a23, a33;
DOUBLE dum[3][3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1tram");
#endif
/*----------------------------------------------------------------------*
 | local axes: x = tangent at r-axes                                    |
 |             y = cross product z*x                                    |
 |             z = normal at surface (t-axes)                           |
 |                                                                      |
 *----------------------------------------------------------------------*/
      x1r=xjm[0][0];
      x2r=xjm[0][1];
      x3r=xjm[0][2];
      x1s=xjm[1][0];
      x2s=xjm[1][1];
      x3s=xjm[1][2];

      x1n=x3s*x2r - x2s*x3r;
      x2n=x1s*x3r - x3s*x1r;
      x3n=x2s*x1r - x1s*x2r;

      x1c=x3r*x2n - x2r*x3n;
      x2c=x1r*x3n - x3r*x1n;
      x3c=x2r*x1n - x1r*x2n;
/*------------- reduce to unit vectors (divide by the scalar lenght) ---*/
      dm = sqrt (x1r*x1r + x2r*x2r + x3r*x3r);
      dm = 1./dm;
      x1r = x1r*dm;
      x2r = x2r*dm;
      x3r = x3r*dm;
      dm = sqrt (x1c*x1c + x2c*x2c + x3c*x3c);
      dm = 1./dm;
      x1c = x1c*dm;
      x2c = x2c*dm;
      x3c = x3c*dm;
      dm = sqrt (x1n*x1n + x2n*x2n + x3n*x3n);
      dm = 1./dm;
      x1n = x1n*dm;
      x2n = x2n*dm;
      x3n = x3n*dm;
/*---------------------------------- set matrix of direction cosines ---*/
      dum[0][0]=x1r;
      dum[0][1]=x2r;
      dum[0][2]=x3r;
      dum[1][0]=x1c;
      dum[1][1]=x2c;
      dum[1][2]=x3c;
      dum[2][0]=x1n;
      dum[2][1]=x2n;
      dum[2][2]=x3n;
/*---------------- evaluate final transformation matrix g and g(inv) ---*/
      a11=dum[0][0];
      a21=dum[0][1];
      a31=dum[0][2];
      a12=dum[1][0];
      a22=dum[1][1];
      a32=dum[1][2];
      a13=dum[2][0];
      a23=dum[2][1];
      a33=dum[2][2];


      gi[0][0] = a11*a11;
      gi[1][0] = a12*a12;
      gi[2][0] = a13*a13;
      gi[3][0] = a11*a12;
      gi[4][0] = a12*a13;
      gi[5][0] = a11*a13;

      gi[0][1] = a21*a21;
      gi[1][1] = a22*a22;
      gi[2][1] = a23*a23;
      gi[3][1] = a21*a22;
      gi[4][1] = a22*a23;
      gi[5][1] = a21*a23;

      gi[0][2] = a31*a31;
      gi[1][2] = a32*a32;
      gi[2][2] = a33*a33;
      gi[3][2] = a31*a32;
      gi[4][2] = a32*a33;
      gi[5][2] = a31*a33;

      gi[0][3] = a11*a21 * 2.     ;
      gi[1][3] = a12*a22 * 2.     ;
      gi[2][3] = a13*a23 * 2.     ;
      gi[3][3] = a11*a22 + a21*a12;
      gi[4][3] = a12*a23 + a22*a13;
      gi[5][3] = a11*a23 + a21*a13;

      gi[0][4] = a21*a31 * 2.     ;
      gi[1][4] = a22*a32 * 2.     ;
      gi[2][4] = a23*a33 * 2.     ;
      gi[3][4] = a21*a32 + a31*a22;
      gi[4][4] = a22*a33 + a32*a23;
      gi[5][4] = a21*a33 + a31*a23;

      gi[0][5] = a11*a31 * 2.     ;
      gi[1][5] = a12*a32 * 2.     ;
      gi[2][5] = a13*a33 * 2.     ;
      gi[3][5] = a11*a32 + a31*a12;
      gi[4][5] = a12*a33 + a32*a13;
      gi[5][5] = a11*a33 + a31*a13;

      a11=dum[0][0];
      a12=dum[0][1];
      a13=dum[0][2];
      a21=dum[1][0];
      a22=dum[1][1];
      a23=dum[1][2];
      a31=dum[2][0];
      a32=dum[2][1];
      a33=dum[2][2];

      g[0][0] = a11*a11;
      g[1][0] = a12*a12;
      g[2][0] = a13*a13;
      g[3][0] = a11*a12;
      g[4][0] = a12*a13;
      g[5][0] = a11*a13;

      g[0][1] = a21*a21;
      g[1][1] = a22*a22;
      g[2][1] = a23*a23;
      g[3][1] = a21*a22;
      g[4][1] = a22*a23;
      g[5][1] = a21*a23;

      g[0][2] = a31*a31;
      g[1][2] = a32*a32;
      g[2][2] = a33*a33;
      g[3][2] = a31*a32;
      g[4][2] = a32*a33;
      g[5][2] = a31*a33;

      g[0][3] = a11*a21 * 2.     ;
      g[1][3] = a12*a22 * 2.     ;
      g[2][3] = a13*a23 * 2.     ;
      g[3][3] = a11*a22 + a21*a12;
      g[4][3] = a12*a23 + a22*a13;
      g[5][3] = a11*a23 + a21*a13;

      g[0][4] = a21*a31 * 2.     ;
      g[1][4] = a22*a32 * 2.     ;
      g[2][4] = a23*a33 * 2.     ;
      g[3][4] = a21*a32 + a31*a22;
      g[4][4] = a22*a33 + a32*a23;
      g[5][4] = a21*a33 + a31*a23;

      g[0][5] = a11*a31 * 2.     ;
      g[1][5] = a12*a32 * 2.     ;
      g[2][5] = a13*a33 * 2.     ;
      g[3][5] = a11*a32 + a31*a12;
      g[4][5] = a12*a33 + a32*a13;
      g[5][5] = a11*a33 + a31*a13;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1tram */

/*!----------------------------------------------------------------------
\brief transformation of local material-matrix to global axes

<pre>                                                              al 06/02
This routine transforms of local material-matrix to global axes for a 3D-hex-element.

</pre>
\param       **d   DOUBLE  (o)  material matrix
\param   g[6][6]   DOUBLE  (i)  transformation matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1gld(DOUBLE **d,     /* material matrix                           */
           DOUBLE g[6][6]) /* transformation matrix                     */
{
/*----------------------------------------------------------------------*/
INT i,j,k;
DOUBLE dgt[6][6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1gld");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++) for (j=0; j<6; j++) dgt[i][j] = 0.;
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      for (k=0; k<6; k++)
      {
        dgt[i][j] += d[i][k]*g[j][k];
      }
    }
  }
/* R(I,J) = A(I,K)*B(K,J) ---  R = A*B */
  for (i=0; i<6; i++) for (j=0; j<6; j++) d[i][j] = 0.;
  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      for (k=0; k<6; k++)
      {
        d[i][j] += g[i][k]*dgt[k][j];
      }
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1gld */

/*!----------------------------------------------------------------------
\brief transformation of global stress vector to local axes

<pre>                                                              al 06/02
This routine transforms of global stress vector to local axes for a 3D-hex-element.

</pre>
\param         s   DOUBLE* (o)  stress vector to be transformed
\param   g[6][6]   DOUBLE  (i)  inverse of transformation matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1trss2local(DOUBLE *s,       /* stress vector to be transformed   */
                  DOUBLE gi[6][6]) /* inverse of transformation matrix  */
{
/*----------------------------------------------------------------------*/
INT i,j;
DOUBLE sh[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1trss2local");
#endif
/*--------------- sig(global) --> sig(local)  s(loc) = g(inv)*s(glo) ---*/
  for (i=0; i<6; i++) sh[i] = s[i];
  for (i=0; i<6; i++)
  {
    s[i] = 0.;
    for (j=0; j<6; j++)
    {
      s[i] += gi[i][j]*sh[j];
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1trss2local */

/*!----------------------------------------------------------------------
\brief transformation of local stress vector to global axes

<pre>                                                              al 06/02
This routine transforms of local stress vector to global axes for a 3D-hex-element.

</pre>
\param         s   DOUBLE* (o)  stress vector to be transformed
\param   g[6][6]   DOUBLE  (i)  transformation matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1trss2global(DOUBLE *s,       /* stress vector to be transformed  */
                   DOUBLE g[6][6])  /* transformation matrix            */
{
/*----------------------------------------------------------------------*/
INT i,j;
DOUBLE sh[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1trss2global");
#endif
/*--------------- sig(local) --> sig(global)  s(glo) = g(inv)*s(loc) ---*/
  for (i=0; i<6; i++) sh[i] = s[i];
  for (i=0; i<6; i++)
  {
    s[i] = 0.;
    for (j=0; j<6; j++)
    {
      s[i] += g[i][j]*sh[j];
    }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1trss2global */
/*----------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief consistent mass matrix for brick1

<pre>                                                              al 06/02
This routine evaluates mass matrix for a 3D-hex-element.

</pre>
\param  *funct   DOUBLE  (i)   shape functions
\param  *emass   DOUBLE  (o)   mass vector
\param  *emasdg  DOUBLE  (o)   factor for lumped mass
\param  iel      INT     (i)   number of nodes
\param  ilmp     INT     (i)   flag for lumped/consistent matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1cptp(
            DOUBLE     *funct,
            DOUBLE     *lmass,
            DOUBLE     *consm,
            DOUBLE     *emasdg,
            INT         iel,
            INT         ilmp,
            DOUBLE      fac
            )
{
/*----------------------------------------------------------------------*/
INT i,j,k,l,cc;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("c1cptp");
#endif
/*----------------------------------------------------------------------*/
switch (ilmp)
{
/*------------------------------------------- consistent mass (full) ---*/
case 2:
  cc = 0;
  for (i=0; i<iel; i++)
  {
    for (k=0; k<3; k++)
    {
      for (j=0; j<iel; j++)
      {
        for (l=0; l<3; l++)
        {
           if(k==l) consm[cc] += fac*funct[i]*funct[j];
           cc++;
        }
      }
    }
  }
break;/*----------------------------------------------------------------*/
/*------------------------------------------------------ lumped mass ---*/
case 1:
  cc = 0;
  for (i=0; i<iel; i++)
  {
    for (l=0; l<3; l++)
    {
            (*emasdg) += fac*funct[i]*funct[i];
            lmass[cc] += fac*funct[i]*funct[i];
            cc++;
    }
  }
break;/*----------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
default:
   dserror("unknown type of mass matrix for hex element");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1cptp */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
