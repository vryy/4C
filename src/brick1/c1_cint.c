/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_disd' which calclate displacement
       derivatives for a 3D hex element

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
\param           *ele ELEMENT  (i)   element data
\param          *data C1_DATA  (i)   hex element data
\param           *mat MATERIAL (i)   material data
\param  *estif_global ARRAY    (o)   element stiffness matrix
\param         *force DOUBLE   (o)   vector for internal forces
\param         *init  INT      (i)   flag for initialization (alloc mem...)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_cint(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             MATERIAL  *mat,
             ARRAY     *estif_global, 
             double    *force,  /* global vector for internal forces (initialized!) */
             int        init
             )
{
int                 i,j,k;            /* some loopers */
int                 nir,nis,nit;      /* num GP in r/s/t direction */
int                 lr, ls, lt;       /* loopers over GP */
int                 ip;
int                 iel;              /* numnp to this element */
int                 nd, nd1;
int                 dof;
int                 istore = 0;/* controls storing of new stresses to wa */
int                 newval = 0;/* controls evaluation of new stresses    */
const int           numdf =3;
const int           numeps=6;

double              fac;
double              e1,e2,e3;         /*GP-coords*/
double              facr,facs,fact;   /* weights at GP */
double              xnu;              /* value of shell shifter */
double              weight;
double disd[9];
double F[6]; /* element stress vector   (stress-resultants) */
double fie[81];
double fielo[81];
double strain[6];
double xyze[60];
double edis[60];  
double  g[6][6]; /* transformation matrix s(glob)= g*s(loc)   */
double gi[6][6]; /* inverse of g          s(loc) = gi*s(glob) */

/*-------------------------  for postprocessing - stress calculation ---*/
/* 
gpstrs[0][ 0.. 7] [stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz ]
gpstrs[0][ 6..11] [stress-rr stress-ss stress-tt stress-rs stress-st stress-tr ]
gpstrs[0][12..14] [stress-11 stress-22 stress-33 ]
gpstrs[0][15..23] [ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3 ]
*/
double gpstrs[27][26];
/*nostrs[nn] [ 0.. 5]: x-coord.    y-coord.    z-coord.     stress-rr   stress-ss   stress-tt   stress-rs   stress-st   stress-tr*/  
/*nostrs[nn] [ 6..11]: x-coord.    y-coord.    z-coord.     stress-xx   stress-yy   stress-zz   stress-xy   stress-yz   stress-xz*/ 
/*nostrs[nn] [12..23]: stress-11   stress-22   stress-33    ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3*/
double nostrs[20][26];
double srst[6];
double s123[12];
double gpcod[3]; /* natural coordinates of g.p.*/
/*-------------------------------------------  for eas elements only ---*/
int    l1, l3, ihyb, cc;
double xjm0[3][3];
double ehdis[3][10],fi[6][6],ff[6][6];
double det0, det1;
double bn1[3][10];
double disd1[9];
double fieh[30];
double epsh[6];

static double **estiflo;       
static ARRAY    estiflo_a; /* local element stiffness matrix ke for eas */   



static double **estif9;       
static ARRAY    estif9_a;   /* element stiffness matrix ke for eas */   

static ARRAY    D_a;      /* material tensor */     
static double **D;         
static ARRAY    funct_a;  /* shape functions */    
static double  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static double **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static double **xjm;         
static ARRAY    bop_a;    /* B-operator */   
static double **bop;       
static ARRAY    bnop_a;   /* BN-operator */   
static double **bn;       
static double **estif;    /* element stiffness matrix ke */

double det;

int    iform;             /* index for nonlinear formulation of element */
int    calstr;            /* flag for stress calculation                */
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
  calstr = 1;
  newval = 0;
}
else if(init==0)/*?!*/
{
  newval = 1; /* foam plasticity with large deformations ?! */
}
/*------------------------------------------- integration parameters ---*/
c1intg(ele,data,1);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
amzero(&estiflo_a);


for (i=0; i<81; i++) fielo[i] = 0.0;
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
      c1_call_mat(ele, mat,bop,xjm,ip,F,strain,D,disd,g,gi,istore,newval);
      /*-----------------------  calculate stresses - postprocessing ---*/
      if(calstr==1) 
      {
        for (i=0; i<6; i++) s123[i]=F[i];
        for (i=0; i<6; i++) srst[i]=F[i];
        c1trss2local(srst, gi);

        c1_cstr (srst ,s123, ele->e.c1->stress[0].gpstrs[ip]);
        /*------------- global coordinates of actual  gaussian point ---*/
        c1gcor  (funct,xyze, iel  ,ele->e.c1->stress[0].gpcoor[ip]);
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
    c1_sext(ele->e.c1->stress[0].npstrs, funct, deriv,xjm,xyze,
            ele->e.c1->stress[0].gpstrs, data->xgrr, data->xgss, data->xgtt, nir,nis,nit,iel);
    goto end;
  }
/*----------------------------------------------------------------------*/
  if(ihyb>0)
  {
    c1rkefi(ele, estif9, estiflo, fieh, fielo, l1);
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
\param    F   DOUBLE  (i)   force vector integral (stress-resultants) 
\param  fac   DOUBLE  (i)   multiplier for numerical integration      
\param  bop   DOUBLE  (i)   b-operator matrix                         
\param   nd   INT     (i)   total number degrees of freedom of element
\param  fie   DOUBLE  (o)   internal force vector                     

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1fi( double  *F,   /*  force vector integral (stress-resultants)  */
           double   fac, /*  multiplier for numerical integration       */
           double **bop, /*  b-operator matrix                          */
           int      nd,  /*  total number degrees of freedom of element */
           double  *fie) /*  internal force vector                      */
{
/*----------------------------------------------------------------------*/
int i,j,k;
double n11,n22,n33,n12,n23,n31;
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
\param       xjm   DOUBLE  (i)  jacobian matrix r,s,t-direction          
\param   g[6][6]   DOUBLE  (o)  transformation matrix s(glob)=g*s(loc)   
\param  gi[6][6]   DOUBLE  (o)  inverse of g          s(loc) =gi*s(glob) 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1tram( double **xjm,   /* jacobian matrix r,s,t-direction         */
             double g[6][6], /* transformation matrix s(glob)=g*s(loc)  */
             double gi[6][6])/* inverse of g          s(loc) =gi*s(glob)*/
{
/*----------------------------------------------------------------------*/
int i,j,k;
double dm;
double x1r, x2r, x3r, x1s, x2s, x3s;
double x1n, x2n, x3n, x1c, x2c, x3c;
double a11, a21, a31, a12, a22, a32, a13, a23, a33;
double dum[3][3];
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
void c1gld(double **d,     /* material matrix                           */
           double g[6][6]) /* transformation matrix                     */
{
/*----------------------------------------------------------------------*/
int i,j,k;
double dgt[6][6];
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
\param        *s   DOUBLE  (o)  stress vector to be transformed   
\param   g[6][6]   DOUBLE  (i)  inverse of transformation matrix  

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1trss2local(double *s,       /* stress vector to be transformed   */
                  double gi[6][6]) /* inverse of transformation matrix  */
{
/*----------------------------------------------------------------------*/
int i,j;
double sh[6];
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
\param        *s   DOUBLE  (o)  stress vector to be transformed   
\param   g[6][6]   DOUBLE  (i)  transformation matrix  

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1trss2global(double *s,       /* stress vector to be transformed  */
                   double g[6][6])  /* transformation matrix            */
{
/*----------------------------------------------------------------------*/
int i,j;
double sh[6];
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
#endif
/*! @} (documentation module close)*/
