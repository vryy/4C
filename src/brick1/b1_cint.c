#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | Prototypes                                           a.lipka 6/01    |
 *----------------------------------------------------------------------*/
void c1cintg(int  nir, int  nis, int  nit);
void c1hexa(int  iel, double r, double s, double t);
void c1jaco(ELEMENT *actele);
void c1bop(int  iel);
void c1matli(double ym,double pv);
void c1mat(ELEMENT *actele);
void c1vke(ELEMENT *actele,ARRAY estiff);

double funct[8];
double deriv[3][8];
double xgr[3], xgs[3], xgt[3];        /* coordinates              */
double wgr[3], wgs[3], wgt[3];        /* weighting factors        */
double b[6][3][8];
double d[6][6];
double det;
double xjm[3][3];
double s[24][24];
double fac;
/*----------------------------------------------------------------------*
 | HEX-ELEMENT    INTEGRATION OVER THE ELEMENT VOLUME   a.lipka 5/01    |
 *----------------------------------------------------------------------*/
void c1cint(ELEMENT *actele,ARRAY estiff)
{
/*----------------------------------------------------------------------*/
  int i,k;
  int nir, nis, nit;
  int lr, ls, lt;
  int iel;
  double facr, facs, fact;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1cint");
  #endif
/*----------------------------------------------------------------------*/
  nir = actele->e.b1->nGP[0];                                  /* gauss-integration points */
  nis = actele->e.b1->nGP[1];
  nit = actele->e.b1->nGP[2];
  iel = actele->numnp;
  
/*------------------------------------------ initialize stiffness matrix*/
  for(i=0;i<24;i++)
    for(k=0;k<24;k++)
        s[i][k] = 0.0;
/*----------------------------------------------------------------------*/
  
  c1cintg(nir, nis, nit);                   
/*-------------------------------------------------- LOOP IN R-DIRECTION*/
  for (lr=0; lr<nir; lr++)
  {/*lr*/
/*-------------------------------------------------- LOOP IN S-DIRECTION*/
  for (ls=0; ls<nis; ls++)
  {/*ls*/
/*-------------------------------------------------- LOOP IN T-DIRECTION*/
  for (lt=0; lt<nit; lt++)
  {/*lt*/
/*----------------------------------------------------------------------*/
    c1hexa (iel,xgr[lr],xgs[ls],xgt[lt]);
/*--------------------------------------------- EVALUATE JACOBIAN MATRIX*/
    c1jaco (actele);
/*------------------------------------- FORMULATION OF OPERATOR MATRIX B*/
    c1bop  (iel);
/*------------------------------------------ EVALUATE INTEGRATION FACTOR*/
    facr = wgr[lr];
    facs = wgs[ls];
    fact = wgt[lt];
    fac=facr*facs*fact*det;
/*---------------------- CONTROL PROGRAM FOR FORMULATION OF MATERIAL LAW*/
    c1mat (actele);
/*------------------------------------------- USUAL STIFFNESS-MATRIX KE */
    c1vke (actele, estiff);

/*----------------------------------------------------------------------*/
 }}}/*lrst*/
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1cint */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | COORDINATES AND WEIGHTING FACTORS OF GAUSS-INTEGRATION-POINTS FOR    |
 | NUMERICAL INTEGRATION                                                |
 ----------------------------------------------------------------------*/
void c1cintg(int  nir, int  nis, int  nit)
{
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1cintg");
  #endif
/*----------------------------------------------------------------------*/
 if(nir==2&&nis==2&&nit==2)
 {
   xgr[0]= -0.5773502691896;
   xgr[1]=  0.5773502691896;
   wgr[0]=  1.0            ;
   wgr[1]=  1.0            ;
   xgs[0]= -0.5773502691896;
   xgs[1]=  0.5773502691896;
   wgs[0]=  1.0            ;
   wgs[1]=  1.0            ;
   xgt[0]= -0.5773502691896;
   xgt[1]=  0.5773502691896;
   wgt[0]=  1.0            ;
   wgt[1]=  1.0            ;
 }
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1cintg */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR EVALUATION OF SHAPE FUNCTIONS AND THEIR NATURAL          |
 | DERIVATIVES WITH RESPECT TO  R/S/T  FOR  H E X A H E D R O N S       |
 ----------------------------------------------------------------------*/
void c1hexa(int  iel, double r, double s, double t)
{
/*----------------------------------------------------------------------*/
  double rp, rm, sp, sm, tp, tm;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1hexa");
  #endif
/*------------------------------------------- FORM BASIC FUNCTION VALUES*/
  rp=1.0 + r;
  rm=1.0 - r;
  sp=1.0 + s;
  sm=1.0 - s;
  tp=1.0 + t;
  tm=1.0 - t;
/*----------------------------------------------------------------------*/
 if(iel==8)
 {
    /*---------------------------------------------------------------+
    |  L I N E A R     SHAPE FUNCTIONS AND THEIR NATURAL DERIVATIVES |
    |  SERENDIPITY     ( 8-NODED ELEMENT )                           |
    +---------------------------------------------------------------*/

   funct[0]=0.125*rp*sm*tm;
   funct[1]=0.125*rp*sp*tm;
   funct[2]=0.125*rm*sp*tm;
   funct[3]=0.125*rm*sm*tm;
   funct[4]=0.125*rp*sm*tp;
   funct[5]=0.125*rp*sp*tp;
   funct[6]=0.125*rm*sp*tp;
   funct[7]=0.125*rm*sm*tp;
/*------------------------------------------------ DERIVATIVE EVALUATION*/
   deriv[0][0]= 0.125*sm*tm;
   deriv[0][1]= 0.125*sp*tm;
   deriv[0][2]=-deriv[0][1];
   deriv[0][3]=-deriv[0][0];
   deriv[0][4]= 0.125*sm*tp;
   deriv[0][5]= 0.125*sp*tp;
   deriv[0][6]=-deriv[0][5];
   deriv[0][7]=-deriv[0][4];
   deriv[1][0]=-0.125*tm*rp;
   deriv[1][1]=-deriv[1][0];
   deriv[1][2]= 0.125*tm*rm;
   deriv[1][3]=-deriv[1][2];
   deriv[1][4]=-0.125*tp*rp;
   deriv[1][5]=-deriv[1][4];
   deriv[1][6]= 0.125*tp*rm;
   deriv[1][7]=-deriv[1][6];
   deriv[2][0]=-0.125*rp*sm;
   deriv[2][1]=-0.125*rp*sp;
   deriv[2][2]=-0.125*rm*sp;
   deriv[2][3]=-0.125*rm*sm;
   deriv[2][4]=-deriv[2][0];
   deriv[2][5]=-deriv[2][1];
   deriv[2][6]=-deriv[2][2];
   deriv[2][7]=-deriv[2][3];
 }
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1hexa */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR EVALUATION OF SHAPE FUNCTIONS AND THEIR NATURAL          |
 | DERIVATIVES WITH RESPECT TO  R/S/T  FOR  H E X A H E D R O N S       |
 ----------------------------------------------------------------------*/
void c1jaco(ELEMENT *actele)
{
/*----------------------------------------------------------------------*/
  int i,j,k,iel;
  double c;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1jaco");
  #endif
/*----------------------------------------------------------------------*/
  iel = actele->numnp;
/*--------------------------- DETERMINE JACOBIAN MATRIX AT POINT (R,S,T)*/
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      c = 0.0;
      for (k=0; k<iel; k++)
      {
        c+=deriv[i][k]* (actele->node[k]->x[j]);
      }
      xjm[i][j]=c;
    }
  }
/*---------- COMPUTE DETERMINANT OF THE JACOBIAN MATRIX AT POINT (R,S,T)*/
   det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
         xjm[0][1]*xjm[1][2]*xjm[2][0]+
         xjm[0][2]*xjm[1][0]*xjm[2][1]-
         xjm[0][2]*xjm[1][1]*xjm[2][0]-
         xjm[0][0]*xjm[1][2]*xjm[2][1]-
         xjm[0][1]*xjm[1][0]*xjm[2][2] ;
  if (det<=0.0) 
      dserror("JACOBIAN DETERMINANT AT ELEMENT");
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1jaco */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR FORMULATION OF OPERATOR-MATRIX B FOR BRICK ELEMENT       |
 ----------------------------------------------------------------------*/
void c1bop(int iel)
{
/*----------------------------------------------------------------------*/
  int i;
  double dum;  
  double x1r, x2r, x3r, x1s, x2s, x3s, x1t, x2t, x3t;
  double xi11, xi12, xi13, xi21, xi22, xi23, xi31, xi32, xi33;
  double hr, hs, ht;
  double h1, h2, h3;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1bop");
  #endif
/*----------------------------------------------------------------------*/
  x1r = xjm[0][0];
  x2r = xjm[0][1];
  x3r = xjm[0][2];
  x1s = xjm[1][0];
  x2s = xjm[1][1];
  x3s = xjm[1][2];
  x1t = xjm[2][0];
  x2t = xjm[2][1];
  x3t = xjm[2][2];
 
  dum=1.0/det;

  xi11=dum*(x2s*x3t - x2t*x3s);
  xi12=dum*(x3r*x2t - x2r*x3t);
  xi13=dum*(x2r*x3s - x3r*x2s);
  xi21=dum*(x3s*x1t - x3t*x1s);
  xi22=dum*(x1r*x3t - x3r*x1t);
  xi23=dum*(x3r*x1s - x1r*x3s);
  xi31=dum*(x1s*x2t - x1t*x2s);
  xi32=dum*(x2r*x1t - x1r*x2t);
  xi33=dum*(x1r*x2s - x2r*x1s);
/*----------------------------------------------------------------------*/
  for (i=0; i<iel; i++)         
  {
    hr   = deriv[0][i];
    hs   = deriv[1][i];
    ht   = deriv[2][i];

    h1 = xi11*hr + xi12*hs + xi13*ht;
    h2 = xi21*hr + xi22*hs + xi23*ht;
    h3 = xi31*hr + xi32*hs + xi33*ht;

    b[0][0][i] = h1 ;
    b[0][1][i] = 0.0;
    b[0][2][i] = 0.0;
    b[1][0][i] = 0.0;
    b[1][1][i] = h2 ;
    b[1][2][i] = 0.0;
    b[2][0][i] = 0.0;
    b[2][1][i] = 0.0;
    b[2][2][i] = h3 ;
    b[3][0][i] = h2 ;
    b[3][1][i] = h1 ;
    b[3][2][i] = 0.0;
    b[4][0][i] = 0.0;
    b[4][1][i] = h3 ;
    b[4][2][i] = h2 ;
    b[5][0][i] = h3 ;
    b[5][1][i] = 0.0;
    b[5][2][i] = h1 ;
  }   
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1bop */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR FORMULATION OF MATERIAL LAW                              |
 ----------------------------------------------------------------------*/
void c1mat(ELEMENT *actele)
{
/*----------------------------------------------------------------------*/
  int    nmat;
  double ym,pv;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1mat");
  #endif
/*------------------------------------------- SELECT PROPER MATERIAL LAW*/
  nmat = actele->mat;
  if(mat[0].mattyp==m_pl_por_mises)
  {
    ym = mat[nmat-1].m.pl_por_mises->youngs;
    pv = mat[nmat-1].m.pl_por_mises->possionratio;
  }
  c1matli (ym, pv);                                   /*LINEAR ISOTROPIC*/
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1mat */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR FORMULATION OF MATERIAL LAW   - LINEAR ISOTROPIC         |
 ----------------------------------------------------------------------*/
void c1matli(double ym,double pv)
{
/*----------------------------------------------------------------------*/
  int    i, j;
  double d1,d2,d3;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1matli");
  #endif
/*------------------------------------------- INITIALIZE MATERIAL MATRIX*/
  for (i=0; i<6; i++)         
    for (j=0; j<6; j++)         
      d[i][j] = 0.0;
/*--------------------------------------- EVALUATE BASIC MATERIAL VALUES*/
  d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
  d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
  d3=ym/((1.0 + pv)*2.0);
/*----------------------------------------- SET VALUES IN MATERIAL-MATRIX*/
  d[0][0]=d1;
  d[0][1]=d2;
  d[0][2]=d2;
  d[1][0]=d2;
  d[1][1]=d1;
  d[1][2]=d2;
  d[2][0]=d2;
  d[2][1]=d2;
  d[2][2]=d1;
  d[3][3]=d3;
  d[4][4]=d3;
  d[5][5]=d3;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1matli */
/*----------------------------------------------------------------------*
 | HEX-ELEMENT                                          a.lipka 5/01    |
 | PROGRAM FOR EVALUATION OF USUAL STIFFNESS MATRIX                     |
 ----------------------------------------------------------------------*/
void c1vke(ELEMENT *actele, ARRAY estiff)
{
/*----------------------------------------------------------------------*/
  int    i, j, k, iel, i1, i2, i3, il;
  int    iss,issl,kkk,iaa,ill,iab;
  double b11, b21, b31, b41, b51, b61,
         b12, b22, b32, b42, b52, b62,
         b13, b23, b33, b43, b53, b63;
  double g11, g21, g31, g41, g51, g61,
         g12, g22, g32, g42, g52, g62,
         g13, g23, g33, g43, g53, g63;
  double d1 , d2 , d3 , d4 , d5 , d6 ;
  double db11, db21, db31, db41, db51, db61;
  double db[6][3];
  int c1,c2,c3,r1;
 

  FILE *filep      ;
  char filename[50];
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("c1vke");
  #endif
/*----------------------------------------------------------------------*/
  iel = actele->numnp;
/*-------------------------------- LOOP OVER ALL NODAL POINTS AT ELEMENT*/
  iab=0;
  for (j=0; j<iel; j++)
  {/*lnp01*/         
     /*if("blockdiagonalmatrix")  kkk=j;*/
     /*if("full         matrix")*/kkk=0;
     
     b11=b[0][0][j];
     b21=b[1][0][j];
     b31=b[2][0][j];
     b41=b[3][0][j];
     b51=b[4][0][j];
     b61=b[5][0][j];
     b12=b[0][1][j];
     b22=b[1][1][j];
     b32=b[2][1][j];
     b42=b[3][1][j];
     b52=b[4][1][j];
     b62=b[5][1][j];
     b13=b[0][2][j];
     b23=b[1][2][j];
     b33=b[2][2][j];
     b43=b[3][2][j];
     b53=b[4][2][j];
     b63=b[5][2][j];
/*-------------------------------------LOOP OVER ALL MATERIAL COMPONENTS*/
    for (k=0; k<6; k++)
    {
      d1=d[k][0]*fac;
      d2=d[k][1]*fac;
      d3=d[k][2]*fac;
      d4=d[k][3]*fac;
      d5=d[k][4]*fac;
      d6=d[k][5]*fac;
/*-----------------------------------------------------------DB - MATRIX*/
      db[k][0] = d1*b11 + d2*b21 + d3*b31 + d4*b41 + d5*b51 + d6*b61;
      db[k][1] = d1*b12 + d2*b22 + d3*b32 + d4*b42 + d5*b52 + d6*b62;
      db[k][2] = d1*b13 + d2*b23 + d3*b33 + d4*b43 + d5*b53 + d6*b63;
    }
/*------------------------------------------------TRANSPOSED  B - MATRIX*/
      for (i=kkk; i<iel; i++)
      {/*lnp02*/
        g11=b[0][0][i];
        g21=b[1][0][i];
        g31=b[2][0][i];
        g41=b[3][0][i];
        g51=b[4][0][i];
        g61=b[5][0][i];
        g12=b[0][1][i];
        g22=b[1][1][i];
        g32=b[2][1][i];
        g42=b[3][1][i];
        g52=b[4][1][i];
        g62=b[5][1][i];
        g13=b[0][2][i];
        g23=b[1][2][i];
        g33=b[2][2][i];
        g43=b[3][2][i];
        g53=b[4][2][i];
        g63=b[5][2][i];

        c1= i*3;
        c2=c1+1;
        c3=c2+1;
/*----------------------------------------- LOOP OVER DEGREES OF FREEDOM*/
        for (il=0; il<3; il++)
        {
          r1 = il + iab;
 
          db11=db[0][il];
          db21=db[1][il];
          db31=db[2][il];
          db41=db[3][il];
          db51=db[4][il];
          db61=db[5][il];
 
          estiff.a.da[r1][c1] +=     g11*db11 + g21*db21 + g31*db31 + g41*db41
                              + g51*db51 + g61*db61;
          estiff.a.da[r1][c2] +=     g12*db11 + g22*db21 + g32*db31 + g42*db41
                              + g52*db51 + g62*db61;
          estiff.a.da[r1][c3] +=     g13*db11 + g23*db21 + g33*db31 + g43*db41
                              + g53*db51 + g63*db61;
        }
/*----------------------------------------------------------------------*/
      }/*lnp02*/
/*----------------------------------------------------------------------*/
      iab=iab+3;
  }/*lnp01*/
/*----------------------------------------------------------------------*/
  /* matrix: dof - dof connectivity */
  sprintf(filename,"kel_proc%d_1.txt",par.myrank);
  
 
  if ((filep = fopen(filename, "w")) != NULL)
  {
    fprintf(filep, "    ");
    for(i=0;i<24;i++)  fprintf(filep, " %8d",i);
    fprintf(filep, "\n");
    for(i=0;i<222;i++)  fprintf(filep, "_",i);
    fprintf(filep, "\n");

    for(i=0;i<24;i++)
    {
      fprintf(filep, "%4d |",i);
      for(k=0;k<24;k++)
        {
        fprintf(filep, " %8.4f",estiff.a.da[i][k]);
        }
      fprintf(filep, "\n");
    }
    fclose(filep);
  }

/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1vke */
/*----------------------------------------------------------------------*/
