/*!----------------------------------------------------------------------
\file
\brief contains the routines forstress evaluation for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief compute principal stresses

<pre>                                                              al 06/02
This routine evaluates principal stresses at given gauss point 
for a 3D hex element.

</pre>
\param **srst   DOUBLE  (i)   stresses at given gauss point           
\param  *s123   DOUBLE  (o)   principal stresses and direction at g.p.

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | stress  (1) = sig-rr                                                 |
 | stress  (2) = sig-ss                                                 |
 | stress  (3) = sig-tt                                                 |
 | stress  (4) = sig-rs                                                 |
 | stress  (5) = sig-st                                                 |
 | stress  (6) = sig-tr                                                 |
 | stress  (7) = sig-  i   )                                            |
 | stress  (8) = sig- ii   ) --> principal stresses                     |
 | stress  (9) = sig-iii   )                                            |
 | stress (10) = alpha  (r,i)  )                                        |
 | stress (11) = alpha  (s,i)  )                                        |
 | stress (12) = alpha  (t,i)  )                                        |
 | stress (13) = alpha (r,ii)  )                                        |
 | stress (14) = alpha (s,ii)  ) -->   directions                       |
 | stress (15) = alpha (t,ii)  )                                        |
 | stress (16) = alpha(r,iii)  )                                        |
 | stress (17) = alpha(s,iii)  )                                        |
 | stress (18) = alpha(t,iii)  )                                        |
 *----------------------------------------------------------------------*/
void c1pstr(double    *srst,/* stresses at given gauss point            */
            double    *s123 /* principal stresses and direction at g.p. */
            )
{
/*----------------------------------------------------------------------*/
int i, j, k, ci, cj;
int i1[3];
double a[9], v[9], vn;
double pi, p;
double c180=180.;
double wr,ws,wt;
double strmin  = 1.0E-9;
double strzero = 0.;
#ifdef DEBUG 
dstrc_enter("c1pstr");
#endif
/*--------- initialize a symmetric matrix with upper triangular part ---*/
  a[0]=srst[0];/* a[0][0] */
  a[1]=srst[3];/* a[1][0] */
  a[2]=srst[5];/* a[2][0] */
  a[3]=srst[3];/* a[0][1] */
  a[4]=srst[1];/* a[1][1] */
  a[5]=srst[4];/* a[2][1] */
  a[6]=srst[5];/* a[0][2] */
  a[7]=srst[4];/* a[1][2] */
  a[8]=srst[2];/* a[2][2] */
/*--- calculate the principal stresses/directions by jacobi method ----*/
  for (i=0; i<9; i++) v[i] = 0.;
  c1jacb (a,v);

  pi=acos(-1.);
  
  for (i=0; i<3; i++) i1[i]=i+1;
  
  ci=0;
  cj=0;
  for (i=0; i<2; i++)
  {
    p=a[ci];
    cj=ci+4;
    for (j=i+1; j<3; j++)
    {
       if(p<a[cj]) 
       {
         a[ci]=a[cj];
         a[cj]=p;
         i1[i]=j+1;
         i1[j]=i+1;
         i1[3-i-j] = 4 - i - j;
       }
       cj+=4;
    }
    ci+=4;
  }
/*------------------------------------------- initialize vector s123 ---*/
  for (i=0; i<12; i++) s123[i] = 0.;

  ci=-4;
  for (i=0; i<3; i++)
  {
    ci+=4;
    k=i1[i];
    cj=(k-1)*3;
    vn=sqrt(v[0+(k-1)*3]*v[0+(k-1)*3]+v[1+(k-1)*3]*v[1+(k-1)*3]+v[2+(k-1)*3]*v[2+(k-1)*3]);
    v[0+(k-1)*3] = v[0+(k-1)*3]/vn;
    v[1+(k-1)*3] = v[1+(k-1)*3]/vn;
    v[2+(k-1)*3] = v[2+(k-1)*3]/vn;
    wr=acos(v[0+(k-1)*3])*c180/pi;
    ws=acos(v[1+(k-1)*3])*c180/pi;
    wt=acos(v[2+(k-1)*3])*c180/pi;
    s123[i]=a[ci];
    if(fabs(s123[i])<=strmin) continue;
    s123[3*(i+1)+0]=wr;
    s123[3*(i+1)+1]=ws;
    s123[3*(i+1)+2]=wt;
  }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1pstr */
/*----------------------------------------------------------------------*/

/*!----------------------------------------------------------------------
\brief calculate element stresses for postprocessing

<pre>                                                              al 06/02
This routine calculates element stresses for postprocessing 
for a 3D hex element.

</pre>
\param **srst   DOUBLE  (i)   stresses at given gauss point           
\param  *s123   DOUBLE  (i)   principal stresses and direction at g.p.
\param  *pstrs  DOUBLE  (o)   postprocessing stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_cstr(double    *srst,
             double    *s123,
             double    *pstrs
            )
{
/*----------------------------------------------------------------------*/
int i;
double aux[6];
double strmin  = 1.0E-9;
double strzero = 0.;
#ifdef DEBUG 
dstrc_enter("c1_cstr");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++) if(fabs(s123[i])<strmin) s123[i]=strzero;
  for (i=0; i<6; i++) if(fabs(srst[i])<strmin) srst[i]=strzero;
  pstrs[ 0]   =  srst[ 0];
  pstrs[ 1]   =  srst[ 1];
  pstrs[ 2]   =  srst[ 2];
  pstrs[ 3]   =  srst[ 3];
  pstrs[ 4]   =  srst[ 4];
  pstrs[ 5]   =  srst[ 5];
  pstrs[ 6]   =  s123[ 0];
  pstrs[ 7]   =  s123[ 1];
  pstrs[ 8]   =  s123[ 2];
  pstrs[ 9]   =  s123[ 3];
  pstrs[10]   =  s123[ 4];
  pstrs[11]   =  s123[ 5];
  pstrs[12]   =  s123[ 6];
  pstrs[13]   =  s123[ 7];
  pstrs[14]   =  s123[ 8];
  pstrs[15]   =  s123[ 9];
  pstrs[16]   =  s123[10];
  pstrs[17]   =  s123[11];
/*--------------------------------- evaluation of principal stresses ---*/
  for (i=0; i<6; i++) aux[i] = s123[i];
  c1pstr (aux,s123); 
  for (i=0; i<12; i++) if(fabs(s123[i])<strmin) s123[i]=strzero;
/*----------------------------------------------------------------------*/
  pstrs[12]   =  s123[ 0];
  pstrs[13]   =  s123[ 1];
  pstrs[14]   =  s123[ 2];
  pstrs[15]   =  s123[ 3];
  pstrs[16]   =  s123[ 4];
  pstrs[17]   =  s123[ 5];
  pstrs[18]   =  s123[ 6];
  pstrs[19]   =  s123[ 7];
  pstrs[20]   =  s123[ 8];
  pstrs[21]   =  s123[ 9];
  pstrs[22]   =  s123[10];
  pstrs[23]   =  s123[11];
  pstrs[24]   =  srst[ 6];
  pstrs[25]   =  srst[ 7];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1_cstr */

/*!----------------------------------------------------------------------
\brief calculate global coordinates referring to the natural ones

<pre>                                                              al 06/02
This routine calculates global coordinates referring to the natural ones 
for a 3D hex element.

</pre>
\param *funct  DOUBLE  (i)   value of form functions           
\param  *xyze  DOUBLE  (i)   element coordinates               
\param    iel     INT  (o)   number of nodes                   
\param *gpcod  DOUBLE  (o)   global coordinates of actual point

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1gcor( 
             double     *funct, /* value of form functions              */
             double     *xyze,  /* element coordinates                  */
             int         iel,   /* number of nodes                      */
             double     *gpcod  /* global coordinates of actual point   */
            )
{
/*----------------------------------------------------------------------*/
int i,j,k,l,pc;
double dum;
#ifdef DEBUG 
dstrc_enter("c1gcor");
#endif
/*----------------------------------------------------------------------*/
  for (k=0; k<3; k++)
  {
     pc=0;
     gpcod[k] = 0.;
     for (i=0; i<iel; i++)
     {
       gpcod[k] += xyze[pc+k]*funct[i];
       pc+=3;
     }
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1gcor */
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double c1rsn (
             int node,
             int irs,
             int iel  /* number of nodes */
             )
{
/*----------------------------------------------------------------------*/
static double  xh8[24] = { 1., 1.,-1.,-1., 1., 1.,-1.,-1.,
                          -1., 1., 1.,-1.,-1., 1., 1.,-1.,
                          -1.,-1.,-1.,-1., 1., 1., 1., 1.};
static double xh20[36] = { 
              1.,  0., -1.,  0.,  1.,  0., -1.,  0.,  1.,  1., -1., -1., 
              0.,  1.,  0., -1.,  0.,  1.,  0., -1., -1.,  1.,  1., -1., 
             -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0. };
/*----------------------------------------------------------------------*/
int inode;
double ret_val;
#ifdef DEBUG 
dstrc_enter("c1rsn");
#endif
/*----------------------------------------------------------------------*/
  if(node<=8)
  {
    inode = node-1 + (irs-1)*8;
    ret_val=xh8[inode];
  }
  else if(node<=20)
  {
    inode=((node-8)-1)+ (irs-1)*12;
    ret_val=xh20[inode];
  }
  else
  {
   dserror("unknown number of nodes in hex element");
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return(ret_val);
} /* end of c1rsn */

/*!----------------------------------------------------------------------
\brief calculate nth order legendre polynomial of degree 'n' at 'z'

<pre>                                                              al 06/02
This routine calculates nth order legendre polynomial of degree 'n' at 'z'
for a 3D hex element.

</pre>
\param       i     INT  (i)
\param       n     INT  (i)
\param     *zr     INT  (i)
\param       z  DOUBLE  (i)
\param  *value  DOUBLE  (o)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1lgpl (
              int         i,
              int         n,
              double    *zr,  
              double      z,
              double *value
              )
{
/*----------------------------------------------------------------------*/
int j;
double zi, zj;
#ifdef DEBUG 
dstrc_enter("c1lgpl");
#endif
/*----------------------------------------------------------------------*/
  zi = zr[i];
  *value = 1.;
  for (j=0; j<n; j++)
  {
    zj = zr[j];
    if(j==i) continue;
    *value *= (z-zj)/(zi-zj);
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1lgpl */
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void c1hxsm (
              int nir,
              int nis,
              int nit,  
              double rk,
              double sk,
              double tk,
              double f[8][27],
              double *fp,
              double *xgr,
              double *xgs,
              double *xgt
              )
{
/*----------------------------------------------------------------------*/
int i, j, k, ns, kkk, ngp;
double xlr, xls, xlt;
#ifdef DEBUG 
dstrc_enter("c1rsn");
#endif
/*----------------------------------------------------------------------*/
  kkk=6;
/*----------------------------------------------------------------------*/
  for (i=0; i<8; i++) fp[i] = 0.;
/*----------------------------------------------------------------------*/
  ngp=0;
/*----------------------------------------------------------------------*/
  for (i=0; i<nir; i++) {
    c1lgpl (i,nir,xgr,rk,&xlr);
    for (j=0; j<nis; j++) {
      c1lgpl (j,nis,xgs,sk,&xls);
      for (k=0; k<nit; k++) {
        c1lgpl (k,nit,xgt,tk,&xlt);
        for (ns=0; ns<kkk; ns++)  fp[ns] += xlr*xls*xlt*f[ns][ngp];
        ngp = ngp + 1;
  }}}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1hxsm */

/*!----------------------------------------------------------------------
\brief calculate element stresses for postprocessing

<pre>                                                              al 06/02
This routine calculates element stresses for postprocessing 
for a 3D hex element.

</pre>
\param **srst   DOUBLE  (i)   stresses at given gauss point           
\param  *s123   DOUBLE  (i)   principal stresses and direction at g.p.
\param  *pstrs  DOUBLE  (o)   postprocessing stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_nstr(double    *srst,     /* element vector with stress resultants */
             double    *s123,
             double    *pstrs
            )
{
/*----------------------------------------------------------------------*/
int i;
double aux[6];
double strmin  = 1.0E-9;
double strzero = 0.;
#ifdef DEBUG 
dstrc_enter("c1_nstr");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<6; i++) if(fabs(s123[i])<strmin) s123[i]=strzero;
  for (i=0; i<6; i++) if(fabs(srst[i])<strmin) srst[i]=strzero;
  pstrs[ 0]   =  srst[ 0];
  pstrs[ 1]   =  srst[ 1];
  pstrs[ 2]   =  srst[ 2];
  pstrs[ 3]   =  srst[ 3];
  pstrs[ 4]   =  srst[ 4];
  pstrs[ 5]   =  srst[ 5];
  pstrs[ 6]   =  s123[ 0];
  pstrs[ 7]   =  s123[ 1];
  pstrs[ 8]   =  s123[ 2];
  pstrs[ 9]   =  s123[ 3];
  pstrs[10]   =  s123[ 4];
  pstrs[11]   =  s123[ 5];
  pstrs[12]   =  s123[ 6];
  pstrs[13]   =  s123[ 7];
  pstrs[14]   =  s123[ 8];
  pstrs[15]   =  s123[ 9];
  pstrs[16]   =  s123[10];
  pstrs[17]   =  s123[11];
/*--------------------------------- evaluation of principal stresses ---*/
  for (i=0; i<6; i++) aux[i] = s123[i];
  c1pstr (aux,s123); 
  for (i=0; i<12; i++) if(fabs(s123[i])<strmin) s123[i]=strzero;
/*----------------------------------------------------------------------*/
  pstrs[12]   =  s123[ 0];
  pstrs[13]   =  s123[ 1];
  pstrs[14]   =  s123[ 2];
  pstrs[15]   =  s123[ 3];
  pstrs[16]   =  s123[ 4];
  pstrs[17]   =  s123[ 5];
  pstrs[18]   =  s123[ 6];
  pstrs[19]   =  s123[ 7];
  pstrs[20]   =  s123[ 8];
  pstrs[21]   =  s123[ 9];
  pstrs[22]   =  s123[10];
  pstrs[23]   =  s123[11];
  pstrs[24]   =  srst[ 6];
  pstrs[25]   =  srst[ 7];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1_nstr */

/*!----------------------------------------------------------------------
\brief extrapolation of stress from gauss points to nodal points

<pre>                                                              al 06/02
This routine extrapolates stresses from gauss points to nodal points 
for a 3D hex element.

</pre>
\param   **nostrs  DOUBLE  (o)
\param     *funct  DOUBLE  (i)
\param    **deriv  DOUBLE  (i)
\param      **xjm  DOUBLE  (i)
\param      *xyze  DOUBLE  (i)
\param **gpstress  DOUBLE  (i)
\param       *xgr  DOUBLE  (i)
\param       *xgs  DOUBLE  (i)
\param       *xgt  DOUBLE  (i)
\param        nir     INT  (i)
\param        nis     INT  (i)
\param        nit     INT  (i)
\param        iel     INT  (i)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_sext(
            double   **nostrs,
            double     *funct,
            double    **deriv,
            double      **xjm,
            double      *xyze,
            double **gpstress,
            double       *xgr,
            double       *xgs,
            double       *xgt,
            int           nir,
            int           nis, 
            int           nit,
            int           iel /* number of nodes */
            )
{
/*----------------------------------------------------------------------*/
int nn,i,j,k,l,ngp;
double g[6][6]; 
double gi[6][6];
double cnp1, cnp2, cnp3, det;
double s123[12];
double fgp[8][27];
double fnp[8];
double dum;
#ifdef DEBUG 
dstrc_enter("c1_sext");
#endif
/*----------------------------------------------------------------------*/
  ngp = nir*nis*nit;
/*----------------------------------------------------------------------*/
  for (i=0; i<8; i++){
    for (j=0; j<27; j++){
       fgp[i][j] = 0.;}}

  for (i=0; i<ngp; i++){
    for (j=0; j<6; j++){
              fgp[j][i] = gpstress[i][j+6];
  }}
/*----------------------------------------------------------------------*/
  for (nn=1; nn<=iel; nn++)
  {
    cnp1 = c1rsn (nn,1,iel);
    cnp2 = c1rsn (nn,2,iel);
    cnp3 = c1rsn (nn,3,iel);

    c1hxsm (nir,nis,nit,cnp1,cnp2,cnp3,fgp,fnp,xgr,xgs,xgt);
    /*------------------------- retransform stresses to local system ---*/
    c1_funct_deriv(funct,deriv,cnp1,cnp2,cnp3,iel,1);

    c1_jaco (deriv,xjm,&det,xyze,iel);

    c1tram (xjm,g,gi);

    for (i=0; i<12; i++) s123[i] = 0.;
    for (i=0; i< 6; i++) s123[i] = fnp[i];

    c1trss2local (fnp,gi);
   /*------------------------- store nodal stresses in array stresk ---*/
    /*nostrs[nn] [ 0.. 5]: x-coord.    y-coord.    z-coord.     stress-rr   stress-ss   stress-tt   stress-rs   stress-st   stress-tr*/  
    /*nostrs[nn] [ 6..11]: x-coord.    y-coord.    z-coord.     stress-xx   stress-yy   stress-zz   stress-xy   stress-yz   stress-xz*/ 
    /*nostrs[nn] [12..23]: stress-11   stress-22   stress-33    ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3*/
    c1_nstr (fnp,s123,nostrs[nn-1]);
    }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_sext */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
