/*!----------------------------------------------------------------------
\file
\brief contains the routines for stress evaluation for a 3D hex element

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
\brief compute principal stresses

<pre>                                                              al 06/02
This routine evaluates principal stresses at given gauss point 
for a 3D hex element.

</pre>
\param  srst   DOUBLE*  (i)   stresses at given gauss point           
\param  s123   DOUBLE*  (o)   principal stresses and direction at g.p.

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1pstr(DOUBLE    *srst,/* stresses at given gauss point           */
            DOUBLE    *s123 /* principal stresses and direction at g.p.*/
            )
{
/*----------------------------------------------------------------------*/
INT i, j, k, ci, cj;
INT i1[3];
DOUBLE a[9], v[9], vn;
DOUBLE pi, p;
DOUBLE c180=180.;
DOUBLE wr,ws,wt;
DOUBLE strmin  = 1.0E-9;
DOUBLE strzero = 0.;
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
for a 3D hex element (integretion point values).

pstrs[ 0.. 5] [stress-rr stress-ss stress-tt stress-rs stress-st stress-tr ]
pstrs[ 6..11] [stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz ]
pstrs[12..14] [stress-11 stress-22 stress-33 ]
pstrs[15..23] [ang-r1  ang-s1  ang-t1  ang-r2  ang-s2  ang-t2  ang-r3  ang-s3  ang-t3 ]
pstrs[24..26] [x y z ] -global coordinates of integration points

</pre>
\param srst   DOUBLE*   (i)   stresses at given gauss point           
\param s123   DOUBLE*   (i)   principal stresses and direction at g.p.
\param pstrs  DOUBLE*   (o)   postprocessing stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_cstr(DOUBLE    *srst,
             DOUBLE    *s123,
             DOUBLE    *pstrs
            )
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE aux[6];
DOUBLE strmin  = 1.0E-9;
DOUBLE strzero = 0.;
DOUBLE seqv;
#ifdef DEBUG 
dstrc_enter("c1_cstr");
#endif
/*------------------------------------------- mises effective stress ---*/
seqv = ( srst[0]*srst[0] + srst[1]*srst[1] + srst[2]*srst[2])
     + 2.0*(srst[3]*srst[3] + srst[4]*srst[4] + srst[5]*srst[5]);
seqv=sqrt(seqv);
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
  pstrs[12]   =  s123[ 0]; /*sig-  i   )                        */
  pstrs[13]   =  s123[ 1]; /*sig- ii   ) --> principal stresses */
  pstrs[14]   =  s123[ 2]; /*sig-iii   )                        */
  pstrs[15]   =  s123[ 3]; /*alpha  (r,i)  )                    */
  pstrs[16]   =  s123[ 4]; /*alpha  (s,i)  )                    */
  pstrs[17]   =  s123[ 5]; /*alpha  (t,i)  )                    */
  pstrs[18]   =  s123[ 6]; /*alpha (r,ii)  )                    */
  pstrs[19]   =  s123[ 7]; /*alpha (s,ii)  ) -->   directions   */
  pstrs[20]   =  s123[ 8]; /*alpha (t,ii)  )                    */
  pstrs[21]   =  s123[ 9]; /*alpha(r,iii)  )                    */
  pstrs[22]   =  s123[10]; /*alpha(s,iii)  )                    */
  pstrs[23]   =  s123[11]; /*alpha(t,iii)  )                    */
  pstrs[24]   =  seqv;     /*equivalent stress                  */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1_cstr */

/*!----------------------------------------------------------------------
\brief calculate global coordinates of integration points 
       referring to the natural ones

<pre>                                                              al 06/02
This routine calculates global coordinates referring to the natural ones 
for a 3D hex element.

</pre>
\param funct  DOUBLE*  (i)   value of form functions           
\param  xyze  DOUBLE*  (i)   element coordinates               
\param   iel      INT  (i)   number of nodes                   
\param gpcod  DOUBLE*  (o)   global coordinates of actual point

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1gcor( 
             DOUBLE     *funct, /* value of form functions              */
             DOUBLE     *xyze,  /* element coordinates                  */
             INT         iel,   /* number of nodes                      */
             DOUBLE     *gpcod  /* global coordinates of actual point   */
            )
{
/*----------------------------------------------------------------------*/
INT i,k,pc;
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

/*!----------------------------------------------------------------------
\brief returns local coordinates of element nodes in rst-coordinates 

<pre>                                                              al 06/02
This routine returns local coordinates of element nodes in rst-coordinates 
for a 3D hex element.

</pre>
\param   node     INT  (i)   element node           
\param    irs     INT  (i)   flag for r,s or t             

\warning There is nothing special to this routine
\return DOUBLE local coordinates of element node                                               
\sa calling: ---; called by: c1_cstr()

*----------------------------------------------------------------------*/
DOUBLE c1rsn (
             INT node,
             INT irs
             )
{
/*----------------------------------------------------------------------*/
static DOUBLE  xh8[24] = { 1., 1.,-1.,-1., 1., 1.,-1.,-1.,
                          -1., 1., 1.,-1.,-1., 1., 1.,-1.,
                          -1.,-1.,-1.,-1., 1., 1., 1., 1.};
static DOUBLE xh20[36] = { 
              1.,  0., -1.,  0.,  1.,  0., -1.,  0.,  1.,  1., -1., -1., 
              0.,  1.,  0., -1.,  0.,  1.,  0., -1., -1.,  1.,  1., -1., 
             -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0. };
/*----------------------------------------------------------------------*/
INT inode;
DOUBLE ret_val;
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
\param       n     INT  (i) order legendre polynomial
\param      zr     INT* (i)
\param       z  DOUBLE  (i) z-coordinate
\param   value  DOUBLE* (o) value at z

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1hxsm()

*----------------------------------------------------------------------*/
void c1lgpl (
              INT         i,
              INT         n,
              DOUBLE    *zr,  
              DOUBLE      z,
              DOUBLE *value
              )
{
/*----------------------------------------------------------------------*/
INT j;
DOUBLE zi, zj;
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
/*!----------------------------------------------------------------------
\brief subroutine of c1_sext

<pre>                                                              al 06/02
subroutine of c1_sext

</pre>
\param   nir,nis,nit     INT  (i) num GP in r/s/t direction
\param      rk,sk,tk     INT  (i) r,s,t -coordinates
\param      f[8][27]     DOUBLE  (i) original values on g.p.
\param            fp     DOUBLE* (o) extrapolated values
\param           xgr     DOUBLE* (i) coordinates of g.p.

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_sext()

*----------------------------------------------------------------------*/
void c1hxsm (
              INT nir,
              INT nis,
              INT nit,  
              DOUBLE rk,
              DOUBLE sk,
              DOUBLE tk,
              DOUBLE f[8][27],
              DOUBLE *fp,
              DOUBLE *xgr,
              DOUBLE *xgs,
              DOUBLE *xgt
              )
{
/*----------------------------------------------------------------------*/
INT i, j, k, ns, kkk, ngp;
DOUBLE xlr, xls, xlt;
#ifdef DEBUG 
dstrc_enter("c1rsn");
#endif
/*----------------------------------------------------------------------*/
  kkk=6+1;
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
This routine calculates principal element stresses and fill stress vector
for postprocessing for a 3D hex element.

pstrs[ 0.. 5]: x-coord. y-coord. z-coord.  stress-rr   stress-ss   stress-tt   stress-rs   stress-st   stress-tr  
pstrs[ 6..11]: x-coord. y-coord. z-coord. stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz 
pstrs[12..23]: stress-11 stress-22 stress-33 ang-r1 ang-s1 ang-t1 ang-r2 ang-s2 ang-t2 ang-r3 ang-s3 ang-t3
pstrs[    26]: equivalent stress
</pre>
\param  srst   DOUBLE*  (i)   element stresses at given gauss point           
\param  s123   DOUBLE*  (i)   principal stresses and direction at g.p.
\param  pstrs  DOUBLE*  (o)   postprocessing stresses

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_nstr(DOUBLE    *srst,     
             DOUBLE    *s123,
             DOUBLE    *pstrs
            )
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE aux[6];
DOUBLE strmin  = 1.0E-9;
DOUBLE strzero = 0.;
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
  pstrs[26]   =  srst[ 6];/* equivalent stress */
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

nostrs[element node] [ 0.. 5]: stress-rr stress-ss stress-tt stress-rs stress-st stress-tr  
nostrs[element node] [ 6..11]: stress-xx stress-yy stress-zz stress-xy stress-yz stress-xz
nostrs[element node] [12..23]: stress-11 stress-22 stress-33 ang-r1 ang-s1 ang-t1 ang-r2 ang-s2 ang-t2 ang-r3 ang-s3 ang-t3

</pre>
\param     nostrs  DOUBLE**  (o) element stresses extrapolated to the nodes
\param      funct  DOUBLE*   (i) shape functions
\param      deriv  DOUBLE**  (i) derivatives of the shape functions
\param        xjm  DOUBLE**  (i) the Jacobian matrix
\param       xyze  DOUBLE*   (i) element-node coordinates
\param   gpstress  DOUBLE**  (i) element stresses of integration points
\param        xgr  DOUBLE*   (i) local rst-coordinates of integration points
\param        xgs  DOUBLE*   (i) ..
\param        xgt  DOUBLE*   (i) ..
\param        nir     INT    (i) number of integration points in rst-direction
\param        nis     INT    (i) ..
\param        nit     INT    (i) ..
\param        iel     INT    (i) number of element nodes

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_sext(
            DOUBLE nostrs[20][26],
            DOUBLE     *funct,
            DOUBLE    **deriv,
            DOUBLE      **xjm,
            DOUBLE      *xyze,
            DOUBLE gpstress[27][26],
            DOUBLE       *xgr,
            DOUBLE       *xgs,
            DOUBLE       *xgt,
            INT           nir,
            INT           nis, 
            INT           nit,
            INT           iel 
            )
{
/*----------------------------------------------------------------------*/
INT nn,i,j,ngp;
DOUBLE g[6][6]; 
DOUBLE gi[6][6];
DOUBLE cnp1, cnp2, cnp3, det;
DOUBLE s123[12];
DOUBLE fgp[8][27];
DOUBLE fnp[8];
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
              fgp[6][i] = gpstress[i][24]; /* equivalent stress */
    for (j=0; j<6; j++){
              fgp[j][i] = gpstress[i][j+6];
  }}
/*----------------------------------------------------------------------*/
  for (nn=1; nn<=iel; nn++)
  {
    cnp1 = c1rsn (nn,1);
    cnp2 = c1rsn (nn,2);
    cnp3 = c1rsn (nn,3);

    c1hxsm (nir,nis,nit,cnp1,cnp2,cnp3,fgp,fnp,xgr,xgs,xgt);
    /*------------------------- retransform stresses to local system ---*/
    c1_funct_deriv(funct,deriv,cnp1,cnp2,cnp3,iel,1);

    c1_jaco (deriv,xjm,&det,xyze,iel);

    c1tram (xjm,g,gi);

    for (i=0; i<12; i++) s123[i] = 0.;
    for (i=0; i< 7; i++) s123[i] = fnp[i];

    c1trss2local (fnp,gi);
   /*------------------------- store nodal stresses in array stresk ---*/
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
