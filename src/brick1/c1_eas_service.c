/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_disd' which calclate displacement
       derivatives for a 3D hex element

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
\brief determine the transformation matrix

<pre>                                                              al 06/02
This routine calcuates the transformation matrix for orthogonal
systems and its inverse jacobian matrix at point 0,0,0                                       |
for an 3D-hex-element.

</pre>
\param  fi[6][6]   DOUBLE  (i)   
\param  ff[6][6]   DOUBLE  (i)   
\param    **xjm0   DOUBLE  (o)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1t0( DOUBLE   fi[6][6],    
           DOUBLE   ff[6][6],    
           DOUBLE     **xjm0)  
{
/*----------------------------------------------------------------------*/
INT i,j,cc;
INT fdim=6;
DOUBLE ffi[36];
DOUBLE fii[36];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1t0");
#endif
/*----------------------------------------------------------------------*/
  ff[0][0] = xjm0[0][0]*xjm0[0][0];
  ff[0][1] = xjm0[1][0]*xjm0[1][0];
  ff[0][2] = xjm0[2][0]*xjm0[2][0];
  ff[0][3] = xjm0[0][0]*xjm0[1][0] + xjm0[0][0]*xjm0[1][0];
  ff[0][5] = xjm0[0][0]*xjm0[2][0] + xjm0[0][0]*xjm0[2][0];
  ff[0][4] = xjm0[1][0]*xjm0[2][0] + xjm0[1][0]*xjm0[2][0];
  ff[1][0] = xjm0[0][1]*xjm0[0][1];
  ff[1][1] = xjm0[1][1]*xjm0[1][1];
  ff[1][2] = xjm0[2][1]*xjm0[2][1];
  ff[1][3] = xjm0[0][1]*xjm0[1][1] + xjm0[0][1]*xjm0[1][1];
  ff[1][5] = xjm0[0][1]*xjm0[2][1] + xjm0[0][1]*xjm0[2][1];
  ff[1][4] = xjm0[1][1]*xjm0[2][1] + xjm0[1][1]*xjm0[2][1];
  ff[2][0] = xjm0[0][2]*xjm0[0][2];
  ff[2][1] = xjm0[1][2]*xjm0[1][2];
  ff[2][2] = xjm0[2][2]*xjm0[2][2];
  ff[2][3] = xjm0[0][2]*xjm0[1][2] + xjm0[0][2]*xjm0[1][2];
  ff[2][5] = xjm0[0][2]*xjm0[2][2] + xjm0[0][2]*xjm0[2][2];
  ff[2][4] = xjm0[1][2]*xjm0[2][2] + xjm0[1][2]*xjm0[2][2];

  ff[3][0] = xjm0[0][0]*xjm0[0][1];
  ff[3][1] = xjm0[1][0]*xjm0[1][1];
  ff[3][2] = xjm0[2][0]*xjm0[2][1];
  ff[3][3] = xjm0[0][0]*xjm0[1][1] + xjm0[0][1]*xjm0[1][0];
  ff[3][5] = xjm0[0][0]*xjm0[2][1] + xjm0[0][1]*xjm0[2][0];
  ff[3][4] = xjm0[1][0]*xjm0[2][1] + xjm0[1][1]*xjm0[2][0];
  ff[5][0] = xjm0[0][0]*xjm0[0][2];
  ff[5][1] = xjm0[1][0]*xjm0[1][2];
  ff[5][2] = xjm0[2][0]*xjm0[2][2];
  ff[5][3] = xjm0[0][0]*xjm0[1][2] + xjm0[0][2]*xjm0[1][0];
  ff[5][5] = xjm0[0][0]*xjm0[2][2] + xjm0[0][2]*xjm0[2][0];
  ff[5][4] = xjm0[1][0]*xjm0[2][2] + xjm0[1][2]*xjm0[2][0];
  ff[4][0] = xjm0[0][1]*xjm0[0][2];
  ff[4][1] = xjm0[1][1]*xjm0[1][2];
  ff[4][2] = xjm0[2][1]*xjm0[2][2];
  ff[4][3] = xjm0[0][1]*xjm0[1][2] + xjm0[0][2]*xjm0[1][1];
  ff[4][5] = xjm0[0][1]*xjm0[2][2] + xjm0[0][2]*xjm0[2][1];
  ff[4][4] = xjm0[1][1]*xjm0[2][2] + xjm0[1][2]*xjm0[2][1];

  for (i=0; i<36; i++ ) fii[i] = 0.;
  for (i=0; i<36; i+=7) fii[i] = 1.;
  cc=0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) ffi[cc++]=ff[j][i];

  solveq (ffi,fii,&fdim,&fdim);

  cc=0;
  for (i=0; i<6; i++) for (j=0; j<6; j++) fi[j][i]=fii[cc++];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1t0 */

/*!----------------------------------------------------------------------
\brief evaluation of bop with extended deformations

<pre>                                                              al 06/02
This routine calcuates bop with extended deformations                                       |
for an 3D-hex-element.

</pre>
\param        **bop9  DOUBLE  (o) b-operator matrix modified   
\param    bn1[3][10]  DOUBLE  (i)   
\param      fi[6][6]  DOUBLE  (i)   
\param      disd1[9]  DOUBLE  (i)   
\param  ehdis[3][10]  DOUBLE  (i)   
\param          det0  DOUBLE  (i)   
\param          det1  DOUBLE  (i)   
\param            e1  DOUBLE  (i)   
\param            e2  DOUBLE  (i)   
\param            e3  DOUBLE  (i)   
\param           iel     INT  (i)   
\param            l1     INT  (i)   
\param            l3     INT  (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1bop9(
            DOUBLE          **bop9,
            DOUBLE      bn1[3][10],
            DOUBLE        fi[6][6],
            DOUBLE        disd1[9],
            DOUBLE    ehdis[3][10],
            DOUBLE            det0,
            DOUBLE            det1,
            DOUBLE              e1,
            DOUBLE              e2,
            DOUBLE              e3,
            INT                iel,
            INT                 l1,
            INT                 l3 
           )  
{
/*----------------------------------------------------------------------*/
INT i,j,k,m,cb,ce;
INT fdim=6;
DOUBLE rdet, dum;
DOUBLE eas[6][30],aux[9];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1bop9");
#endif
/*----------------------------------------------------------------------*/
      cb=iel*3;
      ce=(iel+l3)*3;
      for (i=0; i<6; i++) for (j=cb; j<ce; j++) bop9[i][j] = 0.;
/*----------------------------------------------------------------------*
 |   det1= det def grad at isoparametric origin (ne 1 for updated lagr) |
 |   det0= det of jacobian at isoparametric origin                      |
 |   det = det of jacobian at gauss point                               |
 *----------------------------------------------------------------------*/

      rdet=det0/det1;
      for (i=0;i<6;i++) for (j=0;j<l1;j++) eas[i][j] = 0. ;
      for (i=0;i<9;i++)   aux[i] = 0. ;
      for (i=0;i<9;i++) disd1[i] = 0. ;
/* enhanced strains for l3 = 3 (l1=9) */
      eas[3][3] = e1 * rdet;
      eas[5][6] = e1 * rdet;
      eas[4][7] = e2 * rdet;
/* enhanced strains for l3 >= 5 (l1 >= 15) */
      eas[3][17]= e1 * e3 * rdet;
      eas[5][20]= e1 * e2 * rdet;
      eas[4][19]= e1 * e2 * rdet;
      
    for (k=0; k<l3; k++) { 
      for (m=0; m<3; m++) {
        for (i=0; i<6; i++) {
          dum=0.;
          for (j=0; j<6; j++) dum += fi[j][i] * eas[j][m+k*3];
            bop9[i][(k+iel)*3 + m] += dum; }}}
      
      for (i=0; i<l3; i++) {
        for (j=0; j<6; j++) {
          for (k=0; k<3; k++) {
            aux[j] += bop9[j][k + (i+iel)*3] * ehdis[k][i];
      }}}
/* */
      disd1[4]=aux[3];
      disd1[6]=aux[4];
      disd1[8]=aux[5];

      for (i=0; i<6; i++) for (j=cb; j<ce; j++) bop9[i][j] = 0.;

/* enhanced strains for l3 = 3 (l1=9) */
      eas[0][0] = e1   * rdet;
      eas[1][4] = e2   * rdet;
      eas[2][8] = e3   * rdet;
      eas[3][1] = e2   * rdet;
      eas[5][2] = e3   * rdet;
      eas[4][5] = e3   * rdet;
/* enhanced strains for l3 = 5 (l1=15) */
      eas[0][ 9]= e1*e3* rdet;
      eas[0][10]= e1*e2* rdet;
      eas[1][11]= e1*e2* rdet;
      eas[1][12]= e2*e3* rdet;
      eas[2][13]= e1*e3* rdet;
      eas[2][14]= e2*e3* rdet;
/* enhanced strains for l3 = 7 (l1=21) */
      eas[3][15]= e2*e3* rdet;
      eas[3][17]= e1*e3* rdet;
      eas[5][16]= e2*e3* rdet;
      eas[5][20]= e1*e2* rdet;
      eas[4][18]= e1*e3* rdet;
      eas[4][19]= e1*e2* rdet;
/* enhanced strains for l3 = 10 (l1=30) */
        eas[3][21]= e2*e1* rdet;
        eas[5][23]= e1*e3* rdet;
        eas[4][22]= e2*e3* rdet;
        eas[0][24]= e1*e2*e3* rdet;
        eas[1][25]= e1*e2*e3* rdet;
        eas[2][26]= e1*e2*e3* rdet;
        eas[3][27]= e1*e2*e3* rdet;
        eas[5][29]= e1*e2*e3* rdet;
        eas[4][28]= e1*e2*e3* rdet;

    for (k=0; k<l3; k++) { 
      for (m=0; m<3; m++) {
        for (i=0; i<6; i++) {
          dum=0.;
          for (j=0; j<6; j++) dum += fi[j][i] * eas[j][m+k*3];
            bop9[i][(k+iel)*3 + m] += dum; }}}

      for (i=0; i<9; i++) aux[i] = 0.;

      for (i=0; i<l3; i+=3) {
        for (j=0; j<6; j++) {
          for (k=0; k<3; k++) {
            aux[j] += bop9[j][k + (i+iel)*3] * ehdis[k][i];
      }}}

      disd1[0]=aux[0];
      disd1[1]=aux[1];
      disd1[2]=aux[2];
      disd1[3]=aux[3]-disd1[4];
      disd1[5]=aux[4]-disd1[6];
      disd1[7]=aux[5]-disd1[8];

      for (k=0; k<l3; k++) {
       bn1[0][k] = bop9[0][0 + (k+iel)*3];
       bn1[1][k] = bop9[1][1 + (k+iel)*3];
       bn1[2][k] = bop9[2][2 + (k+iel)*3];
      }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1bop9 */

/*!----------------------------------------------------------------------
\brief include initial displacements to derivative operator

<pre>                                                              al 06/02
This routine includes initial displacements to derivative operator                                       |
for an 3D-hex-element.

</pre>
\param      **bop  DOUBLE  (o) b-operator  matrix   
\param  bn[3][10]  DOUBLE  (o) bn-operator matrix  
\param    disd[9]  DOUBLE  (i)   
\param        iel     INT  (i)   
\param         l3     INT  (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1bdish(
            DOUBLE       **bop,
            DOUBLE   bn[3][10],
            DOUBLE     disd[9],
            INT            iel,
            INT             l3 
           )  
{
/*----------------------------------------------------------------------*/
INT l,k,node_start;
DOUBLE rl11, rl12, rl13, rl21, rl22, rl23, rl31, rl32, rl33;
DOUBLE h1, h2, h3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1bdish");
#endif
/*----------------------------------------------------------------------*/
  rl11 = disd[0];
  rl12 = disd[3];
  rl13 = disd[7];
  rl21 = disd[4];
  rl22 = disd[1];
  rl23 = disd[5];
  rl31 = disd[8];
  rl32 = disd[6];
  rl33 = disd[2];

  for (k=iel; k<iel+l3; k++)
  {
    l = k - iel;
    node_start = k*3;

    h1 = bn[0][l];
    h2 = bn[1][l];
    h3 = bn[2][l];
    bop[0][node_start+0]  += rl11*h1;
    bop[0][node_start+1]  += rl21*h1;
    bop[0][node_start+2]  += rl31*h1;
    bop[1][node_start+0]  += rl12*h2;
    bop[1][node_start+1]  += rl22*h2;
    bop[1][node_start+2]  += rl32*h2;
    bop[2][node_start+0]  += rl13*h3;
    bop[2][node_start+1]  += rl23*h3;
    bop[2][node_start+2]  += rl33*h3;
    bop[3][node_start+0]  += rl11*h2 + rl12*h1;
    bop[3][node_start+1]  += rl21*h2 + rl22*h1;
    bop[3][node_start+2]  += rl31*h2 + rl32*h1;
    bop[4][node_start+0]  += rl12*h3 + rl13*h2;
    bop[4][node_start+1]  += rl22*h3 + rl23*h2;
    bop[4][node_start+2]  += rl32*h3 + rl33*h2;
    bop[5][node_start+0]  += rl11*h3 + rl13*h1;
    bop[5][node_start+1]  += rl21*h3 + rl23*h1;
    bop[5][node_start+2]  += rl31*h3 + rl33*h1;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1bdish */

/*!----------------------------------------------------------------------
\brief modify stiffness matrix and internal forces

<pre>                                                              al 06/02
This routine modifies stiffness matrix and internal forces 
for enhanced strain elements,  update displacement parameters 
for an 3D-hex-element.

</pre>
\param      *ele  ELEMENT  (i) element stiffness-matrix           
\param  **estif9   DOUBLE  (o) modified element stiffness-matrix 
\param   **estif   DOUBLE  (i) element residuals                
\param     *fieh   DOUBLE  (o)                                   
\param      *fie   DOUBLE  (o)   
\param        l1      INT  (i)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1rkefi(
            ELEMENT    *ele,
            DOUBLE **estif9,
            DOUBLE  **estif,
            DOUBLE    *fieh,
            DOUBLE     *fie,
            INT          l1
           )  
{
/*----------------------------------------------------------------------*/
INT i,j;
INT cc;
INT dim1 = 1;
INT dim24=24;
DOUBLE estif1[576];
DOUBLE fiehi[30];
DOUBLE saa[ 576];
DOUBLE sab[ 720];
DOUBLE sba[ 720];
DOUBLE sbai[720];
DOUBLE sbb[ 900];
DOUBLE sbbi[900];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1rkefi");
#endif
/*----------------------------------------------------------------------*/
  cc=0;
  for (j=0; j<l1; j++) {
    for (i=0; i<l1; i++) {
      sbb[cc++]=estif9[i+24][j+24];}}

  cc=0;
  for (j=0; j<l1; j++) {
    for (i=0; i<24; i++) {
      sab[cc++]=-estif9[i][j+24];
      }}

  cc=0;
  for (i=0; i<24; i++) {
    for (j=0; j<l1; j++) {
      sba[cc++]= estif9[i][j+24];
      }}


  cc=0;
  for (j=0; j<24; j++) {
    for (i=0; i<24; i++) {
      estif1[cc++]=estif9[i][j];
      }}
  cc=0;
  for (i=0; i<l1*l1; i++    ) sbbi[i]=0.;
  for (j=0; j<l1*l1; j+=l1+1) sbbi[j]=1.;

/*----------------------------------------------------------------------*
 |  determine inverse of stiffness matrix corresponding                 |
 |  to displacement parameters [he] with dim(l1 x l1)                   |
 *---------------------------------------------------------------------*/
  solveq (sbb,sbbi,&l1,&l1);

  mxmab(sbbi,sba,sbai,&l1,&l1,&dim24);
  mxmab(sbbi,fieh,fiehi,&l1,&l1,&dim1);

/*----------------------------------------------------------------------*
 |  save [ihe] [te] for update displacement parameters                  |
 |  ehdis = ehdis - [ihe] [te] [(edis - disl) - fieh ]                  |
 *---------------------------------------------------------------------*/

/*       call c1var ('numelm',numelm,'get')
       call c1wah1(wah,lwah,nch,nel,iel,sbai,fiehi,'put',l1)*/
  for (i=0; i<l1   ; i++) ele->e.c1->elewa[0].eas[0].hih[i] = fiehi[i];
  for (i=0; i<24*l1; i++) ele->e.c1->elewa[0].eas[0].hil[i] =  sbai[i];

/*----------------------------------------------------------------------*
 |  compute modified internal forces                                    |
 |  fie = fie - [te] [ihe] fieh                                         |
 *---------------------------------------------------------------------*/

  cc=0;
  for (j=0; j<l1; j++) {
    for (i=0; i<24; i++) {
      fie[i] += sab[cc++] * fiehi[j];
      }}
/*      do 35 i=1,3*iel
      do 35 j=1,l1
 35    fie(i) = fie(i) + sab(i,j)*fiehi(j)

/*----------------------------------------------------------------------*
 |  compute modified stiffness matrix                                   |
 |  estif = estif - [te] [ihe] [te]                                     |
 *---------------------------------------------------------------------*/

  mxmab(sab,sbai,saa,&dim24,&l1,&dim24);

  for (i=0; i<576; i++) estif1[i]+=saa[i];

  cc=0;
  for (j=0; j<24; j++) {
    for (i=0; i<24; i++) {
      estif[i][j]+=estif1[cc++];
      }}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1rkefi */

/*!----------------------------------------------------------------------
\brief evaluation of residual 

<pre>                                                              al 06/02
This routine evaluates residuals  
for an 3D-hex-element.

</pre>
\param      *F   DOUBLE  (i) force vector integral (stress-resultants)
\param     fac   DOUBLE  (i) multiplier for numerical integration     
\param  **bop9   DOUBLE  (i) b-operator matrix                        
\param     iel      INT  (i) number nodes of element                  
\param   *fieh   DOUBLE  (o) internal force vector                    
\param      l3      INT  (i) 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1res( 
           DOUBLE     *F,    /*  force vector integral (stress-resultants)*/
           DOUBLE    fac,    /*  multiplier for numerical integration     */
           DOUBLE **bop9,    /*  b-operator matrix                        */
           INT       iel,    /*  number nodes of element                  */
           DOUBLE  *fieh,    /*  internal force vector                    */
           INT        l3)    /*                                           */
{
/*----------------------------------------------------------------------*/
INT i,j,k;
INT ca,ce;
DOUBLE n11,n22,n33,n12,n23,n31;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1res");
#endif
/*---------------------------- set values of force vector components ---*/
  n11 = F[0]*fac;
  n22 = F[1]*fac;
  n33 = F[2]*fac;
  n12 = F[3]*fac;
  n23 = F[4]*fac;
  n31 = F[5]*fac;
/*---  loop over all nodal points ----- and over enhanced parameters ---*/
  ca = iel*3+2;
  ce = (iel+l3)*3;
  
  for (j=ca; j<ce; j+=3)
  {
    k=j-1;
    i=j-2;
    fieh[i-ca+2]+=  bop9[0][i]*n11 + bop9[1][i]*n22 + bop9[2][i]*n33 +
                    bop9[3][i]*n12 + bop9[4][i]*n23 + bop9[5][i]*n31;
    fieh[k-ca+2]+=  bop9[0][k]*n11 + bop9[1][k]*n22 + bop9[2][k]*n33 +
                    bop9[3][k]*n12 + bop9[4][k]*n23 + bop9[5][k]*n31;
    fieh[j-ca+2]+=  bop9[0][j]*n11 + bop9[1][j]*n22 + bop9[2][j]*n33 +
                    bop9[3][j]*n12 + bop9[4][j]*n23 + bop9[5][j]*n31;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1res */

/*!----------------------------------------------------------------------
\brief update of strain paramenters 

<pre>                                                              al 06/02
This routine updates strain paramenters  
for an 3D-hex-element.

</pre>
\param          *ele  ELEMENT  (i)
\param      edis[60]   DOUBLE  (i)
\param  ehdis[3][10]   DOUBLE  (i)
\param            l1      INT  (i)
\param            l3      INT  (i)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1upenh(  
             ELEMENT         *ele,
             DOUBLE      edis[60],  
             DOUBLE  ehdis[3][10],
             INT               l1,
             INT               l3
            )
{
/*----------------------------------------------------------------------*/
INT i,j,k;
INT cc;
DOUBLE disl[24];
DOUBLE ehdisl[30];
DOUBLE hih[30];
DOUBLE hil[30][24];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1upenh");
#endif
/*----------------------------------------------------------------------*/
    cc=0;
    for (i=0; i<24; i++) disl[i]=ele->e.c1->elewa[0].eas[0].disl[cc++];
    cc=0;
    for (i=0; i<l1; i++) ehdisl[i]=ele->e.c1->elewa[0].eas[0].ehdis[cc++];

    cc=0;
    for (i=0; i<l1; i++) hih[i]=ele->e.c1->elewa[0].eas[0].hih[cc++];
    cc=0;
    for (i=0; i<24; i++){ for (j=0; j<l1; j++){ 
                  hil[j][i]=ele->e.c1->elewa[0].eas[0].hil[cc++];}}

    /*----------------------------------- update enhanced parameters ----|
     |   ehdis = ehdisl - [ hi ] ( [ l ] ddisd + fieh )                  |
     |                                                                   |
     |   terms [ hi ]*(fieh)  and [- hi ]*[ l ] from last iteration      |
     |   stored in wah in vectors hil and hih resp. ( function:  c1rkefi)|
     *------------------------------------------------------------------*/    
   
    
    cc=0;
    for (i=0; i<l3; i++) {
      for (j=0; j<3; j++) { 
         ehdis[j][i] = ehdisl[cc] - hih[cc];
         cc++; }}

    cc=0;
    for (i=0; i<l3; i++) {
      for (j=0; j<3; j++) { 
        for (k=0; k<24; k++) { 
         ehdis[j][i] -= hil[cc][k]*(edis[k]-disl[k]);}
         cc++;}}
    /*---------------------------------------- put new values in wah ---*/
    cc=0;
    for (i=0; i<l3; i++) {
      for (j=0; j<3; j++) { 
               ele->e.c1->elewa[0].eas[0].ehdis[cc++] = ehdis[j][i];}}
    
    for (i=0; i<24; i++) ele->e.c1->elewa[0].eas[0].disl[i] = edis[i];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1upenh */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
