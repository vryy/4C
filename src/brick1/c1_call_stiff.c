/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_keku' which calclate usual stiffness matrix
       for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief usual stiffness matrix total lagrangian formulation

<pre>                                                              al 06/02
This routine calcuates the usual stiffness matrix for an 3D-hex-element.

</pre>
\param    **s    DOUBLE  (o)   element stiffness matrix 
\param   **bs    DOUBLE  (i)   derivative operator  
\param    **d    DOUBLE  (i)   constitutive matrix
\param    fac    DOUBLE  (i)   integration factor
\param     nd    INT     (i)   total number degrees of freedom of element
\param   neps    INT     (i)   actual number of strain components 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_keku(double  **s, 
             double  **bs, 
             double  **d, 
             double    fac, 
             int       nd,
             int       neps)
{
int            i, j, k, l, m;
double         dum;
double         db[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_keku");
#endif
/*----------------------------------------------------------------------*/
   for (j=0; j<nd; j++)
   {
     for (k=0; k<neps; k++)
     {
      db[k] = 0.0 ;                                                              
       for (l=0; l<neps; l++)
       {
       db[k] = db[k] + d[k][l]*bs[l][j]*fac ;
       }
     }
     for (i=0; i<nd; i++)
     {
       dum = 0.0 ;                                                                
       for (m=0; m<neps; m++)
       {
        dum = dum + bs[m][i]*db[m] ;
       }
        s[i][j] = s[i][j] + dum ;
     }
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_keku */


/*!----------------------------------------------------------------------
\brief geometric stiffness matrix (initial stress)

<pre>                                                              al 06/02
This routine calcuates the usual stiffness matrix for an 3D-hex-element.

</pre>
\param    **s    DOUBLE  (o)   element stiffness matrix 
\param   **bn    DOUBLE  (i)   b-operator matrix  
\param     *f    DOUBLE  (i)   force vector integral (stress-resultants)
\param    fac    DOUBLE  (i)   integration factor
\param    iel    INT     (i)   number of nodes at actual element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1vkg(  double  **s,   /* element stiffness-matrix                 */
             double  **bn,  /* b-operator matrix                        */
             double   *f,   /* force vector integral (stress-resultants)*/
             double    fac, /* multiplier for numerical integration     */
             int       iel) /* number of nodes at actual element         */
{
int            i, j, k, l, m, n, i1, i2, i3, sc1, sc2;
int            node_start1, node_start2;
double         g11, g12, g13, b11, b21, b31, n11, n22, n33, n12, n23, n31;
double         nb[9][3];
double       dum1, dum2, dum3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1vkg");
#endif
/*---------------------------- set values of force vector components ---*/
  n11=f[0]*fac;
  n22=f[1]*fac;
  n33=f[2]*fac;
  n12=f[3]*fac;
  n23=f[4]*fac;
  n31=f[5]*fac;
/*--------------------------------------- loop over all nodal points ---*/
  sc1 = 0;
  sc2 = 0;
  for (j=0; j<iel; j++)
  {
    b11 = bn[0][j];
    b21 = bn[1][j];
    b31 = bn[2][j];
/*------------------------------------ initialize matrix  n*b  (9*3) ---*/
    for (m=0; m<9; m++) for (n=0; n<3; n++) nb[m][n] = 0.;

    nb[0][0]=n11*b11 + n12*b21 + n31*b31;
    nb[1][0]=n12*b11 + n22*b21 + n23*b31;
    nb[2][0]=n31*b11 + n23*b21 + n33*b31;
    nb[3][1]=nb[0][0];
    nb[4][1]=nb[1][0];
    nb[5][1]=nb[2][0];
    nb[6][2]=nb[0][0];
    nb[7][2]=nb[1][0];
    nb[8][2]=nb[2][0];
/*------------------------------------------- transposed  b - matrix ---*/
    for (i=0; i<iel; i++)
    {
      g11 = bn[0][i];
      g12 = bn[1][i];
      g13 = bn[2][i];
/*------------------------------------- loop over degrees of freedom ---*/
      for (l=0; l<3; l++)
      {
        sc1 = j*3+l;
        sc2 = i*3;
        s[sc1][sc2+0] +=  g11*nb[0][l] + g12*nb[1][l] + g13*nb[2][l];
        s[sc1][sc2+1] +=  g11*nb[3][l] + g12*nb[4][l] + g13*nb[5][l];
        s[sc1][sc2+2] +=  g11*nb[6][l] + g12*nb[7][l] + g13*nb[8][l];
/*-------------------------- end of loop over all degrees of freedom ---*/
      }
/*-------------------------- end of loop over remaining nodal points ---*/
    }
/*-------------------------------- end of loop over all nodal points ---*/
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1vkg */
/*----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief reorder stiffness matrix for 'gid'  hex20

<pre>                                                              al 06/02
This routine reorders stiffness matrix for 'gid'  hex20 for an 3D-hex-element.

</pre>
\param   **si    DOUBLE  (i)  element stiffness-matrix - convent.
\param   **so    DOUBLE  (o)  element stiffness-matrix - gid.    

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1kgid( double  **si,  /* element stiffness-matrix - convent.      */
             double  **so)  /* element stiffness-matrix - gid.          */
{
/*----------------------------------------------------------------------*/
  int i,j;
  int reor[20] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1kgid");
#endif
/*----------------------------------------------------------------------*/
   for (i=0; i<20; i++) {
      for (j=0; j<20; j++){
        
        so[ reor[i]*3+0 ][ reor[j]*3+0 ]=si[i*3+0][j*3+0];
        so[ reor[i]*3+0 ][ reor[j]*3+1 ]=si[i*3+0][j*3+1];
        so[ reor[i]*3+0 ][ reor[j]*3+2 ]=si[i*3+0][j*3+2];
        so[ reor[i]*3+1 ][ reor[j]*3+0 ]=si[i*3+1][j*3+0];
        so[ reor[i]*3+1 ][ reor[j]*3+1 ]=si[i*3+1][j*3+1];
        so[ reor[i]*3+1 ][ reor[j]*3+2 ]=si[i*3+1][j*3+2];
        so[ reor[i]*3+2 ][ reor[j]*3+0 ]=si[i*3+2][j*3+0];
        so[ reor[i]*3+2 ][ reor[j]*3+1 ]=si[i*3+2][j*3+1];
        so[ reor[i]*3+2 ][ reor[j]*3+2 ]=si[i*3+2][j*3+2];
        
        /*
        so[ reor[i]*3+0 ][ reor[j]*3+0 ]=si[i*3+0][j*3+0];
        so[ reor[i]*3+1 ][ reor[j]*3+1 ]=si[i*3+1][j*3+1];
        so[ reor[i]*3+2 ][ reor[j]*3+2 ]=si[i*3+2][j*3+2];
        */
        }}
/*        so[i][j]=si[reor[i]][reor[j]];}}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1kgid */
/*----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief reorder element force vector for 'gid'  hex20

<pre>                                                              al 06/02
This routine reorders force vector for 'gid'  hex20 for an 3D-hex-element.

</pre>
\param   **fi    DOUBLE  (i)  element force vector - convent.
\param   **fo    DOUBLE  (o)  element force vector - gid.    

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1fgid( double  *fi,  /* element force vector - convent.       */
             double  *fo)  /* element force vector - gid.           */
{
/*----------------------------------------------------------------------*/
  int i,j;
  int reor[20] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1fgid");
#endif
/*----------------------------------------------------------------------*/
   for (i=0; i<20; i++)
   {
    fo[ reor[i]*3+0 ]=fi[i*3+0];
    fo[ reor[i]*3+1 ]=fi[i*3+1];
    fo[ reor[i]*3+2 ]=fi[i*3+2];
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1fgid */
/*----------------------------------------------------------------------*/
#endif
 /*! @} (documentation module close)*/





