/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1mat_trans_up' which blows up the plane
       stress/strain conditions to 3D --> 3D material law
 contains the routine 'w1mat_trans_down' which condeses back the 3D
       stresses/strains/d to 2D conditions
 contains the routine 'w1_vec_switch' which changes to rows of a vector
 contains the routine 'w1_matrix_switch' which changes to rows and columns
       of a square matrix
 
*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------|
|      topic: blowing up plane stress/strain conditions           sh 7/02|
|             to 3D --> 3D-Material Law                                  |
|-----------------------------------------------------------------------*/
void w1mat_trans_up (double     ym,
                     double     pv,
                     ELEMENT   *ele,                                        
                     WALL_TYPE  wtype,
                     double   **bop,
                     double    *gop,
                     double    *alpha,
                     int        ip,
                     double    *stress,    /*actuel stress condensed        */
                     double    *stress3D,  /*actuel stress condensed [6]    */
                     double    *strain3D,  /*strains to be calculated [6]   */
                     double    *sig3D,     /*stresses from last update [6]  */ 
                     double    *eps3D,     /*strains from last update [6]   */
                     double    *qn3D,      /*backstress vektor [6]          */
                     int        newval)    /*controls evaluation of new stresses*/           
{
/*----------------------------------------------------------------------*/
int i,j;
double disd[5];
double sig[4];      /*stresses from last update -> WA*/
double eps[4];      /*strains from last update -> WA*/
double strain[4];   /*actual strains from displacements*/
double qn[4];       /*backstress vektor from last update -> WA*/
double sigi[4];     /*stress from last iteration step (not condensed*/
double epsi[4];     /*strains from last iteration step (with e33)*/
double di[4];       /*components d41,d42,d43,d44 of th const. tensor from
                      last iteration step*/
WALL_TYPE local_wtype;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1mat_trans_up");
#endif
/*------------ original global elastic matrix for current point -> D ---*/
  /* look at horst's work (pp.30): in case of plane stress ->
     switch to plane strain and condense stress and material 
     tensor at the end */
  local_wtype=plane_strain;
/*--------------------------------- compute displacement derivatives ---*/        
  w1_disd (ele,bop,gop,alpha,wtype,disd) ;                  
/*------------------------------------- get actual strains -> strain ---*/
  w1_eps (disd,local_wtype,strain);
/*----------------------------- get old values -> sig, eps,epstn,yip ---*/
  for (i=0; i<4; i++)
  {
    sig[i] = ele->e.w1->elewa[0].ipwa[ip].sig[i];
    eps[i] = ele->e.w1->elewa[0].ipwa[ip].eps[i];
    qn[i]  = ele->e.w1->elewa[0].ipwa[ip].qn[i];
  }

  if(newval==1)
  {
    for (i=0; i<4; i++)  stress[i] = sig[i];
    goto end;
  }

/*extend the strain vector by eps_zz if plane_stress */  
if(wtype==plane_stress)      
{  
  for (i=0; i<4; i++)
  {
    sigi[i]   = ele->e.w1->elewa[0].ipwa[ip].sigi[i];
    epsi[i]   = ele->e.w1->elewa[0].ipwa[ip].epsi[i];
    di[  i]   = ele->e.w1->elewa[0].ipwa[ip].di[  i];
  }
  if (fabs(di[3]) - 0.0001 < 0.)
  {
    w1iwadi (ym, pv, di);
  }
  w1de33 (sigi,epsi,di,strain);
} 

/*-- change row 3<->4 and enlarge to 6 entries -> 3D ---------*/
/*-- Sort: 11[0],22[1],33[2],12[3],23[4],13[5]  3D  (brick) --*/
/*--       11[0],22[1],12[2],33[3]                  (wall)  --*/

for (i=0; i<4; i++)
{
  sig3D[i]    = sig[i];
  strain3D[i] = strain[i];
  eps3D[i]    = eps[i];
  qn3D[i]     = qn[i];
}

for (i=4; i<6; i++)
{
  sig3D[i]    = 0.;
  strain3D[i] = 0.;
  eps3D[i]    = 0.;
  qn3D[i]     = 0.;
}

w1_vec_switch(sig3D,2,3);      
w1_vec_switch(strain3D,2,3);      
w1_vec_switch(eps3D,2,3);      
w1_vec_switch(qn3D,2,3);      


/*sig3D[0] = sig[0];    
sig3D[1] = sig[1];
sig3D[2] = sig[3];
sig3D[3] = sig[2];
sig3D[4] = 0.0;
sig3D[5] = 0.0;

strain3D[0] = strain[0];    
strain3D[1] = strain[1];
strain3D[2] = strain[3];
strain3D[3] = strain[2];
strain3D[4] = 0.0;
strain3D[5] = 0.0;

eps3D[0] = eps[0];    
eps3D[1] = eps[1];
eps3D[2] = eps[3];
eps3D[3] = eps[2];
eps3D[4] = 0.0;
eps3D[5] = 0.0;

qn3D[0] = qn[0];    
qn3D[1] = qn[1];
qn3D[2] = qn[3];
qn3D[3] = qn[2];
qn3D[4] = 0.0;
qn3D[5] = 0.0;*/

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1mat_trans_up */


/*-----------------------------------------------------------------------|
|      topic: kondense 3D conditions                              sh 7/02|
|             to plane stress/strain conditions                          |
|-----------------------------------------------------------------------*/
void w1mat_trans_down (double   **d, /*current material matrix 3D -> 2D */
                       ELEMENT   *ele,
                       WALL_TYPE  wtype,
                       int        ip,
                       int        yipc,
                       double    *stressc, /*condensed*/
                       double    *sig,
                       double    *eps,
                       double    *stress,
                       double    *strain,  /*actual strain [4]          */
                       double    *qn)
{
/*----------------------------------------------------------------------*/
int i,j;
double tau[4];
double tauc[4];
double qnc[4];
double di[4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1mat_trans_down");
#endif
/*----------------------------------------------------------------------*/
/*Initialize*/
for (i=0; i<4; i++) 
{
  tau[i]    = 0.0;
  tauc[i]   = 0.0;
  qnc[i]    = 0.0;
}
/*--------Sorting from [11,22,33,12,23,13] -> [11,22,12,33 | 23,13]  ---*/
/*------------    switch columns and row position 3<->4           ------*/
  for (i=0; i<6; i++) tau[i]  = stress[i] - qn[i];

  w1_vec_switch(tau,2,3);      
  w1_vec_switch(stress,2,3);      
  w1_vec_switch(strain,2,3);      
  w1_vec_switch(qn,2,3);      
  w1_vec_switch(sig,2,3);      
  w1_vec_switch(eps,2,3);      
  w1_matrix_switch(d,2,3,6);      

  for (i=0; i<6; i++) tauc[i] = tau[i] ;
  for (i=0; i<6; i++) qnc[i]  = qn[i];
/******************************************/
if (yipc>0)     /*stresses are available from last update*/
{
  if(wtype==plane_stress)
  {
    for (i=0; i<4; i++)
    {
     ele->e.w1->elewa[0].ipwa[ip].sigi[i] = sig[i];
     ele->e.w1->elewa[0].ipwa[ip].epsi[i] = eps[i];
     ele->e.w1->elewa[0].ipwa[ip].di[  i] = d[i][3];
     strain[i] = eps[i];
    }
    w1concep (d);
  }  
}
else /*(yipc < 0)*/
{
  if(wtype==plane_stress)
  {
    w1consig (d,tau,tauc);
    w1consig (d,qn,qnc);
    for (i=0; i<4; i++)
    {
     ele->e.w1->elewa[0].ipwa[ip].sigi[i] = tau[i];
     ele->e.w1->elewa[0].ipwa[ip].epsi[i] = strain[i];
     ele->e.w1->elewa[0].ipwa[ip].di[  i] = d[i][3];
    }
    w1concep (d);
  }
}

for (i=0 ; i<4; i++)
{
  stress[ i] = tau[ i] + qn[ i];
  stressc[i] = tauc[i] + qnc[i];
}  

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of w1mat_trans_down */




/*-----------------------------------------------------------------------*
|     changes 2 rows of a vector                               sh 8/02   |
|  vec[x,x,a,x,b,x,...] -> vec[x,x,b,x,a,x,...]                          |
*-----------------------------------------------------------------------*/
void w1_vec_switch(double *vec,      /* vector do be modified           */
                   int a,            /* row to be changed to b          */
                   int b)            /* row to be changed to a          */
{
double help;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_vec_switch");
#endif
/*----------------------------------------------------------------------*/
help   = vec[a];
vec[a] = vec[b];
vec[b] = help;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_vec_switch */



/*-----------------------------------------------------------------------*
|     changes 2 rows & columns of a square matrix              sh 8/02   |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
|     [a,a,a,c,a]        [b,b,b,c,b]                                     |
|  mat[-,a,-,b,-] ->  mat[-,b,-,a,-]                                     |
|     [b,c,b,b,b]        [a,c,a,a,a]                                     |
|     [-,a,-,b,-]        [-,b,-,a,-]                                     |
*-----------------------------------------------------------------------*/
void w1_matrix_switch(double **mat,   /* matrix do be modified          */
                      int a,          /* row & colum to be changed to b */
                      int b,          /* row & colum to be changed to a */
                      int l)          /* length of row/column of matrix */
{
double help;
int i,j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_matrix_switch");
#endif
/*----------------------------------------------------------------------*/
  for (i=0; i<l; i++)   /*rows*/
  {
    help      = mat[i][a];
    mat[i][a] = mat[i][b];
    mat[i][b] = help;
  }

  for (i=0; i<l; i++)   /*colums*/
  {
    help      = mat[a][i];
    mat[a][i] = mat[b][i];
    mat[b][i] = help;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_matrix_switch */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
