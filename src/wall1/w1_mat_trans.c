/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1mat_trans_up' which blows up the plane
       stress/strain conditions to 3D --> 3D material law
 contains the routine 'w1mat_trans_down' which condeses back the 3D
       stresses/strains/d to 2D conditions
 contains the routine 'w1_vec_switch' which changes to rows of a vector
 contains the routine 'w1_matrix_switch' which changes to rows and columns
       of a square matrix

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief blowing up plane stress/strain conditions

<pre>                                sh 7/02 modified by         sh 03/04
This routine performs the blowing up of plane stress/strain conditions
for a 3D-material formulation
</pre>
\param  DOUBLE       ym           (i)  youngs modulus
\param  DOUBLE       pv           (i)  poissons ratio
\param  ELEMENT     *ele          (i)  element array
\param  WALL_TYPE    wtype        (i)  wtype=plane_stress/plane_strain
\param  DOUBLE     **bop          (i)  B-Operator
\param  DOUBLE      *gop          (i)  for incomp. modes
\param  DOUBLE      *alpha        (i)  for incomp. modes
\param  INT          ip           (i)  actual integration point
\param  DOUBLE       strain[4]    (o)  strains (extended from) [11,22,12,33]

\warning  NOTE: everything in here is with sorting [11,22,12,33]
\return void
\sa calling: ---; called by: w1_mat_plast_mises_3D()   [w1_mat_plast_mises3D.c]
                             w1_mat_plast_epc3D()      [w1_mat_plast_epc3D.c]

*----------------------------------------------------------------------*/
void w1mat_trans_up (DOUBLE     ym,
                     DOUBLE     pv,
                     ELEMENT   *ele,
                     WALL_TYPE  wtype,
                     DOUBLE   **bop,
                     DOUBLE    *gop,
                     DOUBLE    *alpha,
                     INT        ip,
                     DOUBLE     strain[4])
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE disd[5];
DOUBLE sigi[4];     /*stress from last iteration step (not condensed*/
DOUBLE epsi[4];     /*strains from last iteration step (with e33)*/
DOUBLE di[4];       /*components d41,d42,d43,d44 of th const. tensor from
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


/*extend the strain vector by eps_zz if plane_stress */
if(wtype==plane_stress)
{
  for (i=0; i<4; i++)
  {
    sigi[i]   = ele->e.w1->elewa[0].ipwa[ip].sigi[i];
    epsi[i]   = ele->e.w1->elewa[0].ipwa[ip].epsi[i];
    di[  i]   = ele->e.w1->elewa[0].ipwa[ip].di[  i];
  }
  if (fabs(di[3]) - 0.0001 < 0.)  /*get elastic constitutive matrix DI*/
  {
    w1iwadi (ym, pv, di);
  }
  w1de33 (sigi,epsi,di,strain);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1mat_trans_up */


/*!----------------------------------------------------------------------
\brief kondense 3D conditions to plane stress/strain conditions

<pre>                                sh 7/02 modified by         sh 03/04
This routine condenses stresses and D from 3D to
plane stress/strain conditions
</pre>
\param  DOUBLE     **d           (i/o) to be condensed
\param  ELEMENT     *ele          (i)  element array
\param  WALL_TYPE    wtype        (i)  wtype=plane_stress/plane_strain
\param  INT          ip           (i)  actual integration point
\param  INT          yipc         (i)  flag for storing in Ele-WA
\param  DOUBLE      *stressc      (o)  condensed stresses [11,22,12,33]
\param  DOUBLE      *sig          (i)  sig from last update [11,22,12,33]
\param  DOUBLE      *eps          (i)  eps from last update [11,22,12,33]
\param  DOUBLE      *stress       (i)  stress [11,22,12,33] to be condensed
\param  DOUBLE      *strain       (i)  strain  [11,22,12,33]
\param  DOUBLE      *qn           (i)  qn [11,22,12,33] to be condensed

\warning NOTE: everything in here is with sorting [11,22,12,33]
\return void
\sa calling: ---; called by: w1_mat_plast_mises_3D()   [w1_mat_plast_mises3D.c]
                             w1_mat_plast_epc3D()      [w1_mat_plast_epc3D.c]

*----------------------------------------------------------------------*/
void w1mat_trans_down (DOUBLE   **d,
                       ELEMENT   *ele,
                       WALL_TYPE  wtype,
                       INT        ip,
                       INT        yipc,
                       DOUBLE    *stressc,
                       DOUBLE    *sig,
                       DOUBLE    *eps,
                       DOUBLE    *stress,
                       DOUBLE    *strain,
                       DOUBLE    *qn)
{
/*----------------------------------------------------------------------*/
INT i;
DOUBLE tau[4];
DOUBLE tauc[4];
DOUBLE qnc[4];
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
/*----------------------------------------------------------------------*/
  for (i=0; i<4; i++) tau[i]  = stress[i] - qn[i];
  for (i=0; i<4; i++) tauc[i] = tau[i] ;
  for (i=0; i<4; i++) qnc[i]  = qn[i];

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
    for (i=0; i<4; i++)
    {
      strain[i] = eps[i];
    }
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
     ele->e.w1->elewa[0].ipwa[ip].sigi[i] = stress[i];
     ele->e.w1->elewa[0].ipwa[ip].epsi[i] = strain[i];
     ele->e.w1->elewa[0].ipwa[ip].di[  i] = d[i][3];
    }
    w1concep (d);
  }
}

for (i=0 ; i<4; i++)
{
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
void w1_vec_switch(DOUBLE *vec,      /* vector do be modified           */
                   INT a,            /* row to be changed to b          */
                   INT b)            /* row to be changed to a          */
{
DOUBLE help;
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
void w1_matrix_switch(DOUBLE **mat,   /* matrix do be modified          */
                      INT a,          /* row & colum to be changed to b */
                      INT b,          /* row & colum to be changed to a */
                      INT l)          /* length of row/column of matrix */
{
DOUBLE help;
INT i;
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
