/*======================================================================*/
/*!
\file
\brief Tensor manipulations

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 01/07
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Identity 2-tensor in 3 dimensions.

\param  at     DOUBLE[][]    (i)  input tensor
\return void

\author bborn
\date 01/07
*/
void so3_tns3_id(DOUBLE it[3][3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_tns3_it");
#endif

  memset(it, 0, sizeof(it));
  for (i=0; i<3; i++)
  {
    it[i][i] = 1.0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Trace of a 2-tensor in 3 dimensions.

About:
  Trace operator of a tensor A is defined as
    tr(A) = SUM(A_ii,i=1,3)

\param  at     DOUBLE[][]    (i)  input tensor
\param  tr     DOUBLE*       (o)  trace
\return void

\author bborn
\date 01/07
*/
void so3_tns3_tr(DOUBLE at[3][3],
                 DOUBLE *tr)
{
  INT i;  /* index */
  DOUBLE trace;  /* auxiliar trace */

#ifdef DEBUG
  dstrc_enter("so3_tns3_tr");
#endif

  trace = 0.0;
  for (i=0; i<3; i++)
  {
    trace += at[i][i];
  }
  *tr = trace;

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Determinant of a 2-tensor in 3 dimensions.

\param  at     DOUBLE[][]    (i)  input tensor
\param  det    DOUBLE*       (o)  determinant
\return void

\author bborn
\date 01/07
*/
void so3_tns3_det(DOUBLE at[3][3],
                  DOUBLE *det)
{

#ifdef DEBUG
  dstrc_enter("so3_tns3_det");
#endif

  /* determionant of 3x3-matrix by Sarrus' rule */
  *det = at[0][0] * at[1][1] * at[2][2]
       + at[0][1] * at[1][2] * at[2][0]
       + at[0][2] * at[1][0] * at[2][1]
       - at[0][2] * at[1][1] * at[2][0]
       - at[0][0] * at[1][2] * at[2][1]
       - at[0][1] * at[1][0] * at[2][2];

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Invariants of a 2-tensor A in 3dim

1st invariant
    I_A = tr(A)
2nd invariant
    II_A = 1/2*(tr(A)^2 - tr(A^2))
3rd invariant
    III_A = det(A)

\param  at     DOUBLE[][]    (i)  tensor A
\param  ai     DOUBLE*       (o)  I_A
\param  aii    DOUBLE*       (o)  II_A
\param  aiii   DOUBLE*       (o)  III_A
\return void

\author bborn
\date 01/07
*/
void so3_tns3_inva(DOUBLE at[3][3],
                   DOUBLE *ai,
                   DOUBLE *aii,
                   DOUBLE *aiii)
{
  INT i;  /* index */
  DOUBLE aix, aiix;  /* auxiliar inva. */

#ifdef DEBUG
  dstrc_enter("so3_tns3_inva");
#endif

  /* 1st invariant */
  if (ai != NULL)
  {
    so3_tns3_tr(at, ai);
  }

  /* 2nd invariant */
  if (aii != NULL)
  {
    if (ai != NULL)
    {
      aiix = (*ai)*(*ai);
    }
    else
    {
      so3_tns3_tr(at, &aix);
      aiix = aix*aix;
    }
    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
      {
        aiix -= at[i][j]*at[j][i];
      }
    }
    aiix *= 0.5;
    *aii = aiix;
  }

  /* 3rd invariant */
  if (aiii != NULL)
  {
    so3_tns3_det(at, aiii);
  }


#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Dot product of two 2-tensors in 3 dimensions

  We talk about the trivial
    C = A . B

\param  at     DOUBLE[][]    (i)  input tensor
\param  bt     DOUBLE[][]    (i)  input tensor
\param  ct     DOUBLE[][]    (o)  product tensor
\return void

\author bborn
\date 01/07
*/
void so3_tns3_dotprod(DOUBLE at[3][3],
                      DOUBLE bt[3][3],
                      DOUBLE ct[3][3])
{
  INT i, j, k;  /* indices */
  DOUBLE ctjk;  /* auxiliar resulting component */

#ifdef DEBUG
  dstrc_enter("so3_tns3_dotprod");
#endif

  for (j=0; j<3; j++)
  {
    for (k=0; k<3; k++)
    {
      ctjk = 0.0;
      for (i=0; i<3; i++)
      {
        ctjk += at[j][i] * bt[i][k];
      }
      ct[j][k] = ctjk;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Dot product of two 2-tensors in 3 dimensions, 
       left input tensor is transposed

We talk about the trivial
    C = A^T . B

\param  at     DOUBLE[][]    (i)  input tensor
\param  bt     DOUBLE[][]    (i)  input tensor
\param  ct     DOUBLE[][]    (o)  product tensor
\return void

\author bborn
\date 01/07
*/
void so3_tns3_dotprod_tl(DOUBLE at[3][3],
                         DOUBLE bt[3][3],
                         DOUBLE ct[3][3])
{
  INT i, j, k;  /* indices */
  DOUBLE ctjk;  /* auxiliar resulting component */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns3_dotprod_tl");
#endif

  for (j=0; j<3; j++)
  {
    for (k=0; k<3; k++)
    {
      ctjk = 0.0;
      for (i=0; i<3; i++)
      {
        ctjk += at[i][j] * bt[i][k];
      }
      ct[j][k] = ctjk;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Dot product of two 2-tensors in 3 dimensions, 
       right input tensor is transposed

We talk about the trivial
    C = A . B^T

\param  at     DOUBLE[][]    (i)  input tensor A
\param  bt     DOUBLE[][]    (i)  input tensor B
\param  ct     DOUBLE[][]    (o)  product tensor C
\return void

\author bborn
\date 01/07
*/
void so3_tns3_dotprod_tr(DOUBLE at[3][3],
                         DOUBLE bt[3][3],
                         DOUBLE ct[3][3])
{
  INT i, j, k;  /* indices */
  DOUBLE ctjk;  /* auxiliar resulting component */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns3_dotprod_tr");
#endif

  for (j=0; j<3; j++)
  {
    for (k=0; k<3; k++)
    {
      ctjk = 0.0;
      for (i=0; i<3; i++)
      {
        ctjk += at[j][i] * bt[k][i];
      }
      ct[j][k] = ctjk;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Polar decomposition of a 2-tensor in 3 dimensions.

About:
  Polar decomposition is often applied to the material deformation
tensor F, i.e.
     F = R . U = v . R
in which R is the rotation matrix (two-point tensor), U the material
stretch tensor (refers to the undeformed configuration) and v the
spatial stretch tensor (refers to the deformed configuration).
  This polar decomposition is also applied to the isoparametric
Jacobian tensor (mapping quantities in parameter space to material
configuration). However, the local variables are denoted for the
case of a deformation tensor.

References:
[1] A. Hoger & D.E. Carlson, "Determination of the stretch and
        rotation in the polar decomposition of the deformation 
        gradient", Quart. Appl. Math., 42(2):113-117, 1984.
[2] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
        esp. Section 2.6
[3] I.N. Bronstein & K.A. Semendjajew, "Taschenbuch der
        Mathematik", Teubner, 25.ed, 1991.
       

\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void

\author bborn
\date 01/07
*/
void so3_tns3_plrdcmp(DOUBLE ft[3][3],  /* input tensor */
                      DOUBLE rt[3][3],  /* rotation matrix */
                      DOUBLE ut[3][3],  /* right stretch tensor */
                      DOUBLE vt[3][3])  /* left stretch tensor */
{
  DOUBLE it[3][3];  /* identity 2-tensor I */
  DOUBLE ct[3][3];  /* right Cauchy-Green tensor C */
  DOUBLE ci, cii, ciii;  /* invarients of C */
  DOUBLE c2t[3][3];  /* squared Cauchy-Green tensor */
  DOUBLE c2i;  /* invariants of C^2 */
  DOUBLE ui, uii, uiii;  /* invariants of U */
  DOUBLE invut[3][3];  /* inverted material stretch tensor */
  INT i, j, k;  /* indices */
  DOUBLE xi, eta, zeta;
  DOUBLE denom, pc2, pc, pi;  /* coefficients */
  DOUBLE p, q, r;  /* coefficients in quartic polynomial of tr(U) */
  DOUBLE pp, qq;
  DOUBLE disc, rho, phi, rhort;
  DOUBLE x1, x2, x3;  /* roots of normalised cubic resolvent */
  DOUBLE z1, z2, z3;  /* roots of cubic resolvent */
  DOUBLE z1rt, z2rt, z3rt;  /* radicals of roots of cubic resolvent */
  DOUBLE rtx[3][3];  /* auxiliar rotation matrix in case of solely
                      * interested in vt */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns_plrdcmp");
#endif

  /* set identity tensor */
  so3_tns3_id(it);

  /* right Cauchy-Green deformation tensor */
  /* C = F^T . F */
  so3_tns3_dotprod_tl(ft, ft, ct);

  /* squared right Cauchy-Green deformation tensor */
  /* C^2 = C . C */
  so3_tns3_dotprod(ct, ct, c2t);

  /* invariants of right Cauchy-Green tensor */
  /* 1st invariant: I_C = tr(C) */
  so3_tns3_tr(ct, &ci);
  /* 2nd invariant: II_C = 1/2 ( tr(C)^2 - tr(C^2) ) */
  so3_tns3_tr(c2t, &c2i);
  cii = 0.5*(ci*ci - c2i);
  /* 3rd invariant: III_C = det(C) */
  so3_tns3_det(ct, &ciii);


  /*--------------------------------------------------------------------*/
  /* determination of I_U acc. to [1], 
   * ==> BUT FAILS ==> next block working alternative */
#if 0
  /* auxiliar variables to get trace of material stretch tensor U */
  xi = 32.0*(2.0*ci*ci*ci - 9.0*ci*cii + 27.0*ciii)/27.0;
  eta = 1024.0*(4.0*cii*cii*cii - ci*ci*cii*cii + 4.0*ci*ci*ci*ciii 
                - 18.0*ci*cii*ciii + 27.0*ciii*ciii)/27.0;
  zeta = -2.0*ci/3.0 + pow(xi+sqrt(eta), 1.0/3.0) 
    + pow(xi-sqrt(eta), 1.0/3.0);
  
  /* invariants of material stretch tensor U */
  /* 1st invariant: I_U = tr(U) */
  if (zeta == -2.0*ci)
  {
    ui = sqrt(ci + 2.0*sqrt(cii));
  }
  else
  {
    ui = 0.5*(sqrt(2.0*ci+zeta) 
              + sqrt(2.0*ci-zeta+16.0*sqrt(ciii)/sqrt(2.0*ci+zeta)));
  }
#endif

#if 1
  /*--------------------------------------------------------------------*/
  /* 1st invariant of the material stretch tensor U */
  /* Summary:
   *     The 1st invariant I_U depends non-linearly 
   *     on I_C, II_C and III_C.
   */
  /* Derivation:
   * This relation is found by Cayley-Hamilton's theorem applied to
   * the material stretch tensor U
   *     U^3 - I_U*U^2 + II_U*U - III_U*I = 0                        (1)
   * in which U^2 = U . U and U^3 = U . U . U and I is the indentity
   * tensor in 3dim.
   * The trace of this equation is
   *     tr(U^3) - I_U*tr(U^2) + II_U*tr(U) - 3*III_U = 0
   * Equivalently with tr(U^2)=tr(C)=I_C, tr(U) = I_U it is obtained
   *     tr(U^3) - I_U*I_C + II_U*I_U - 3*III_U = 0                  (2)
   * The problem is the unknown tr(U^3). II_U and III_U can be referred
   * back to I_U and the invariants if C. 
   *     II_U = 1/2*(tr(U)^2 - tr(U^2)) = 1/2*(I_U^2 - I_C)          (3)
   *     III_U = det(U) = sqrt(det(U)*det(U)) 
   *           = sqrt(det(U . U)) = sqrt(det(C)) = sqrt(III_C)       (4)
   *
   * We can resolve the problem by
   * considering Eq.(*) multiplied with U, ie
   *     U^4 - I_U*U^3 + II_U*U^2 - III_U*U = 0
   * The trace of this matrix equation reveals
   *     tr(U^4) - I_U*tr(U^3) + II_U*tr(U^2) - III_U*tr(U) = 0
   * or with tr(U^4)=tr(C^2)=I_C^2 - 2*II_C, 
   *     I_C^2 - 2*II_C - I_U*tr(U^3) + II_U*I_C - III_U*I_U = 0     (5)
   *
   * Eqs (2) and (5) can be used to eliminate tr(U^3), and we end up at
   *     I_U^4 - 2*I_C*I_U^2 - 8*sqrt(III_C)*I_U 
   *                                        + I_C^2 - 4*II_C = 0     (6)
   * in which advantage is taken of Eqs (3) and (4).
   * This quartic polynomial of I_U needs to be solved.
   * It is not done by the formulas presented in [1], because they seem
   * the fail even for the simple test case in which the deformation
   * gradient is a diagonal matrix with >1 entries on the diagonal. 
   * The solution found in [3] is applied.
   */

  /* solution of quartic polynomial of y = I_U = tr(U) 
   *     y^4 + p*y^2 + q*y + r = 0
   * with */
  p = -2.0*ci;  /* p = -2*I_C */
  q = -8.0*sqrt(ciii);  /* q = -8*(III_C)^{1/2} likely always < 0 */
  r = ci*ci - 4.0*cii;  /* r = I_C^2 - 2*II_C */

  /* associated cubic resolvent
   *     z^3 + 2*p*z^2 + (p^2-4*r)*z - q^2 = 0
   */

  /* normalised cubic resolvent with z = x - 2*p/3 (or x = z + 2*p/3)
   *     x^3 + pp*x + qq = 0
   * and */
  pp = -(p*p + 12.0*r)/3.0;
  qq = (-2.0*p*p*p + 72.0*p*r - 27.0*q*q)/27.0;
  /* solution with Cardan's formulae */
  disc = (4.0*pp*pp*pp + 27.0*qq*qq)/108.0;  /* discriminant */
  /* discriminant==0  ==>  3 real roots in x */
  if (abs(disc) < EPS12)
  {
    /* triple real root */
    if ( (abs(pp) < EPS12) && (abs(qq) < EPS12) )
    {
      x1 = 0.0;  /* triple real solution */
      z1 = x1 - 2.0*p/3.0;  /* triple real solution */
      z1rt = sqrt(z1);
      if (abs(-z1rt*z1rt*z1rt - q) < EPS12)
      {
        ui = 1.5*z1rt;
      }
      else
      {
        dserror("Trouble with radicals\n");
      }
    }
    /* 1 real root and 1 real double x root */
    else
    {
      /* should always be the case, but you never know ... */
      if (pp < 0.0)  
      {
        /* roots of normal form of cubic resolvent */
        x1 = -2.0*sqrt(-pp/3.0);  /* single root */
        x2 = sqrt(-pp/3.0);  /* double root */
        /* roots of cubic resolvent */
        z1 = x1 - 2.0*p/3.0;
        z2 = x2 - 2.0*p/3.0;
        /* radicals of roots of cubic resolvent */
        z1rt = sqrt(z1);  /* single */
        z2rt = sqrt(z2);  /* double */
        if (abs(-z1rt*z2rt*z2rt - q) < EPS12)
        {
          ui = 0.5*(z1rt + 2.0*z2rt);
        }
        else
        {
          dserror("Trouble with radicals\n");
        }
      }
      else
      {
        dserror("Error in finding roots\n");
      }
    }
  }
  /* discriminant<0  ==>  3 real roots in x  */
  else if (disc < 0.0) 
  {
    rho = sqrt(-pp*pp*pp/27.0);
    phi = acos(-qq/2.0/rho);
    rhort = 2.0 * pow(rho, 1.0/3.0);
    /* roots of normal form of cubic resolvent */
    x1 = rhort * cos(phi/3.0);
    x2 = rhort * cos(phi/3.0 + 2.0*PI/3.0);
    x3 = rhort * cos(phi/3.0 + 4.0*PI/3.0);
    /* roots of cubic resolvent */
    z1 = x1 - 2.0*p/3.0;
    z2 = x2 - 2.0*p/3.0;
    z3 = x3 - 2.0*p/3.0;
    /* radicals of roots of cubic resolvent */
    z1rt = sqrt(z1);
    z2rt = sqrt(z2);
    z3rt = sqrt(z3);
    if (abs(-z1rt*z2rt*z3rt - q) < EPS12)
    {
      ui = 0.5*(z1rt + z2rt + z3rt);
    }
    else
    {
      dserror("Trouble with radicals\n");
    }
  }
  /* discriminant>0  ==>  1 real and 2 complex roots in x */
  else
  {
    /* 1 real and 2 conjugated complex x roots */
    dserror("Discriminant is positive!\n");
  }
#endif


  /*--------------------------------------------------------------------*/
  /* 2nd and 3rd invariant of material stretch tensor U */
  /* 2nd invariant: II_U = 1/2 * (I_U^2 - I_C) */
  uii = 0.5*(ui*ui - ci);
  /* 3rd invariant: III_U = det(U) = sqrt(III_C) */
  uiii = sqrt(ciii);

  /*--------------------------------------------------------------------*/
  /* inverse of material stretch tensor U^{-1} */
  /* Hoger & Carlson [1] identified
   *     U^{-1} = [ III_U^2*(III_U+I_U*I_C) 
   *                + I_U^2*(I_U*III_C + III_U*II_C) ]^{-1}
   *            * [ I_U*(I_U*II_U - III_U)*C^2
   *                - (I_U*II_U - III_U)*(III_U +I_U*I_C)*C
   *                + { II_U*III_U*(III_U+I_U*I_C)
   *                    + I_U^2*(II_U*II_C+III_C) }*I ]
   *            = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
   */
  denom = uiii*uiii*(uiii + ui*ci)
    + ui*ui*(ui*ciii + uiii*cii);
  pc2 = ui*(ui*uii - uiii);
  pc = -(ui*uii - uiii)*(uiii + ui*ci);
  pi = uii*uiii*(uiii + ui*ci)
    + ui*ui*(uii*cii + ciii);
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      invut[i][j] = pc2*c2t[i][j] + pc*ct[i][j] + pi*it[i][j];
      invut[i][j] /= denom;
    }
  }

  /*--------------------------------------------------------------------*/
  /* rotation tensor R = F . U^{-1} */
  if (rt != NULL)
  {
    so3_tns3_dotprod(ft, invut, rt);
  }
  else if (vt != NULL)
  {
    so3_tns3_dotprod(ft, invut, rtx);
  }

  /*--------------------------------------------------------------------*/
  /* material stretch tensor U */
  /* Hoger & Carlson [1] wrote
   *     U = [ II_U * { II_U * (II_U + I_C) } + III_C ]^{-1}
   *       * [ -(I_U*II_U - III_U)*C^2
   *           + (I_U*II_U - III_U)*(II_U + I_C)*C
   *           + { I_U*III_U + III_U * ( II_U*(II_U+I_C)+II_C ) }*I ]
   *       = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
   */
  /* alternative:
   * U could be calculated based on R later: U = R^T . F
   */
  if (ut != NULL)
  {
    denom = uii*(uii*(uii+ci) + cii) + ciii;
    pc2 = -(ui*uii - uiii);
    pc = (ui*uii - uiii)*(uii + ci);
    pi = ui*ciii + uiii*(uii*(uii+ci) + cii);
    for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
      {
        ut[i][j] = pc2*c2t[i][j] + pc*ct[i][j] + pi*it[i][j];
        ut[i][j] /= denom;
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* spatial stretch tensor V could be  V = F . R^T */
  if (vt != NULL)
  {
    if (rt != NULL)
    {
      so3_tns3_dotprod_tr(ft, rt, vt);
    }
    else
    {
      so3_tns3_dotprod_tr(ft, rtx, vt);
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void so3_tns3_plrdcmp */



/*======================================================================*/
/*!
\brief Reshape symmetric tensor to equivalent vector

The tensor 
         [ A_11  A_12  A_31 ]
     A = [ A_12  A_22  A_23 ]
         [ A_31  A_23  A_33 ]
is stored as
     Av^T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]

\param  at     DOUBLE[][]    (i)  tensor A
\param  av     DOUBLE[]      (o)  vector Av
\return void

\author bborn
\date 01/07
*/
void so3_tns3_tsym2v(DOUBLE at[3][3],
                     DOUBLE av[6])
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns3_tsym2v");
#endif

  av[0] = at[0][0];
  av[1] = at[1][1];
  av[2] = at[2][2];
  av[3] = at[0][1];
  av[4] = at[1][2];
  av[5] = at[2][0];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Reshape equivalent vector to symmetric tensor

The vector
     Av^T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]
is stored as tensor
         [ A_11  A_12  A_31 ]
     A = [ A_12  A_22  A_23 ]
         [ A_31  A_23  A_33 ]

\param  av     DOUBLE[]      (i)  vector Av 
\param  at     DOUBLE[][]    (o)  tensor A
\return void

\author bborn
\date 01/07
*/
void so3_tns3_v2tsym(DOUBLE av[6],
                     DOUBLE at[3][3])
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns3_v2tsym");
#endif


  at[0][0] = av[0];  at[0][1] = av[3];  at[0][2] = av[5];
  at[1][0] = av[3];  at[1][1] = av[1];  at[1][2] = av[4];
  at[2][0] = av[5];  at[2][1] = av[4];  at[2][2] = av[2];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Spectral decomposition

The tensor 
         [ A_11  A_12  A_13 ]
     A = [ A_21  A_22  A_23 ]
         [ A_31  A_32  A_33 ]

\param  at     DOUBLE[][]    (i)  tensor A
\return void

\author bborn
\date 01/07
*/
void so3_tns3_spcdcmp(DOUBLE at[3][3],  /* input tensor */
                      INT *err,
                      DOUBLE ew[3])     /* eigen values */
/*                      DOUBLE ev[3][3])  \/* eigen vectors *\/ */
{
  DOUBLE ai, aii, aiii;  /* invariants of A */
  DOUBLE p, q, disc, rho, phi, rhort;
  DOUBLE y1, y2, y3;
  DOUBLE lam1, lam2, lam3;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_tns3_spcdcmp");
#endif

  /* initialise error  status */
  *int = 0;

  /* invariants of A */
  so3_tns3_inva(at, &ai, &aii, &aiii);

  /* characteristic polynomial (already in normal form)
   *    lam^3 - I_A*lam^2 + II_A*lam - III_A = 0
   * with reduced equation y=lam-II_A/3 (or lam=y+II_A/3)
   *    y^3 + p*y + q = 0   with   p = (3*II_A - I_A^2)/3
   *                               q = (2*I_A^3 + 9*I_A*II_A + 27*III_A)/27
   */
  p = (3.0*aii - ai*ai)/3.0;
  q = (2.0*ai*ai*ai + 9.0*ai*aii + 27.0*aiii)/27.0;

  /* discriminant */
  disc = (4.0*p*p*p + 27.0*q*q)/108.0;

  /* discriminant==0  ==>  3 real roots in y */
  if (abs(disc) < EPS12)
  {
    /* triple real root */
    if ( (abs(p) < EPS12) && (abs(q) < EPS12) )
    {
      y1 = 0.0;  /* triple real solution */
      lam1 = y1 + aii/3.0;  /* triple real solution */
      lam2 = lam1;
      lam3 = lam1;
    }
    /* 1 real root and 1 double real y root */
    else
    {
      /* should always be the case, but you never know ... */
      if (p < 0.0)  
      {
        /* roots of normal form of cubic resolvent */
        y1 = -2.0*sqrt(-p/3.0);  /* single root */
        y2 = sqrt(-p/3.0);  /* double root */
        /* roots of cubic resolvent */
        lam1 = (3.0*y1 + 2.0*aii)/3.0;
        lam2 = (3.0*y2 + 2.0*aii)/3.0;
        lam3 = lam2;
      }
      else
      {
        dserror("Error in finding roots\n");
        *err = 1;
      }
    }
  }
  /* discriminant<0  ==>  3 different real roots in y  */
  else if (disc < 0.0)
  {
    rho = sqrt(-p*p*p/27.0);
    phi = acos(-q/2.0/rho);
    rhort = 2.0 * pow(rho, 1.0/3.0);
    /* roots of normal form of cubic resolvent */
    y1 = rhort * cos(phi/3.0);
    y2 = rhort * cos((phi + 2.0*PI)/3.0);
    y3 = rhort * cos((phi + 4.0*PI)/3.0);
    /* roots of cubic resolvent */
    lam1 = (3.0*y1 + 2.0*aii)/3.0;
    lam2 = (3.0*y2 + 2.0*aii)/3.0;
    lam3 = (3.0*y3 + 2.0*aii)/3.0;
  }
  else
  {
    /* 1 real and 2 conjugated complex y roots */
    dserror("Discriminant is positive!\n");
    *err = 1;
  }

  if (*err == 0)
  {
    ev[0] = lam1;
    ev[1] = lam2;
    ev[2] = lam3;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
