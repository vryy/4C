/*======================================================================*/
/*!
\file
\brief Vector manipulation: Here the vector is 6-dimensional and
       the vector form of a symmetric 2-tensor in 3D.

            [ a_11 ]
            [ a_22 ]
            [ a_33 ]            [ a_11  a_12  a_13 ]
       av = [ ~~~~ ]   of  at = [ a_12  a_22  a_23 ]
            [ a_13 ]            [ a_13  a_23  a_33 ]
            [ a_23 ]
            [ a_13 ]

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 04/07
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
\brief Zero vector
\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_zero(DOUBLE iv[NUMSTR_SOLID3])  /*!< input vector */
{
  INT i;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_zero");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    iv[i] = 0.0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Identity vector

Explicitly
            [  1  ]
            [  1  ]
            [  1  ]
       iv = [ ~~~ ]
            [  0  ]
            [  0  ]
            [  0  ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_id(DOUBLE iv[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_id");
#endif

  for (i=0; i<3; i++)
  {
    iv[i] = 1.0;
  }
  for (i=3; i<NUMSTR_SOLID3; i++)
  {
    iv[i] = 0.0;
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign vector by another vector
\author bborn
\date 04/07
*/
void so3_mv6_v_ass(const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                  DOUBLE bv[NUMSTR_SOLID3])  /*!< assigned vector */
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_ass");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    bv[i] = av[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign vector by another vector and scale it
\author bborn
\date 04/07
*/
void so3_mv6_v_assscl(const DOUBLE scl,  /*!< scale */
                     const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                     DOUBLE bv[NUMSTR_SOLID3])  /*!< assigned vector */
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_assscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    bv[i] = scl * av[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update vector by another vector and scale it
\author bborn
\date 04/07
*/
void so3_mv6_v_updscl(const DOUBLE scl,  /*!< scale */
                      const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                      DOUBLE bv[NUMSTR_SOLID3])  /*!< updated vector */
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_updscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    bv[i] += scl * av[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Trace of a vector

Explicitly
       tr(av) = a_11 + a_22 + a_33

\param    av      DOUBLE[]    (i)    input vector
\param    tr      DOUBLE*     (o)    trace
\return   void
\author bborn
\date 04/07
*/
void so3_mv6_v_tr(const DOUBLE av[NUMSTR_SOLID3],
                  DOUBLE *tr)
{
  INT i;  /* index */
  DOUBLE trace;  /* auxiliar trace */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_tr");
#endif

  trace = 0.0;
  for (i=0; i<3; i++)
  {
    trace += av[i];
  }
  *tr = trace;

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Determinant of vector
\param  av     DOUBLE[]    (i)    input vector
\param  det    DOUBLE*     (o)    determinant
\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_det(const DOUBLE av[NUMSTR_SOLID3],
                   DOUBLE *det)
{

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_det");
#endif

  /* determionant of 3x3-matrix by Sarrus' rule */
  *det = av[0] * av[1] * av[2]
       + av[3] * av[4] * av[5]
       + av[5] * av[3] * av[4]
       - av[5] * av[1] * av[5]
       - av[0] * av[4] * av[4]
       - av[3] * av[3] * av[2];

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Deviator of vector 

Explicitly
       dev(av) = av - 1/3 tr(av) iv
or
                 [ a_11 ]       [ a_11 + a_22 + a_33 ]
                 [ a_22 ]       [ a_11 + a_22 + a_33 ]
                 [ a_33 ]     1 [ a_11 + a_22 + a_33 ]
       dev(av) = [ ~~~~ ]  -  - [ ~~~~~~~~~~~~~~~~~~ ]
                 [ a_12 ]     3 [         0          ]
                 [ a_23 ]       [         0          ]
                 [ a_13 ]       [         0          ]
\param  av     DOUBLE[]     (i)     input vector
\param  adev   DOUBLE[]     (o)     deviator vector
\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_dev(const DOUBLE av[NUMSTR_SOLID3],
                   DOUBLE adev[NUMSTR_SOLID3])
{
  DOUBLE trace;  /* trace */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_dev");
#endif

  so3_mv6_v_tr(av, &trace);

  adev[0] = av[0] - trace/3.0;
  adev[1] = av[1] - trace/3.0;
  adev[2] = av[2] - trace/3.0;
  adev[3] = av[3];
  adev[4] = av[4];
  adev[5] = av[5];

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Subtract 2 vectors 

Explicitly
            [ a_11 ]     [ b_11 ]
            [ a_22 ]     [ b_22 ]
            [ a_33 ]     [ b_33 ]
       cv = [ ~~~~ ]  -  [ ~~~~ ]
            [ a_12 ]     [ b_12 ]
            [ a_23 ]     [ b_23 ]
            [ a_13 ]     [ b_13 ]

\param  av    DOUBLE[]     (i)    1st input vector
\param  bv    DOUBLE[]     (i)    2nd input vector
\param  cv    DOUBLE[]     (o)    output vector
\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_sub(const DOUBLE av[NUMSTR_SOLID3],
                   const DOUBLE bv[NUMSTR_SOLID3],
                   DOUBLE cv[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_sub");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    cv[i] = av[i] - bv[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Double contraction of 2 vectors 

Explicitly
       av : bv = a_11*b_11 + a_22*b_22 + a_33*b_33
               + 2*a_12*b_12 + 2*a_23*b_23 + 2*a_31*b_31

\param  av    DOUBLE[]     (i)    1st input vector
\param  bv    DOUBLE[]     (i)    2nd input vector
\param  prd   DOUBLE*      (o)    double product
\return void
\author bborn
\date 04/07
*/
void so3_mv6_v_dblctr(const DOUBLE av[NUMSTR_SOLID3],
                      const DOUBLE bv[NUMSTR_SOLID3],
                      DOUBLE* prd)
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v_dblprd");
#endif

  *prd = 0.0;
  for (i=0; i<3; i++)
  {
    *prd += av[i] * bv[i];
  }
  for (i=3; i<NUMSTR_SOLID3; i++)
  {
    *prd += 2.0 * av[i] * bv[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Matrix vector product

Explicitly
    Cv = Am . Bv

\param  am       DOUBLE[][]   (i)    input matrix
\param  bv       DOUBLE[]     (i)    input vector
\param  cv       DOUBLE[]     (o)    output vector
\author bborn
\date 04/07
*/
void so3_mv6_v_assmvp(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      const DOUBLE bv[NUMSTR_SOLID3],
                      DOUBLE cv[NUMSTR_SOLID3])
{
  INT i, j;  /* indices */
  DOUBLE rcsum;  /* intermediate row * column sum */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_sub");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    rcsum = 0.0;  /* intermediate row * column sum */
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      rcsum += am[i][j] * bv[j];
    }
    cv[i] = rcsum;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign vector by another vector and scale it
       Due to the fact that strain vector component have a doubled
       shear components (cf. solid3.h), but stress vectors have not,
       we need a special assignment function to scale shear components
       accordingly.
\param   scl     DOUBLE    (i)    scale
\param   av      DOUBLE[]  (i)    input vector ('stress-vector-like')
\param   bv      DOUBLE[]  (o)    assigned vector ('strain-vector-like')
\return  void
\author bborn
\date 04/07
*/
void so3_mv6_v2_assscl(const DOUBLE scl,
                       const DOUBLE av[NUMSTR_SOLID3],
                       DOUBLE bv[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v2_assscl");
#endif

  for (i=0; i<3; i++)
  {
    bv[i] = scl * av[i];
  }
  for (i=3; i<NUMSTR_SOLID3; i++)
  {
    bv[i] = scl * 2.0 * av[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign vector by another vector and scale it
       Due to the fact that strain vector component have a doubled
       shear components (cf. solid3.h), but stress vectors have not,
       we need a special assignment function to scale shear components
       accordingly.
\param   scl     DOUBLE    (i)    scale
\param   av      DOUBLE[]  (i)    input vector ('strain-vector-like')
\param   bv      DOUBLE[]  (o)    assigned vector ('stress-vector-like')
\return  void
\author bborn
\date 04/07
*/
void so3_mv6_v05_assscl(const DOUBLE scl,
                        const DOUBLE av[NUMSTR_SOLID3],
                        DOUBLE bv[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v05_assscl");
#endif

  for (i=0; i<3; i++)
  {
    bv[i] = scl * av[i];
  }
  for (i=3; i<NUMSTR_SOLID3; i++)
  {
    bv[i] = scl * 0.5 * av[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign vector by another vector and scale it
       Due to the fact that strain vector component have a doubled
       shear components (cf. solid3.h), i.e.
          av = [ a11 a22 a33 | 2*a12 2*a23 2*a31 ]
       but stress vectors have not, i.e.
          av = [ a11 a22 a33 | a12 a23 a31 ]
       we need to scale the last three entries.
\param   av      DOUBLE[]  (io)   input vector ('strain-vector-like')
                                  output vector ('stress-vector-like')
\return  void
\author bborn
\date 05/07
*/
void so3_mv6_v05_updvtov05(DOUBLE av[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_v05_updvtov05");
#endif

  for (i=NDIM_SOLID3; i<NUMSTR_SOLID3; i++)
  {
    av[i] *= 0.5;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Zero matrix
\param   am      DOUBLE[][]   (i)   input matrix
\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_zero(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3])
{

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_zero");
#endif

  memset(am, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Identity matrix

Explicitly
            [  1  0  0  |  0  0  0  ]
            [  0  1  0  |  0  0  0  ]
            [  0  0  1  |  0  0  0  ]
       im = [ ~~~~~~~~~   ~~~~~~~~~ ]
            [  0  0  0  |  1  0  0  ]
            [  0  0  0  |  0  1  0  ]
            [  0  0  0  |  0  0  1  ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_id(DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_id");
#endif

  memset(im, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    im[i][i] = 1.0;
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Identity matrix

Explicitly
            [  s  0  0  |  0  0  0  ]
            [  0  s  0  |  0  0  0  ]
            [  0  0  s  |  0  0  0  ]
       iv = [ ~~~~~~~~~   ~~~~~~~~~ ]
            [  0  0  0  |  s  0  0  ]
            [  0  0  0  |  0  s  0  ]
            [  0  0  0  |  0  0  s  ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_idscl(const DOUBLE scale,
                     DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_idscl");
#endif

  memset(im, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    im[i][i] = scale;
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



/*======================================================================*/
/*!
\brief One matrix

Explicitly
            [  1  1  1  |  0  0  0  ]
            [  1  1  1  |  0  0  0  ]
            [  1  1  1  |  0  0  0  ]
       om = [ ~~~~~~~~~   ~~~~~~~~~ ]
            [  0  0  0  |  0  0  0  ]
            [  0  0  0  |  0  0  0  ]
            [  0  0  0  |  0  0  0  ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_one(DOUBLE om[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_one");
#endif

  memset(om, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      om[i][j] = 1.0;
    }
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update matrix with scaled one-matrix

Explicitly
                   [  1  1  1  |  0  0  0  ]
                   [  1  1  1  |  0  0  0  ]
                   [  1  1  1  |  0  0  0  ]
       am += scl * [ ~~~~~~~~~   ~~~~~~~~~ ]
                   [  0  0  0  |  0  0  0  ]
                   [  0  0  0  |  0  0  0  ]
                   [  0  0  0  |  0  0  0  ]

\param  scl        DOUBLE     (i)    scale
\param  am         DOUBLE[][] (io)   updated matrix
\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_updonescl(const DOUBLE scl,
                         DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_updonescl");
#endif

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      am[i][j] += scl;
    }
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign matrix by another matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_ass(const DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indeces */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_ass");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      bm[i][j] = am[i][j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



/*======================================================================*/
/*!
\brief Assign vector by another vector and scale it
\author bborn
\date 04/07
*/
void so3_mv6_m_assscl(const DOUBLE scl,  /*!< scale */
                      DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_assscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      bm[i][j] = scl * am[i][j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update vector by another vector and scale it
\param  scl       DOUBLE    (i)    scale
\param  am        DOUBLE[]  (i)    input matrix
\param  bm        DOUBLE[]  (o)    output matrix
\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_updscl(const DOUBLE scl,
                      DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                      DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_updscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      bm[i][j] += scl * am[i][j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Subtract 2 matrices 
\return void
\author bborn
\date 04/07
*/
void so3_mv6_m_sub(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                   DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_sub");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      cm[i][j] = am[i][j] - bm[i][j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Dyadic vector product

Explicitly
    Cm = Av . Bv^T

\param  av       DOUBLE[]     (i)    input vector
\param  bv       DOUBLE[]     (i)    input vector
\param  cm       DOUBLE[][]   (o)    output matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_assdyd(const DOUBLE av[NUMSTR_SOLID3],
                 const DOUBLE bv[NUMSTR_SOLID3],
                 DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_dyd");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      cm[i][j] = av[i] * bv[j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update matrix with scaled dyadic vector product

Explicitly
    Cm = Cm  +  scl * Av . Bv^T

\param  scl      DOUBLE       (i)    scale
\param  av       DOUBLE[]     (i)    input vector
\param  bv       DOUBLE[]     (i)    input vector
\param  cm       DOUBLE[][]   (io)   updated matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_upddydscl(const DOUBLE scl,
                         const DOUBLE av[NUMSTR_SOLID3],
                         const DOUBLE bv[NUMSTR_SOLID3],
                         DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_upddydscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      cm[i][j] += scl * av[i] * bv[j];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign matrix by matrix-matrix product (inner product)

Explicitly
    Cm = Am . Bm

\param  am       DOUBLE[][]   (i)    output matrix
\param  bm       DOUBLE[][]   (i)    output matrix
\param  cm       DOUBLE[][]   (o)    assigned matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_mprd(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j, k;  /* indices */
  DOUBLE rcsum;

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_mprd");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      rcsum = 0.0;
      for (k=0; k<NUMSTR_SOLID3; k++)
      {
        rcsum += am[i][k] * bm[k][j];
      }
      cm[i][j] = rcsum;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign matrix by scaled matrix-matrix product (inner product)

Explicitly
    Cm = scl * Am . Bm

\param  scl      DOUBLE       (i)    scale
\param  am       DOUBLE[][]   (i)    output matrix
\param  bm       DOUBLE[][]   (i)    output matrix
\param  cm       DOUBLE[][]   (o)    assigned matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_mprdscl(const DOUBLE scl,
                       DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j, k;  /* indices */
  DOUBLE rcsum;

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_mprdscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      rcsum = 0.0;
      for (k=0; k<NUMSTR_SOLID3; k++)
      {
        rcsum += am[i][k] * bm[k][j];
      }
      cm[i][j] = scl * rcsum;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Assign matrix by scaled matrix-matrix product (inner product)

Explicitly
    Cm += Am . Bm

\param  am       DOUBLE[][]   (i)    output matrix
\param  bm       DOUBLE[][]   (i)    output matrix
\param  cm       DOUBLE[][]   (o)    assigned matrix
\author bborn
\date 04/07
*/
void so3_mv6_m_updmprd(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE bm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                       DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i, j, k;  /* indices */
  DOUBLE rcsum;

#ifdef DEBUG
  dstrc_enter("so3_mv6_m_mprdscl");
#endif

  for (i=0; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      rcsum = 0.0;
      for (k=0; k<NUMSTR_SOLID3; k++)
      {
        rcsum += am[i][k] * bm[k][j];
      }
      cm[i][j] += rcsum;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Identity matrix

Explicitly
            [  s  0  0  |   0   0   0  ]
            [  0  s  0  |   0   0   0  ]
            [  0  0  s  |   0   0   0  ]
       iv = [ ~~~~~~~~~   ~~~~~~~~~~~~ ]
            [  0  0  0  |  2*s  0   0  ]
            [  0  0  0  |   0  2*s  0  ]
            [  0  0  0  |   0   0  2*s ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_m2_idscl(const DOUBLE scale,
                      DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m2_idscl");
#endif

  memset(im, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

  for (i=0; i<NDIM_SOLID3; i++)
  {
    im[i][i] = scale;
  }
  for (i=NDIM_SOLID3; i<NUMSTR_SOLID3; i++)
  {
    im[i][i] = scale * 2.0;
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update m-matrix to m2-matrix
\author bborn
\date 05/07
*/
void so3_mv6_m2_updmtom2(DOUBLE am[NUMSTR_SOLID3][NUMSTR_SOLID3]) 
{
  INT i, j;

#ifdef DEBUG
  dstrc_enter("so3_mv6_m2_updmtom");
#endif

  for (i=NDIM_SOLID3; i<NUMSTR_SOLID3; i++)
  {
    for (j=0; j<NUMSTR_SOLID3; j++)
    {
      am[i][j] *= 2.0;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Identity m05-matrix

Explicitly
            [  s  0  0  |   0      0     0   ]
            [  0  s  0  |   0      0     0   ]
            [  0  0  s  |   0      0     0   ]
       iv = [ ~~~~~~~~~   ~~~~~~~~~~~~~~~~~~ ]
            [  0  0  0  |  0.5*s   0     0   ]
            [  0  0  0  |   0    0.5*s   0   ]
            [  0  0  0  |   0      0   0.5*s ]

\return void
\author bborn
\date 04/07
*/
void so3_mv6_m05_idscl(const DOUBLE scale,
                       DOUBLE im[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_mv6_m05_idscl");
#endif

  memset(im, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));

  for (i=0; i<NDIM_SOLID3; i++)
  {
    im[i][i] = scale;
  }
  for (i=NDIM_SOLID3; i<NUMSTR_SOLID3; i++)
  {
    im[i][i] = scale * 0.5;
  }
  

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
