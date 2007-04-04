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
void so3_vct6_zero(DOUBLE iv[NUMSTR_SOLID3])  /*!< input vector */
{
  INT i;  /* indices */

#ifdef DEBUG
  dstrc_enter("so3_vct6_zero");
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
void so3_vct6_id(DOUBLE iv[NUMSTR_SOLID3])
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_vct6_id");
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
\brief Trace of a vector

Explicitly
       tr(av) = a_11 + a_22 + a_33

\return void
\author bborn
\date 04/07
*/
void so3_vct6_tr(const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                 DOUBLE *tr)  /*!< trace */
{
  INT i;  /* index */
  DOUBLE trace;  /* auxiliar trace */

#ifdef DEBUG
  dstrc_enter("so3_vct6_tr");
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
\return void
\author bborn
\date 04/07
*/
void so3_vct6_det(const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                  DOUBLE *det)  /*!< determinant */
{

#ifdef DEBUG
  dstrc_enter("so3_vct6_det");
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
                 [ a_12 ]     3 [        a_12        ]
                 [ a_23 ]       [        a_23        ]
                 [ a_13 ]       [        a_13        ]

\return void
\author bborn
\date 04/07
*/
void so3_vct6_dev(const DOUBLE av[NUMSTR_SOLID3],  /*!< input vector */
                  DOUBLE adev[NUMSTR_SOLID3])  /*!< deviator vector */
{
  DOUBLE trace;  /* trace */

#ifdef DEBUG
  dstrc_enter("so3_vct6_dev");
#endif

  so3_vct6_tr(av, &trace);

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

\return void
\author bborn
\date 04/07
*/
void so3_vct6_sub(const DOUBLE av[NUMSTR_SOLID3],  /*!< 1st input vector */
                  const DOUBLE bv[NUMSTR_SOLID3],  /*!< 2nd input vector */
                  DOUBLE cv[NUMSTR_SOLID3])  /*!< output vector */
{
  INT i;  /* index */

#ifdef DEBUG
  dstrc_enter("so3_vct6_sub");
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
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
