/*----------------------------------------------------------------------*/
/*! \file
\brief Utilities for the inverse analysis

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_UTILS_H
#define INVANA_UTILS_H

#if 0
#define INVANA_DEBUG_PRINTING 1
#else
#define INVANA_DEBUG_PRINTING 0
#endif

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "chol_factor_base.H"
#include "../drt_timestepping/timintmstep.H"


class Epetra_CrsMatrix;

namespace INVANA
{
  //! Helper functions useable throughout the inverse analysis

  //! compute some norm on a vector stored in multivector format for storage reasons only
  void MVNorm(
      const Epetra_MultiVector& avector, const Epetra_Map& uniquemap, int anorm, double* result);

  //! compute dot product of two vectors stored in multivector format for storage reasons only
  void MVDotProduct(const Epetra_MultiVector& avector, const Epetra_MultiVector& bvector,
      const Epetra_Map& uniquemap, double* result);

  Teuchos::RCP<INVANA::CholFactorBase> CreateICT(
      Teuchos::RCP<Epetra_CrsMatrix> covariance,  /// [in] covariance to tbe factorized
      const Teuchos::ParameterList& params);

}  // namespace INVANA

#endif /*INVANA_UTILS_H_*/
