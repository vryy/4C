/*----------------------------------------------------------------------*/
/*!
 * \file invana_utils.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_UTILS_H
#define INVANA_UTILS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"


namespace STR
{

namespace INVANA
{

//! Helper functions useable throughout the inverse analysis

//! compute some norm on a vector stored in multivector format for storage reasons only
void MVNorm(const Epetra_MultiVector& avector, const Epetra_Map& uniquemap, int anorm, double* result);

//! compute dot product of two vectors stored in multivector fomat for storage reasons only
void MVDotProduct(const Epetra_MultiVector& avector, const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result);

}  // namespace INVANA
}  // namespace STR

#endif /*INVANA_UTILS_H_*/
