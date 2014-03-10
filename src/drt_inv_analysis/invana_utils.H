/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_utils.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_UTILS_H
#define INVANA_UTILS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"


namespace STR
{

namespace INVANA
{

//! Helper functions useable throughout the inverse analysis
//! compute some norm on a vector stored in multivector format for storage reasons only
void MVNorm(Teuchos::RCP<Epetra_MultiVector> avector, int anorm, double* result, const Epetra_Map* uniquemap);

//! compute dot product of two vectors stored in multivector fomat for storage reasons only
void MVDotProduct(Teuchos::RCP<Epetra_MultiVector> avector, Teuchos::RCP<Epetra_MultiVector> bvector, double* result, const Epetra_Map* uniquemap);


}  // namespace INVANA
}  // namespace STR

#endif /*STAT_INV_UTILS_H_*/
