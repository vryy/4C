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


#include "invana_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "Epetra_SerialDenseVector.h"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
/* Compute norm of two vectors stored in                    keh 03/14   */
/* multivector format for storage reasons only                          */
/*----------------------------------------------------------------------*/
void STR::INVANA::MVNorm(Teuchos::RCP<Epetra_MultiVector> avector, int anorm, double* result, const Epetra_Map* uniquemap)
{
  Epetra_SerialDenseVector vnorm(avector->NumVectors());

  Teuchos::RCP<Epetra_MultiVector> unique = Teuchos::rcp(new Epetra_MultiVector(*uniquemap, avector->NumVectors(),true));
  LINALG::Export(*avector, *unique);

  if (anorm==2)
    unique->Norm2(vnorm.Values());
  else if (anorm==0)
    unique->NormInf(vnorm.Values());
  else if(anorm==1)
    unique->Norm1(vnorm.Values());
  else
    dserror("provide norm to be computed: 0 (inf-Norm), 1 (1-Norm) or 2 (2-Norm)");

  *result = vnorm.Norm2();
}


/*----------------------------------------------------------------------*/
/* Compute dot product of two vectors stored in             keh 03/14   */
/* multivector format for storage reasons only                          */
/*----------------------------------------------------------------------*/
//! compute dot product of two vectors stored in multivector fomat for storage reasons only
void STR::INVANA::MVDotProduct(Teuchos::RCP<Epetra_MultiVector> avector, Teuchos::RCP<Epetra_MultiVector> bvector, double* result, const Epetra_Map* uniquemap)
{
  dsassert(avector->NumVectors()==bvector->NumVectors(), "give proper multivectors!");

  Epetra_SerialDenseVector anorm(avector->NumVectors());

  Teuchos::RCP<Epetra_MultiVector> uniquea = Teuchos::rcp(new Epetra_MultiVector(*uniquemap, avector->NumVectors(),true));
  LINALG::Export(*avector, *uniquea);
  Teuchos::RCP<Epetra_MultiVector> uniqueb = Teuchos::rcp(new Epetra_MultiVector(*uniquemap, bvector->NumVectors(),true));
  LINALG::Export(*bvector, *uniqueb);

  // do dot product with unique vectors now
  uniquea->Dot(*uniqueb,anorm.Values());

  *result=0.0;
  for (int j=0; j<anorm.Length(); j++) *result+=anorm[j];

}

