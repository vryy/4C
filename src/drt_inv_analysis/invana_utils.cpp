/*----------------------------------------------------------------------*/
/*!
\file invana_utils.cpp
\brief Utilities for the inverse analysis
\level 3
\maintainer Sebastian Kehl
             kehl@mhpc.mw.tum.de
             089 - 289-10361
*/
/*----------------------------------------------------------------------*/

#include "invana_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "Epetra_SerialDenseVector.h"
#include "../linalg/linalg_utils.H"
#include "DcsMatrix.H"
#include "chol_factory.H"
#include "../drt_inpar/inpar_statinvanalysis.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include <fstream>


/*----------------------------------------------------------------------*/
/* Compute norm of two vectors stored in                    keh 03/14   */
/* multivector format for storage reasons only                          */
/*----------------------------------------------------------------------*/
void INVANA::MVNorm(const Epetra_MultiVector& avector, const Epetra_Map& uniquemap, int anorm, double* result)
{
  Epetra_SerialDenseVector vnorm(avector.NumVectors());

  Teuchos::RCP<Epetra_MultiVector> unique = Teuchos::rcp(new Epetra_MultiVector(uniquemap, avector.NumVectors(),true));
  LINALG::Export(avector, *unique);

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
void INVANA::MVDotProduct(const Epetra_MultiVector& avector, const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result)
{
  dsassert(avector.NumVectors()==bvector.NumVectors(), "give proper multivectors!");

  Epetra_SerialDenseVector anorm(avector.NumVectors());

  Teuchos::RCP<Epetra_MultiVector> uniquea = Teuchos::rcp(new Epetra_MultiVector(uniquemap, avector.NumVectors(),true));
  LINALG::Export(avector, *uniquea);
  Teuchos::RCP<Epetra_MultiVector> uniqueb = Teuchos::rcp(new Epetra_MultiVector(uniquemap, bvector.NumVectors(),true));
  LINALG::Export(bvector, *uniqueb);

  // do dot product with unique vectors now
  uniquea->Dot(*uniqueb,anorm.Values());

  *result=0.0;
  for (int j=0; j<anorm.Length(); j++) *result+=anorm[j];

}


/*----------------------------------------------------------------------*/
Teuchos::RCP<INVANA::CholFactorBase> INVANA::CreateICT_lowmem(
    Teuchos::RCP<Epetra_CrsMatrix> covariance, const Teuchos::ParameterList& params)
{
  INVANA::CholFactory cholfactory;

  Teuchos::ParameterList List(params);
  // These four are read by Ifpack_IC::SetParameters
  List.set("fact: ict level-of-fill", params.get<double>("MAP_COV_FILL")); // exact: 1.0
  List.set("fact: absolute threshold",0.0); // exact: 0.0
  List.set("fact: relative threshold",1.0); //exact: 1.0
  List.set("fact: drop tolerance",0.0); //exact: 0.0
  Teuchos::RCP<CholFactorBase> prec = cholfactory.Create(covariance,List);

  // -------- debug print outs
#if INVANA_DEBUG_PRINTING

  // test for type
  Teuchos::RCP<DcsMatrix> Ad = Teuchos::rcp_dynamic_cast<DcsMatrix>(covariance);

  if (Ad != Teuchos::null)
  {
  LINALG::PrintMatrixInMatlabFormat("lower.mtl",prec->H());

  Teuchos::RCP<Epetra_CrsMatrix> fullA = Ad->FillMatrix();
  LINALG::PrintMatrixInMatlabFormat("fullmatrix.mtl",*fullA);
  }
  else
    LINALG::PrintMatrixInMatlabFormat("fullmatrix.mtl",*covariance);

  Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(covariance->RowMap(),false));
  b->Random();
  LINALG::PrintVectorInMatlabFormat("rhs.mtl",*b);

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(covariance->RowMap(),true));
  prec->ApplyInverse(*b,*x);
  LINALG::PrintVectorInMatlabFormat("x.mtl",*x);
#endif
  // -------- end debug

  return prec;
}
