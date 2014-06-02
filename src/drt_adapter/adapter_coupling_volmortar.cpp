/*!----------------------------------------------------------------------
\file adapter_coupling_volmortar.cpp

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  farah 10/13|
 *----------------------------------------------------------------------*/
#include "adapter_coupling_volmortar.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_volmortar/volmortar_coupling.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_particle/binning_strategy.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_io/io_gmsh.H"
/*----------------------------------------------------------------------*
 |  ctor                                                     farah 10/13|
 *----------------------------------------------------------------------*/
ADAPTER::MortarVolCoupl::MortarVolCoupl()
{
  //empty...
}

/*----------------------------------------------------------------------*
 |  setup                                                    farah 10/13|
 *----------------------------------------------------------------------*/
void ADAPTER::MortarVolCoupl::Setup(Teuchos::RCP<DRT::Discretization> slavedis,
                                    Teuchos::RCP<DRT::Discretization> masterdis)
{
  // get problem dimension (2D or 3D) and create (MORTAR::MortarInterface)
  const int dim = DRT::Problem::Instance()->NDim();

  // redistribute discr. with help of binning strategy
  if(slavedis->Comm().NumProc()>1)
  {
    // create vector of discr.
    std::vector<Teuchos::RCP<DRT::Discretization> > dis;
    dis.push_back(slavedis);
    dis.push_back(masterdis);

    /// binning strategy is created and parallel redistribution is performed
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
        Teuchos::rcp(new BINSTRATEGY::BinningStrategy(dis));
  }

  // create coupling instance
  Teuchos::RCP<VOLMORTAR::VolMortarCoupl> coupdis =
      Teuchos::rcp(new VOLMORTAR::VolMortarCoupl(dim,slavedis,masterdis));

  // Evaluate coupling:
  // 1. perform cut/ integration cell identification
  //    OR call ele-based integration
  // 2. integrate cells
  // 3. assemble mortar matrices
  // 4. create projector P
  coupdis->Evaluate();

  // get the P operators
  pmatrixA_ = coupdis->GetPMatrixAB();
  pmatrixB_ = coupdis->GetPMatrixBA();

  return;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping                                             farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::MortarVolCoupl::ApplyVectorMappingAB(Teuchos::RCP<const Epetra_Vector> vec) const
{
  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(pmatrixA_->RowMap(),true);
  pmatrixA_->Multiply(false,*vec,*mapvec);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping                                             farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::MortarVolCoupl::ApplyVectorMappingBA(Teuchos::RCP<const Epetra_Vector> vec) const
{
  Teuchos::RCP<Epetra_Vector> mapvec = LINALG::CreateVector(pmatrixB_->RowMap(),true);
  pmatrixB_->Multiply(false,*vec,*mapvec);

  return mapvec;
}

/*----------------------------------------------------------------------*
 |  ApplyMapping                                             farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMappingAB(Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  return LINALG::MLMultiply(*mat,false,*pmatrixA_,false,false,false,true);
}

/*----------------------------------------------------------------------*
 |  ApplyMapping                                             farah 01/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::MortarVolCoupl::ApplyMatrixMappingBA(Teuchos::RCP<const LINALG::SparseMatrix> mat) const
{
  return LINALG::MLMultiply(*mat,false,*pmatrixB_,false,false,false,true);
}
