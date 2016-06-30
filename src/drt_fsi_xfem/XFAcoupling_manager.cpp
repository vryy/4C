/*!----------------------------------------------------------------------
\file XFAcoupling_manager.cpp
\brief Coupling Manager for eXtended Fluid Ale Coupling

\level 3

<pre>
\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>

*----------------------------------------------------------------------*/
#include "XFAcoupling_manager.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_adapter/ad_ale_fpsi.H"
#include "../linalg/linalg_mapextractor.H"

//#include "../drt_lib/drt_globalproblem.H"
//#include "../drt_io/io.H"

/*-----------------------------------------------------------------------------------------*
| Constructor                                                                 ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
XFEM::XFACoupling_Manager::XFACoupling_Manager(Teuchos::RCP<FLD::XFluid> xfluid,Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale, std::vector<int> idx)
:Coupling_Comm_Manager(*xfluid->Discretization(),*ale->Discretization(),"",0,3),
 ale_(ale),
 xfluid_(xfluid),
 idx_(idx)
{
  if (idx_.size() != 2) dserror("XFACoupling_Manager required two block ( 2 != %d)", idx_.size());
}

/*-----------------------------------------------------------------------------------------*
| Set required displacement & velocity states in the coupling object          ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFACoupling_Manager::SetCouplingStates()
{

  //1 Sets structural conditioned Dispnp onto Ale
  // --> this will be done in one of the following commits! Ale_Struct_coupling_->InsertVector(1,StructureField()->Dispnp(),0,AleField()->WriteAccessDispnp(),XFEM::Coupling_Comm_Manager::full_to_full);
  //2 Get AleDisplacements
  Teuchos::RCP<Epetra_Vector> aledisplacements =  Teuchos::rcp(new Epetra_Vector(*GetMapExtractor(0)->Map(1),true));
  InsertVector(1,ale_->Dispnp(),0,aledisplacements,Coupling_Comm_Manager::partial_to_partial);
  //3 Set Fluid Dispnp
  GetMapExtractor(0)->InsertVector(aledisplacements,1,xfluid_->WriteAccessDispnp());

  //4 new grid velocity
  xfluid_->UpdateGridv();
  return;
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling matrixes to the global systemmatrix                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFACoupling_Manager::AddCouplingMatrix(LINALG::BlockSparseMatrixBase& systemmatrix, double scaling)
{
  //Get Idx of fluid and ale field map extractors
  const int &aidx_other = ALE::UTILS::MapExtractor::cond_other;

  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a      = ale_->BlockSystemMatrix();

  //ALE Condensation
  LINALG::SparseMatrix& aii = a->Matrix(aidx_other,aidx_other);

  systemmatrix.Assign(idx_[1],idx_[1],LINALG::View,aii);

  //  //////////////////////////////////////////////
  //  //////                                  //////
  //  //////    Linearization of FluidField   //////
  //  //////    with respect to ale mesh      //////
  //  //////             motion               //////
  //  //////                                  //////
  //  //////////////////////////////////////////////

  // TODO: THIS IS STILL MISSING, BUT USUALLY DOES NOT HAVE A BIG INFLUENCE ON THE CONVERGENCE!!!
  return;
}

/*-----------------------------------------------------------------------------------------*
| Add the coupling rhs                                                        ager 06/2016 |
*-----------------------------------------------------------------------------------------*/
void XFEM::XFACoupling_Manager::AddCouplingRHS(Teuchos::RCP<Epetra_Vector> rhs,
    const LINALG::MultiMapExtractor& me, double scaling)
{
  Teuchos::RCP<const Epetra_Vector> av = ale_->RHS();
  Teuchos::RCP<Epetra_Vector> aov = ale_->Interface()->ExtractOtherVector(av);
  me.InsertVector(*aov,idx_[1],*rhs);  // add ALE contributions to 'rhs'
  return;
}
