/*!----------------------------------------------------------------------
\file mortar_multifield_coupling.cpp

\brief Class performing coupling (condensation/recovery) for dual mortar
       methods in (volume) monolithic multi-physics applications, i.e. in
       block matrix systems. This also accounts for the correct condensation
       in the off-diagonal matrix blocks

\level 2

\maintainer Alexander Seitz

*-----------------------------------------------------------------------*/
#include "mortar_multifield_coupling.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "mortar_utils.H"

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::PushBackCoupling(
    const Teuchos::RCP<DRT::Discretization>&   dis,
    const int nodeset,
    const std::vector<int> dofs_to_couple
    )
{
  if (!dis->GetCondition("MortarMulti"))
    dserror("this discretization does not have a Mortar-Muti condition");

  Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying =
      Teuchos::rcp(new ADAPTER::CouplingMortar());

  adaptermeshtying->Setup(
      dis,
      dis,
      Teuchos::null,
      dofs_to_couple,
      "MortarMulti",
      dis->Comm(),
      false,
      false,
      nodeset,
      nodeset
  );

  adaptermeshtying->Evaluate();
  p_.push_back(adaptermeshtying->GetPMatrix());
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::CondenseMatrix(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>& mat)
{
  MORTAR::UTILS::MortarMatrixCondensation(mat,p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::CondenseRhs(
    Teuchos::RCP<Epetra_Vector>& rhs)
{
  MORTAR::UTILS::MortarRhsCondensation(rhs,p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::RecoverIncr(
    Teuchos::RCP<Epetra_Vector>& incr)
{
  MORTAR::UTILS::MortarRecover(incr,p_);
}
