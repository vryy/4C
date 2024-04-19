/*-----------------------------------------------------------------------*/
/*! \file
\brief Class performing coupling (condensation/recovery) for dual mortar
       methods in (volume) monolithic multi-physics applications, i.e. in
       block matrix systems. This also accounts for the correct condensation
       in the off-diagonal matrix blocks

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_mortar_multifield_coupling.hpp"

#include "baci_coupling_adapter_mortar.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_blocksparsematrix.hpp"
#include "baci_mortar_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::PushBackCoupling(const Teuchos::RCP<DRT::Discretization>& dis,
    const int nodeset, const std::vector<int> dofs_to_couple)
{
  if (!dis->GetCondition("MortarMulti"))
    FOUR_C_THROW("this discretization does not have a Mortar-Muti condition");

  Teuchos::RCP<CORE::ADAPTER::CouplingMortar> adaptermeshtying =
      Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(GLOBAL::Problem::Instance()->NDim(),
          GLOBAL::Problem::Instance()->MortarCouplingParams(),
          GLOBAL::Problem::Instance()->ContactDynamicParams(),
          GLOBAL::Problem::Instance()->SpatialApproximationType()));

  adaptermeshtying->Setup(dis, dis, Teuchos::null, dofs_to_couple, "MortarMulti", dis->Comm(),
      false, false, nodeset, nodeset);

  adaptermeshtying->Evaluate();
  p_.push_back(adaptermeshtying->GetMortarMatrixP());
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::CondenseMatrix(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>& mat)
{
  MORTAR::UTILS::MortarMatrixCondensation(mat, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::CondenseRhs(Teuchos::RCP<Epetra_Vector>& rhs)
{
  MORTAR::UTILS::MortarRhsCondensation(rhs, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void MORTAR::MultiFieldCoupling::RecoverIncr(Teuchos::RCP<Epetra_Vector>& incr)
{
  MORTAR::UTILS::MortarRecover(incr, p_);
}

FOUR_C_NAMESPACE_CLOSE
