/*-----------------------------------------------------------------------*/
/*! \file
\brief Class performing coupling (condensation/recovery) for dual mortar
       methods in (volume) monolithic multi-physics applications, i.e. in
       block matrix systems. This also accounts for the correct condensation
       in the off-diagonal matrix blocks

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_mortar_multifield_coupling.hpp"

#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_mortar_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::PushBackCoupling(const Teuchos::RCP<Discret::Discretization>& dis,
    const int nodeset, const std::vector<int> dofs_to_couple)
{
  if (!dis->GetCondition("MortarMulti"))
    FOUR_C_THROW("this discretization does not have a Mortar-Muti condition");

  Teuchos::RCP<Core::Adapter::CouplingMortar> adaptermeshtying =
      Teuchos::rcp(new Core::Adapter::CouplingMortar(Global::Problem::Instance()->NDim(),
          Global::Problem::Instance()->mortar_coupling_params(),
          Global::Problem::Instance()->contact_dynamic_params(),
          Global::Problem::Instance()->spatial_approximation_type()));

  adaptermeshtying->Setup(dis, dis, Teuchos::null, dofs_to_couple, "MortarMulti", dis->Comm(),
      Global::Problem::Instance()->FunctionManager(), false, false, nodeset, nodeset);

  adaptermeshtying->Evaluate();
  p_.push_back(adaptermeshtying->GetMortarMatrixP());
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::CondenseMatrix(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& mat)
{
  Mortar::UTILS::MortarMatrixCondensation(mat, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::CondenseRhs(Teuchos::RCP<Epetra_Vector>& rhs)
{
  Mortar::UTILS::MortarRhsCondensation(rhs, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::RecoverIncr(Teuchos::RCP<Epetra_Vector>& incr)
{
  Mortar::UTILS::MortarRecover(incr, p_);
}

FOUR_C_NAMESPACE_CLOSE
