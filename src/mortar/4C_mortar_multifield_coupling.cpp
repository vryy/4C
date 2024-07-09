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
void Mortar::MultiFieldCoupling::push_back_coupling(
    const Teuchos::RCP<Core::FE::Discretization>& dis, const int nodeset,
    const std::vector<int> dofs_to_couple)
{
  if (!dis->get_condition("MortarMulti"))
    FOUR_C_THROW("this discretization does not have a Mortar-Muti condition");

  Teuchos::RCP<Core::Adapter::CouplingMortar> adaptermeshtying =
      Teuchos::rcp(new Core::Adapter::CouplingMortar(Global::Problem::instance()->n_dim(),
          Global::Problem::instance()->mortar_coupling_params(),
          Global::Problem::instance()->contact_dynamic_params(),
          Global::Problem::instance()->spatial_approximation_type()));

  adaptermeshtying->setup(dis, dis, Teuchos::null, dofs_to_couple, "MortarMulti", dis->get_comm(),
      Global::Problem::instance()->function_manager(), false, false, nodeset, nodeset);

  adaptermeshtying->evaluate();
  p_.push_back(adaptermeshtying->get_mortar_matrix_p());
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::condense_matrix(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& mat)
{
  Mortar::UTILS::MortarMatrixCondensation(mat, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::condense_rhs(Teuchos::RCP<Epetra_Vector>& rhs)
{
  Mortar::UTILS::MortarRhsCondensation(rhs, p_);
}

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::recover_incr(Teuchos::RCP<Epetra_Vector>& incr)
{
  Mortar::UTILS::MortarRecover(incr, p_);
}

FOUR_C_NAMESPACE_CLOSE
