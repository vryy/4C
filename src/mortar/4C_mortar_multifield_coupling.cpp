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
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_mortar_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
void Mortar::MultiFieldCoupling::push_back_coupling(
    const Teuchos::RCP<Core::FE::Discretization>& dis, const int nodeset,
    const std::vector<int>& dofs_to_couple, const Teuchos::ParameterList& mortar_params,
    const Teuchos::ParameterList& contact_params, const Teuchos::ParameterList& binning_params,
    const std::map<std::string, Teuchos::RCP<Core::FE::Discretization>>& discretization_map,
    const Core::UTILS::FunctionManager& function_manager,
    Teuchos::RCP<Core::IO::OutputControl> output_control,
    const Core::FE::ShapeFunctionType shape_function_type, const int ndim)
{
  if (!dis->get_condition("MortarMulti"))
    FOUR_C_THROW("this discretization does not have a Mortar-Muti condition");

  Teuchos::RCP<Coupling::Adapter::CouplingMortar> adaptermeshtying =
      Teuchos::rcp(new Coupling::Adapter::CouplingMortar(
          ndim, mortar_params, contact_params, shape_function_type));

  adaptermeshtying->setup(dis, dis, Teuchos::null, dofs_to_couple, "MortarMulti", dis->get_comm(),
      function_manager, binning_params, discretization_map, output_control, shape_function_type,
      false, false, nodeset, nodeset);

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
