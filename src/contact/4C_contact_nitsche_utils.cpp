/*---------------------------------------------------------------------*/
/*! \file
\brief Some helpers for nitsche contact

\level 3


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_nitsche_utils.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType parent_distype>
template <int num_dof_per_node>
void Mortar::ElementNitscheData<parent_distype>::AssembleRHS(Mortar::Element* mele,
    const Core::LinAlg::Matrix<Core::FE::num_nodes<parent_distype> * num_dof_per_node, 1>& rhs,
    std::vector<int>& dofs, Teuchos::RCP<Epetra_FEVector> fc) const
{
  const int nen = Core::FE::num_nodes<parent_distype>;

  if (num_dof_per_node * nen > dofs.size())
    FOUR_C_THROW("num_dof_per_node*nen>dofs.size() %d > %d", num_dof_per_node * nen, dofs.size());

  if (fc != Teuchos::null)
  {
    for (int n = 0; n < nen; ++n)
      fc->SumIntoGlobalValues(
          num_dof_per_node, &dofs.at(n * num_dof_per_node), &rhs.A()[n * num_dof_per_node]);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType parent_distype>
template <int num_dof_per_node>
void Mortar::ElementNitscheData<parent_distype>::AssembleMatrix(Mortar::Element* mele,
    const std::unordered_map<int,
        Core::LinAlg::Matrix<Core::FE::num_nodes<parent_distype> * num_dof_per_node, 1>>& k,
    std::vector<int>& dofs, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc) const
{
  const int nen = Core::FE::num_nodes<parent_distype>;

  if (kc != Teuchos::null)
  {
    for (auto& p : k)
    {
      for (int n = 0; n < nen; ++n)
      {
        if (Core::LinAlg::Matrix<num_dof_per_node, 1>(&(p.second.A()[n * num_dof_per_node]), true)
                .NormInf() < 1e-16)
          continue;
        for (int d = 0; d < num_dof_per_node; ++d)
          kc->FEAssemble(
              p.second(n * num_dof_per_node + d), dofs.at(n * num_dof_per_node + d), p.first);
      }
    }
  }
}


template <Core::FE::CellType parent_distype>
void Mortar::ElementNitscheData<parent_distype>::AssembleRHS(
    Mortar::Element* mele, CONTACT::VecBlockType row, Teuchos::RCP<Epetra_FEVector> fc) const
{
  switch (row)
  {
    case CONTACT::VecBlockType::displ:
      AssembleRHS<Core::FE::dim<parent_distype>>(mele, rhs_, mele->MoData().ParentDof(), fc);
      break;
    case CONTACT::VecBlockType::temp:
      if (mele->MoData().ParentTempDof().size())
        AssembleRHS<1>(mele, tsi_data_.rhs_t_, mele->MoData().ParentTempDof(), fc);
      break;
    case CONTACT::VecBlockType::porofluid:
      if (mele->MoData().ParentPFDof().size())  // not if the parent is an impermeable element
        AssembleRHS<Core::FE::dim<parent_distype> + 1>(
            mele, poro_data_.rhs_p_, mele->MoData().ParentPFDof(), fc);
      break;
    case CONTACT::VecBlockType::scatra:
      if (mele->MoData().ParentScalarDof().size())
        AssembleRHS<1>(mele, ssi_data_.rhs_s_, mele->MoData().ParentScalarDof(), fc);
      break;
    case CONTACT::VecBlockType::elch:
      if (mele->MoData().ParentScalarDof().size())
        AssembleRHS<2>(mele, ssi_elch_data_.rhs_e_, mele->MoData().ParentScalarDof(), fc);
      break;
    default:
      FOUR_C_THROW("unknown row");
  }
}

template <Core::FE::CellType parent_distype>
void Mortar::ElementNitscheData<parent_distype>::AssembleMatrix(Mortar::Element* mele,
    CONTACT::MatBlockType block, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc) const
{
  switch (block)
  {
    case CONTACT::MatBlockType::displ_displ:
      AssembleMatrix<Core::FE::dim<parent_distype>>(mele, k_, mele->MoData().ParentDof(), kc);
      break;
    case CONTACT::MatBlockType::displ_temp:
      AssembleMatrix<Core::FE::dim<parent_distype>>(
          mele, tsi_data_.k_dt_, mele->MoData().ParentDof(), kc);
      break;
    case CONTACT::MatBlockType::temp_displ:
      if (mele->MoData().ParentTempDof().size())
        AssembleMatrix<1>(mele, tsi_data_.k_td_, mele->MoData().ParentTempDof(), kc);
      break;
    case CONTACT::MatBlockType::temp_temp:
      if (mele->MoData().ParentTempDof().size())
        AssembleMatrix<1>(mele, tsi_data_.k_tt_, mele->MoData().ParentTempDof(), kc);
      break;
    case CONTACT::MatBlockType::displ_porofluid:
      AssembleMatrix<Core::FE::dim<parent_distype>>(
          mele, poro_data_.k_dp_, mele->MoData().ParentDof(), kc);
      break;
    case CONTACT::MatBlockType::porofluid_displ:
      if (mele->MoData().ParentPFDof().size())  // not if the parent is an impermeable element
        AssembleMatrix<Core::FE::dim<parent_distype> + 1>(
            mele, poro_data_.k_pd_, mele->MoData().ParentPFDof(), kc);
      break;
    case CONTACT::MatBlockType::porofluid_porofluid:
      if (mele->MoData().ParentPFDof().size())  // not if the parent is an impermeable element
        AssembleMatrix<Core::FE::dim<parent_distype> + 1>(
            mele, poro_data_.k_pp_, mele->MoData().ParentPFDof(), kc);
      break;
    case CONTACT::MatBlockType::displ_scatra:
      AssembleMatrix<Core::FE::dim<parent_distype>>(
          mele, ssi_data_.k_ds_, mele->MoData().ParentDof(), kc);
      break;
    case CONTACT::MatBlockType::scatra_displ:
      if (mele->MoData().ParentScalarDof().size())
        AssembleMatrix<1>(mele, ssi_data_.k_sd_, mele->MoData().ParentScalarDof(), kc);
      break;
    case CONTACT::MatBlockType::scatra_scatra:
      if (mele->MoData().ParentScalarDof().size())
        AssembleMatrix<1>(mele, ssi_data_.k_ss_, mele->MoData().ParentScalarDof(), kc);
      break;
    case CONTACT::MatBlockType::displ_elch:
      AssembleMatrix<Core::FE::dim<parent_distype>>(
          mele, ssi_elch_data_.k_de_, mele->MoData().ParentDof(), kc);
      break;
    case CONTACT::MatBlockType::elch_displ:
      if (mele->MoData().ParentScalarDof().size())
        AssembleMatrix<2>(mele, ssi_elch_data_.k_ed_, mele->MoData().ParentScalarDof(), kc);
      break;
    case CONTACT::MatBlockType::elch_elch:
      if (mele->MoData().ParentScalarDof().size())
        AssembleMatrix<2>(mele, ssi_elch_data_.k_ee_, mele->MoData().ParentScalarDof(), kc);
      break;
    default:
      FOUR_C_THROW("unknown matrix block");
      break;
  }
}


template class Mortar::ElementNitscheData<Core::FE::CellType::hex8>;
template class Mortar::ElementNitscheData<Core::FE::CellType::tet4>;
template class Mortar::ElementNitscheData<Core::FE::CellType::hex27>;
template class Mortar::ElementNitscheData<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
