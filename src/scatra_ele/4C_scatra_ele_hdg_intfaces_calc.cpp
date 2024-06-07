/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of ScaTraHDG internal faces elements

Integrate internal face terms on an internal faces element

\level 3


*----------------------------------------------------------------------*/
#include "4C_scatra_ele_hdg_intfaces_calc.hpp"

#include "4C_fem_discretization_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_ele_action.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ScaTraHDGIntFaceImplInterface*
Discret::ELEMENTS::ScaTraHDGIntFaceImplInterface::Impl(const Core::Elements::Element* ele)
{
  switch (ele->Shape())
  {
    case Core::FE::CellType::quad4:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::quad4>::Instance();
    }
    case Core::FE::CellType::quad8:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::quad8>::Instance();
    }
    case Core::FE::CellType::quad9:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::quad9>::Instance();
    }
    case Core::FE::CellType::tri3:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::tri3>::Instance();
    }
    case Core::FE::CellType::tri6:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::tri6>::Instance();
    }
    case Core::FE::CellType::line2:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::line2>::Instance();
    }
    case Core::FE::CellType::line3:
    {
      return ScaTraHDGIntFaceImpl<Core::FE::CellType::line3>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>*
Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::Instance(Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>>(
            new Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::ScaTraHDGIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::assemble_internal_faces_using_neighbor_data(
    Discret::ELEMENTS::ScaTraHDGIntFace* intface,           ///< internal face element
    std::vector<int>& nds_master,                           ///< nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                            ///< nodal dofset w.r.t. slave element
    Teuchos::ParameterList& params,                         ///< parameter list
    Discret::DiscretizationFaces& discretization,           ///< faces discretization
    Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::evaluate_internal_faces(
    Discret::ELEMENTS::ScaTraHDGIntFace* intface,  ///< internal face element
    Teuchos::ParameterList& params,                ///< parameter list
    Discret::Discretization& discretization,       ///< discretization
    std::vector<int>& patchlm,                     ///< patch local map
    std::vector<int>& lm_masterToPatch,            ///< local map between master dofs and patchlm
    std::vector<int>& lm_slaveToPatch,             ///< local map between slave dofs and patchlm
    std::vector<int>& lm_faceToPatch,              ///< local map between face dofs and patchlm
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    std::vector<Core::LinAlg::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
    std::vector<Core::LinAlg::SerialDenseVector>& elevec_blocks   ///< element vector blocks
)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
