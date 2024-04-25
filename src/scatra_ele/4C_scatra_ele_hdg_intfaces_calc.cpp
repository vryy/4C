/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of ScaTraHDG internal faces elements

Integrate internal face terms on an internal faces element

\level 3


*----------------------------------------------------------------------*/
#include "4C_scatra_ele_hdg_intfaces_calc.hpp"

#include "4C_lib_discret_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_ele_action.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGIntFaceImplInterface* DRT::ELEMENTS::ScaTraHDGIntFaceImplInterface::Impl(
    const DRT::Element* ele)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::quad4>::Instance();
    }
    case CORE::FE::CellType::quad8:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::quad8>::Instance();
    }
    case CORE::FE::CellType::quad9:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::quad9>::Instance();
    }
    case CORE::FE::CellType::tri3:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::tri3>::Instance();
    }
    case CORE::FE::CellType::tri6:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::tri6>::Instance();
    }
    case CORE::FE::CellType::line2:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::line2>::Instance();
    }
    case CORE::FE::CellType::line3:
    {
      return ScaTraHDGIntFaceImpl<CORE::FE::CellType::line3>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
  }
  return nullptr;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>*
DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::Instance(CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>>(
            new DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::ScaTraHDGIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::AssembleInternalFacesUsingNeighborData(
    DRT::ELEMENTS::ScaTraHDGIntFace* intface,               ///< internal face element
    std::vector<int>& nds_master,                           ///< nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                            ///< nodal dofset w.r.t. slave element
    Teuchos::ParameterList& params,                         ///< parameter list
    DRT::DiscretizationFaces& discretization,               ///< faces discretization
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector                ///< systemvector
)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ScaTraHDGIntFaceImpl<distype>::EvaluateInternalFaces(
    DRT::ELEMENTS::ScaTraHDGIntFace* intface,  ///< internal face element
    Teuchos::ParameterList& params,            ///< parameter list
    DRT::Discretization& discretization,       ///< discretization
    std::vector<int>& patchlm,                 ///< patch local map
    std::vector<int>& lm_masterToPatch,        ///< local map between master dofs and patchlm
    std::vector<int>& lm_slaveToPatch,         ///< local map between slave dofs and patchlm
    std::vector<int>& lm_faceToPatch,          ///< local map between face dofs and patchlm
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
    std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   ///< element vector blocks
)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
