/*----------------------------------------------------------------------*/
/*! \file

\brief Integrate internal face terms on an internal faces element

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_elemag_ele_intfaces_calc.hpp"

#include "4C_elemag_ele_action.hpp"
#include "4C_lib_discret_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ElemagIntFaceImplInterface* DRT::ELEMENTS::ElemagIntFaceImplInterface::Impl(
    const CORE::Elements::Element* ele)
{
  switch (ele->Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::quad4>::Instance();
    }
    case CORE::FE::CellType::quad8:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::quad8>::Instance();
    }
    case CORE::FE::CellType::quad9:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::quad9>::Instance();
    }
    case CORE::FE::CellType::tri3:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::tri3>::Instance();
    }
    case CORE::FE::CellType::tri6:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::tri6>::Instance();
    }
    case CORE::FE::CellType::line2:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::line2>::Instance();
    }
    case CORE::FE::CellType::line3:
    {
      return ElemagIntFaceImpl<CORE::FE::CellType::line3>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagIntFaceImpl<distype>* DRT::ELEMENTS::ElemagIntFaceImpl<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::ElemagIntFaceImpl<distype>>(
            new DRT::ELEMENTS::ElemagIntFaceImpl<distype>());
      });

  return singleton_owner.Instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagIntFaceImpl<distype>::ElemagIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagIntFaceImpl<distype>::assemble_internal_faces_using_neighbor_data(
    DRT::ELEMENTS::ElemagIntFace* intface,                  // internal face element
    std::vector<int>& nds_master,                           // nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                            // nodal dofset w.r.t. slave element
    Teuchos::ParameterList& params,                         // parameter list
    DRT::DiscretizationFaces& discretization,               // faces discretization
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systemmatrix,  // systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector                // systemvector
)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagIntFaceImpl<distype>::evaluate_internal_faces(
    DRT::ELEMENTS::ElemagIntFace* intface,   // internal face element
    Teuchos::ParameterList& params,          // parameter list
    DRT::Discretization& discretization,     // discretization
    std::vector<int>& patchlm,               // patch local map
    std::vector<int>& lm_masterToPatch,      // local map between master dofs and patchlm
    std::vector<int>& lm_slaveToPatch,       // local map between slave dofs and patchlm
    std::vector<int>& lm_faceToPatch,        // local map between face dofs and patchlm
    std::vector<int>& lm_masterNodeToPatch,  // local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   // local map between slave nodes and nodes in patch
    std::vector<CORE::LINALG::SerialDenseMatrix>& elemat_blocks,  // element matrix blocks
    std::vector<CORE::LINALG::SerialDenseVector>& elevec_blocks   // element vector blocks
)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
