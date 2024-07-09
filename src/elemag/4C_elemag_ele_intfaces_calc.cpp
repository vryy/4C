/*----------------------------------------------------------------------*/
/*! \file

\brief Integrate internal face terms on an internal faces element

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_elemag_ele_intfaces_calc.hpp"

#include "4C_elemag_ele_action.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagIntFaceImplInterface* Discret::ELEMENTS::ElemagIntFaceImplInterface::impl(
    const Core::Elements::Element* ele)
{
  switch (ele->shape())
  {
    case Core::FE::CellType::quad4:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::quad4>::instance();
    }
    case Core::FE::CellType::quad8:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::quad8>::instance();
    }
    case Core::FE::CellType::quad9:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::quad9>::instance();
    }
    case Core::FE::CellType::tri3:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::tri3>::instance();
    }
    case Core::FE::CellType::tri6:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::tri6>::instance();
    }
    case Core::FE::CellType::line2:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::line2>::instance();
    }
    case Core::FE::CellType::line3:
    {
      return ElemagIntFaceImpl<Core::FE::CellType::line3>::instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->shape(), ele->num_node());
      break;
  }
  return nullptr;
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::ElemagIntFaceImpl<distype>*
Discret::ELEMENTS::ElemagIntFaceImpl<distype>::instance(Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::ElemagIntFaceImpl<distype>>(
            new Discret::ELEMENTS::ElemagIntFaceImpl<distype>());
      });

  return singleton_owner.instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ElemagIntFaceImpl<distype>::ElemagIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ElemagIntFaceImpl<distype>::assemble_internal_faces_using_neighbor_data(
    Discret::ELEMENTS::ElemagIntFace* intface,              // internal face element
    std::vector<int>& nds_master,                           // nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                            // nodal dofset w.r.t. slave element
    Teuchos::ParameterList& params,                         // parameter list
    Core::FE::DiscretizationFaces& discretization,          // faces discretization
    Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix,  // systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector                // systemvector
)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ElemagIntFaceImpl<distype>::evaluate_internal_faces(
    Discret::ELEMENTS::ElemagIntFace* intface,  // internal face element
    Teuchos::ParameterList& params,             // parameter list
    Core::FE::Discretization& discretization,   // discretization
    std::vector<int>& patchlm,                  // patch local map
    std::vector<int>& lm_masterToPatch,         // local map between master dofs and patchlm
    std::vector<int>& lm_slaveToPatch,          // local map between slave dofs and patchlm
    std::vector<int>& lm_faceToPatch,           // local map between face dofs and patchlm
    std::vector<int>& lm_masterNodeToPatch,     // local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,      // local map between slave nodes and nodes in patch
    std::vector<Core::LinAlg::SerialDenseMatrix>& elemat_blocks,  // element matrix blocks
    std::vector<Core::LinAlg::SerialDenseVector>& elevec_blocks   // element vector blocks
)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
