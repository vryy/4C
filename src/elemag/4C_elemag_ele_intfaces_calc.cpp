// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_elemag_ele_intfaces_calc.hpp"

#include "4C_elemag_ele_action.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::Elements::ElemagIntFaceImplInterface* Discret::Elements::ElemagIntFaceImplInterface::impl(
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
Discret::Elements::ElemagIntFaceImpl<distype>*
Discret::Elements::ElemagIntFaceImpl<distype>::instance(Core::Utils::SingletonAction action)
{
  static auto singleton_owner = Core::Utils::make_singleton_owner(
      []()
      {
        return std::unique_ptr<Discret::Elements::ElemagIntFaceImpl<distype>>(
            new Discret::Elements::ElemagIntFaceImpl<distype>());
      });

  return singleton_owner.instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::ElemagIntFaceImpl<distype>::ElemagIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::ElemagIntFaceImpl<distype>::assemble_internal_faces_using_neighbor_data(
    Discret::Elements::ElemagIntFace* intface,                 // internal face element
    std::vector<int>& nds_master,                              // nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                               // nodal dofset w.r.t. slave element
    Teuchos::ParameterList& params,                            // parameter list
    Core::FE::DiscretizationFaces& discretization,             // faces discretization
    std::shared_ptr<Core::LinAlg::SparseMatrix> systemmatrix,  // systemmatrix
    std::shared_ptr<Core::LinAlg::Vector<double>> systemvector  // systemvector
)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::Elements::ElemagIntFaceImpl<distype>::evaluate_internal_faces(
    Discret::Elements::ElemagIntFace* intface,  // internal face element
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
