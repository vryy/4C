/*----------------------------------------------------------------------*/
/*! \file

\brief Integrate internal face terms on an internal faces element

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_elemag_ele_intfaces_calc.H"

#include "baci_elemag_ele_action.H"
#include "baci_lib_discret_faces.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ElemagIntFaceImplInterface* DRT::ELEMENTS::ElemagIntFaceImplInterface::Impl(
    const DRT::Element* ele)
{
  switch (ele->Shape())
  {
    case DRT::Element::quad4:
    {
      return ElemagIntFaceImpl<DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return ElemagIntFaceImpl<DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return ElemagIntFaceImpl<DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return ElemagIntFaceImpl<DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return ElemagIntFaceImpl<DRT::Element::tri6>::Instance();
    }
    case DRT::Element::line2:
    {
      return ElemagIntFaceImpl<DRT::Element::line2>::Instance();
    }
    case DRT::Element::line3:
    {
      return ElemagIntFaceImpl<DRT::Element::line3>::Instance();
    }
    default:
      dserror(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
      break;
  }
  return NULL;
}

template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ElemagIntFaceImpl<distype>::ElemagIntFaceImpl()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagIntFaceImpl<distype>::AssembleInternalFacesUsingNeighborData(
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
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagIntFaceImpl<distype>::EvaluateInternalFaces(
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
