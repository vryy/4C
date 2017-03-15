/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_utils.cpp

\brief Some helpers for nitsche contact

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#include "contact_nitsche_utils.H"
#include "../drt_mortar/mortar_element.H"
#include <Epetra_FECrsMatrix.h>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType parent_distype>
void MORTAR::MortarElementNitscheData<parent_distype>::Assemble(
    MORTAR::MortarElement* mele,
    Teuchos::RCP<Epetra_FEVector> fc,
    Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  const int num_parent_disp_dof =
      DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement
      *
      DRT::UTILS::DisTypeToDim<parent_distype>::dim;
  if (fc!=Teuchos::null)
  fc->SumIntoGlobalValues(num_parent_disp_dof,&mele->MoData().ParentDof()[0],rhs_.A());

  if (kc!=Teuchos::null)
    for (typename std::map<int,LINALG::Matrix<num_parent_disp_dof,1> >::const_iterator
        p=k_.begin();p!=k_.end();++p)
      for (int dof=0;dof<(int)mele->MoData().ParentDof().size();++dof)
        kc->FEAssemble(p->second(dof),mele->MoData().ParentDof()[dof],p->first);
}

template class MORTAR::MortarElementNitscheData<DRT::Element::hex8>;
template class MORTAR::MortarElementNitscheData<DRT::Element::tet4>;
