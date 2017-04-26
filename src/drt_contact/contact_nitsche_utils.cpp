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
template<int num_dof_per_node>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleRHS(
    MORTAR::MortarElement* mele,
    const LINALG::Matrix<
      DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node,1>& rhs,
    Teuchos::RCP<Epetra_FEVector> fc)
{
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement;
  const int nsd = DRT::UTILS::DisTypeToDim<parent_distype>::dim;

  if (num_dof_per_node>nsd)
    dserror("num_dof_per_node > nsd");

  if (fc!=Teuchos::null)
    for (int n=0;n<nen;++n)
      fc->SumIntoGlobalValues(num_dof_per_node,&mele->MoData().ParentDof().at(n*nsd),&rhs.A()[n*num_dof_per_node]);

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType parent_distype>
template<int num_dof_per_node>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleMatrix(
    MORTAR::MortarElement* mele,
    const std::map<int,LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node,1> >& k,
    Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement;
  const int nsd = DRT::UTILS::DisTypeToDim<parent_distype>::dim;

  if (kc!=Teuchos::null)
    for (typename std::map<
        int,LINALG::Matrix<nen*num_dof_per_node,1>
        >::const_iterator p=k.begin();p!=k.end();++p)
      for (int n=0;n<nen;++n)
        for (int d=0;d<num_dof_per_node;++d)
          kc->FEAssemble(p->second(n*num_dof_per_node+d),mele->MoData().ParentDof()[n*nsd+d],p->first);
}


template<DRT::Element::DiscretizationType parent_distype>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleRHS(
    MORTAR::MortarElement* mele,
    DRT::UTILS::VecBlockType row,
    Teuchos::RCP<Epetra_FEVector> fc)
{
  switch(row)
  {
  case DRT::UTILS::block_displ:
    AssembleRHS<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,rhs_,fc);
    break;
  default: dserror("unknown row");
  }
}

template<DRT::Element::DiscretizationType parent_distype>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleMatrix(
    MORTAR::MortarElement* mele,
    DRT::UTILS::MatBlockType block,
    Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  switch(block)
  {
  case DRT::UTILS::block_displ_displ:
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,k_,kc);
    break;
  default: dserror("unknown matrix block"); break;
  }
}


template class MORTAR::MortarElementNitscheData<DRT::Element::hex8>;
template class MORTAR::MortarElementNitscheData<DRT::Element::tet4>;
