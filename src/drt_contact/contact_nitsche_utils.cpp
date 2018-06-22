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
template<int num_dof_per_node,int num_dof_per_node_to_assemble>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleRHS(
    MORTAR::MortarElement* mele,
    const LINALG::Matrix<
    DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node_to_assemble,1>& rhs,
    std::vector<int>& dofs,
    Teuchos::RCP<Epetra_FEVector> fc)
{
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement;
  //const int nsd = DRT::UTILS::DisTypeToDim<parent_distype>::dim;

  if (num_dof_per_node_to_assemble*nen>dofs.size())
    dserror("num_dof_per_node_to_assemble*nen>dofs.size() %d > %d",num_dof_per_node*nen,dofs.size());

  if (fc!=Teuchos::null)
    for (int n=0;n<nen;++n)
      fc->SumIntoGlobalValues(num_dof_per_node_to_assemble,&dofs.at(n*num_dof_per_node),&rhs.A()[n*num_dof_per_node_to_assemble]);

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType parent_distype>
template<int num_dof_per_node,int num_dof_per_node_to_assemble>
void MORTAR::MortarElementNitscheData<parent_distype>::AssembleMatrix(
    MORTAR::MortarElement* mele,
    const std::unordered_map<int,LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node_to_assemble,1> >& k,
    std::vector<int>& dofs,
    Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement;
  //const int nsd = DRT::UTILS::DisTypeToDim<parent_distype>::dim;

  if (kc!=Teuchos::null)
    for (typename std::unordered_map<
        int,LINALG::Matrix<nen*num_dof_per_node_to_assemble,1>
        >::const_iterator p=k.begin();p!=k.end();++p)
      for (int n=0;n<nen;++n)
        for (int d=0;d<num_dof_per_node;++d)
          if (fabs(p->second(n*num_dof_per_node_to_assemble+d)) > 1e-16) kc->FEAssemble(p->second(n*num_dof_per_node_to_assemble+d),dofs.at(n*num_dof_per_node+d),p->first);
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
    AssembleRHS<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,rhs_,mele->MoData().ParentDof(),fc);
    break;
  case DRT::UTILS::block_temp:
    AssembleRHS<DRT::UTILS::DisTypeToDim<parent_distype>::dim,1>(mele,tsi_data_.rhs_t_,mele->MoData().ParentDof(),fc);
    break;
  case DRT::UTILS::block_porofluid:
  {
    if (mele->MoData().ParentPFDof().size()) //not if the parent is an impermeable element
      AssembleRHS<DRT::UTILS::DisTypeToDim<parent_distype>::dim+1>(mele,poro_data_.rhs_p_,mele->MoData().ParentPFDof(),fc);
    break;
  }
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
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,k_,mele->MoData().ParentDof(),kc);
    break;
  case DRT::UTILS::block_displ_temp:
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,tsi_data_.k_dt_,mele->MoData().ParentDof(),kc);
    break;
  case DRT::UTILS::block_temp_displ:
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim,1>(mele,tsi_data_.k_td_,mele->MoData().ParentDof(),kc);
    break;
  case DRT::UTILS::block_temp_temp:
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim,1>(mele,tsi_data_.k_tt_,mele->MoData().ParentDof(),kc);
  case DRT::UTILS::block_displ_porofluid:
    AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim>(mele,poro_data_.k_dp_,mele->MoData().ParentDof(),kc);
    break;
  case DRT::UTILS::block_porofluid_displ:
    if (mele->MoData().ParentPFDof().size()) //not if the parent is an impermeable element
      AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim+1>(mele,poro_data_.k_pd_,mele->MoData().ParentPFDof(),kc);
    break;
  case DRT::UTILS::block_porofluid_porofluid:
    if (mele->MoData().ParentPFDof().size()) //not if the parent is an impermeable element
      AssembleMatrix<DRT::UTILS::DisTypeToDim<parent_distype>::dim+1>(mele,poro_data_.k_pp_,mele->MoData().ParentPFDof(),kc);
  break;
  default: dserror("unknown matrix block"); break;
  }
}


template class MORTAR::MortarElementNitscheData<DRT::Element::hex8>;
template class MORTAR::MortarElementNitscheData<DRT::Element::tet4>;
template class MORTAR::MortarElementNitscheData<DRT::Element::hex27>;
template class MORTAR::MortarElementNitscheData<DRT::Element::nurbs27>;
