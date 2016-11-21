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
template<int num_parent_dof>
void MORTAR::MortarElementNitscheData<num_parent_dof>::Assemble(
    MORTAR::MortarElement* mele,
    Teuchos::RCP<Epetra_FEVector> fc,
    Teuchos::RCP<LINALG::SparseMatrix> kc)
{
  if (fc!=Teuchos::null)
  fc->SumIntoGlobalValues(num_parent_dof,&mele->MoData().ParentDof()[0],rhs_.A());

  if (kc!=Teuchos::null)
    for (typename std::map<int,LINALG::Matrix<num_parent_dof,1> >::const_iterator
        p=k_.begin();p!=k_.end();++p)
      for (int dof=0;dof<(int)mele->MoData().ParentDof().size();++dof)
        kc->FEAssemble(p->second(dof),mele->MoData().ParentDof()[dof],p->first);
}

template class MORTAR::MortarElementNitscheData<8>;
template class MORTAR::MortarElementNitscheData<12>;
template class MORTAR::MortarElementNitscheData<24>;
