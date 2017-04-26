/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_utils.H

\brief Some helpers for nitsche contact

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_NITSCHE_UTILS_H
#define CONTACT_NITSCHE_UTILS_H

#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"


namespace MORTAR
{
class MortarElementNitscheContainer
{
public:
  MortarElementNitscheContainer(){};

  virtual void Clear() =0;

  virtual void AssembleRHS(
      MORTAR::MortarElement* mele,
      DRT::UTILS::VecBlockType row,
      Teuchos::RCP<Epetra_FEVector> fc) = 0;

  virtual void AssembleMatrix(
      MORTAR::MortarElement* mele,
      DRT::UTILS::MatBlockType block,
      Teuchos::RCP<LINALG::SparseMatrix> kc) = 0;

  virtual double* rhs(int dof)=0;
  virtual double* rhs()=0;
  virtual double* k(int col) =0;
  virtual double* k(int col, int dof)=0;

};

template<DRT::Element::DiscretizationType parent_distype>
class MortarElementNitscheData : public MortarElementNitscheContainer
{
public:
  virtual double* rhs(int dof) {return &rhs_(dof);}
  virtual double* rhs() {return rhs_.A();}
  virtual double* k(int col) {return k_[col].A();}
  virtual double* k(int col,int dof){return &k_[col](dof);}

  virtual void AssembleRHS(
      MORTAR::MortarElement* mele,
      DRT::UTILS::VecBlockType row,
      Teuchos::RCP<Epetra_FEVector> fc);

  virtual void AssembleMatrix(
      MORTAR::MortarElement* mele,
      DRT::UTILS::MatBlockType block,
      Teuchos::RCP<LINALG::SparseMatrix> kc);

  template<int num_dof_per_node>
  void AssembleRHS(
      MORTAR::MortarElement* mele,
      const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node,1>& rhs,
      Teuchos::RCP<Epetra_FEVector> fc);

  template<int num_dof_per_node>
  void AssembleMatrix(
      MORTAR::MortarElement* mele,
      const std::map<int,LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement*num_dof_per_node,1> >& k,
      Teuchos::RCP<LINALG::SparseMatrix> kc);

  virtual void Clear(){rhs_.Clear();k_.clear();}


private:
  LINALG::Matrix<
  DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement
  *
  DRT::UTILS::DisTypeToDim<parent_distype>::dim,1> rhs_;
  std::map<int,LINALG::Matrix<
  DRT::UTILS::DisTypeToNumNodePerEle<parent_distype>::numNodePerElement
    *
    DRT::UTILS::DisTypeToDim<parent_distype>::dim,1> > k_;
};

}

#endif // CONTACT_NITSCHE_UTILS_H
