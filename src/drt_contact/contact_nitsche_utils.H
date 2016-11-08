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


namespace MORTAR
{
class MortarElementNitscheContainer
{
public:
  MortarElementNitscheContainer(){};

  virtual void Clear() =0;

  virtual void Assemble(
      MORTAR::MortarElement* mele,
      Teuchos::RCP<Epetra_FEVector> fc,
      Teuchos::RCP<LINALG::SparseMatrix> kc) = 0;

  virtual double* rhs(int dof)=0;
  virtual double* rhs()=0;
  virtual double* k(int col) =0;
  virtual double* k(int col, int dof)=0;

};

template<int num_parent_dof>
class MortarElementNitscheData : public MortarElementNitscheContainer
{
public:
  virtual double* rhs(int dof) {return &rhs_(dof);}
  virtual double* rhs() {return rhs_.A();}
  virtual double* k(int col) {return k_[col].A();}
  virtual double* k(int col,int dof){return &k_[col](dof);}
  virtual void Assemble(
      MORTAR::MortarElement* mele,
      Teuchos::RCP<Epetra_FEVector> fc,
      Teuchos::RCP<LINALG::SparseMatrix> kc);

  virtual void Clear(){rhs_.Clear();k_.clear();}


private:
  LINALG::Matrix<num_parent_dof,1> rhs_;
  std::map<int,LINALG::Matrix<num_parent_dof,1> > k_;
};

}

#endif // CONTACT_NITSCHE_UTILS_H
