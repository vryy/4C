#ifndef DRT_ELEMENTTPYE_H
#define DRT_ELEMENTTPYE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

#include "drt_parobjectfactory.H"
#include "drt_element.H"

namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{

  class Discretization;

class ElementObjectType : public ParObjectType
{
protected:

  // only derived classes might create an instance
  ElementObjectType();

public:

  /// test equality by comparing pointers as there are just singletons allowed
  bool operator==( const ElementObjectType & other ) const
  { return this == &other; }

  virtual int Initialize(DRT::Discretization& dis);

  virtual inline void PreEvaluate(DRT::Discretization& dis,
                                  Teuchos::ParameterList& p,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector1,
                                  Teuchos::RCP<Epetra_Vector>          systemvector2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector3)
  { return; }

};

}

#endif
