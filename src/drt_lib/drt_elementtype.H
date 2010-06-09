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

  class Element;
  class Discretization;

class ElementType : public ParObjectType
{
protected:

  // only derived classes might create an instance
  ElementType();

public:

  /// test equality by comparing pointers as there are just singletons allowed
  bool operator==( const ElementType & other ) const
  { return this == &other; }

  bool operator!=( const ElementType & other ) const
  { return this != &other; }

  virtual int Initialize(DRT::Discretization& dis);

  virtual inline void PreEvaluate(DRT::Discretization& dis,
                                  Teuchos::ParameterList& p,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                  Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector1,
                                  Teuchos::RCP<Epetra_Vector>          systemvector2,
                                  Teuchos::RCP<Epetra_Vector>          systemvector3)
  { return; }

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) = 0;

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) = 0;

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions ) {}

};

}

#endif
