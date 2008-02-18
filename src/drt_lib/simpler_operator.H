/*!----------------------------------------------------------------------
\file simpler_operator.H

\class LINALG::SIMPLER_Operator

\brief An approximate block factorization preconditioner based on the 
       SIMPLE family of methods

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef SIMPLER_OPERATOR_H
#define SIMPLER_OPERATOR_H

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_MultiLevelPreconditioner.h"

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_systemmatrix.H"
#include "linalg_solver.H"
#include "linalg_utils.H"

using namespace std;
using namespace Teuchos;

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class AMG_Operator

\brief A multipurpose experimental multigrid operator

This operator based on the ml advanced programming interface is a
multipurpose development object for amg ideas that shall be tested in
the baci framework

\author gee (gee@lnm.mw.tum.de)
*/
class  SIMPLER_Operator : public virtual Epetra_Operator
{
public:

  /*!
  \brief Standard Constructor
  */
  explicit SIMPLER_Operator(RCP<Epetra_CrsMatrix> A,
                            const ParameterList& velocitylist,
                            const ParameterList& pressurelist,
                            FILE* outfile = NULL);

  /*!
  \brief Destructor
  */
  virtual ~SIMPLER_Operator() {}

  /*!
  \brief Returns the label of this class.
  */
  const char* Label() const { return "LINALG::SIMPLER_Operator"; }

  /*!
  \brief get Comm of this class
  
  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class 
  
  */
  const Epetra_Comm& Comm() const {return(A_->Comm());}


  /*!
  \brief Get fine level OperatorDomainMap
  
  Derived from Epetra_Operator, get fine level OperatorDomainMap
  
  */
  const Epetra_Map& OperatorDomainMap() const {return A_->FullDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap
  
  Derived from Epetra_Operator, get fine level OperatorRangeMap
  
  */
  const Epetra_Map& OperatorRangeMap() const {return A_->FullRangeMap();}
  
  /*!
  \brief ApplyInverse the preconditioner
  
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator. 
   
   
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  /*!
  \brief not implemented
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { 
    dserror("Apply does not make sense for LINALG::SIMPLER_Operator");
    return(-1);
  }
    
  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose) 
  { 
    dserror("SetUseTranspose not impl.");
    return -1;
  }

  /*!
  \brief not implemented
  */
  double NormInf() const 
  {
    dserror("NormInf not impl.");
    return(-1.0);
  }


  /*!
  \brief not implemented
  */
  bool UseTranspose() const 
  {
    dserror("UseTranspose not impl.");
    return false;
  }

  /*!
  \brief not implemented
  */
  bool HasNormInf() const 
  {
    dserror("HasNormInf not impl.");
    return false;
  }

private:

  // don't want copy-ctor and = operator
  SIMPLER_Operator(SIMPLER_Operator& old);
  SIMPLER_Operator operator = (const SIMPLER_Operator& old);

  /*!
  \brief setup phase of preconditioner
  */
  void Setup(RCP<Epetra_CrsMatrix> A);

  /*!
  \brief do one sweep of simple preconditioning
  */
  void Simple(Epetra_Vector& vx, Epetra_Vector& px, Epetra_Vector& vb, Epetra_Vector& pb) const;

  FILE*                                               outfile_;// error file to write output to or NULL
  ParameterList                                       vlist_;  // list for velocity solver
  ParameterList                                       plist_;  // list for pressure solver
  MultiMapExtractor                                   mmex_;   // a  multimapetxractor to handle extracts
  RCP<BlockSparseMatrix<DefaultBlockMatrixStrategy> > A_;      // 2x2 block matrix
  RCP<Epetra_Vector>                                  vx_;     // velocity solution
  RCP<Epetra_Vector>                                  px_;     // pressure solution
  RCP<Epetra_Vector>                                  vb_;     // velocity rhs
  RCP<Epetra_Vector>                                  pb_;     // pressure rhs
  RCP<LINALG::Solver>                                 vsolver_;// velocity solver
  RCP<LINALG::Solver>                                 psolver_;// pressure solver

}; // class  SIMPLER_Operator : public virtual Epetra_Operator
}  // namespace LINALG


#endif // SIMPLER_OPERATOR_H
#endif  // #ifdef CCADISCRET
