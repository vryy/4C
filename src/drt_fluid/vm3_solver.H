#ifndef VM3_SOLVER_H
#define VM3_SOLVER_H

#ifdef CCADISCRET

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_genP.h"
#include "ml_operator.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"


using namespace std;
using namespace Teuchos;
using namespace MLAPI;


/*!
\class VM3_Solver

\brief


\author vg

*/

class  VM3_Solver : public virtual Epetra_Operator
{
public:

  // @{ \name Constructors and destructors


  /*!
  \brief ctor

  */
  explicit VM3_Solver( RefCountPtr<Epetra_CrsMatrix> Aplus,
                       RefCountPtr<Epetra_CrsMatrix> A,
                       ParameterList& mlparams,
                       bool compute);

  /*!
  \brief Destroys an instance of this class

  Destructor

  */
  virtual ~VM3_Solver() {}

  //@}
  // @{ \name Query methods

  /*!
  \brief Returns the label of this class.
  */
  const char* Label() const {return(&(label_[0]));}

  /*!
  \brief get Comm of this class

  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

  */
  const Epetra_Comm & Comm() const {return(A_->Comm());}

  /*!
  \brief Get fine level OperatorDomainMap

  Derived from Epetra_Operator, get fine level OperatorDomainMap

  */
  const Epetra_Map & OperatorDomainMap() const {return A_->OperatorDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap

  Derived from Epetra_Operator, get fine level OperatorRangeMap

  */
  const Epetra_Map & OperatorRangeMap() const {return A_->OperatorRangeMap();}

  //@}
  // @{ \name Apply methods


  /*!
  \brief Compute the preconditioner

   Compute the vm3 solver
  */
  bool Compute();

  /*!
  \brief ApplyInverse the preconditioner

   ApplyInverse the preconditioner. Method is derived from Epetra_Operator.


   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  /*!
  \brief Apply the VM3 solver

   Perform one sweep of the VM3 solver


   \param X  (In)  : Epetra_MultiVector matching the fine level map of this
                     solver
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int Solve(const Epetra_MultiVector& X, Epetra_MultiVector& Y);


  // @{ \name Not implemented methods derived from Epetra_Operator

  /*!
  \brief not implemented
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  { cout << "VM3_Solver::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";throw -1;return(-1);}

  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose) { cout << "VM3_Solver::SetUseTranspose not impl.\n";return -1;}

  /*!
  \brief not implemented
  */
  double NormInf() const {cout << "VM3_Solver::NormInf not impl.\n";
                          return(-1.0);}


  /*!
  \brief not implemented
  */
  bool UseTranspose() const {cout << "VM3_Solver::UseTranspose not impl.\n"; return false;}

  /*!
  \brief not implemented
  */
  bool HasNormInf() const {cout << "VM3_Solver::HasNormInf not impl.\n"; return false;}

  //@}

private:

  // don't want copy-ctor and = operator
  VM3_Solver(VM3_Solver& old);
  VM3_Solver operator = (const VM3_Solver& old);

  // the v-cycle
  int MultiLevelVCycle(MultiVector& b_f, MultiVector& x_f) const;

  // restriction operation (from level 0 to level)
  int Restrict(const MultiVector& x_f, MultiVector& x_r, int level) const;

  // prolongation operation (from level to level 0)
  int Prolong(const MultiVector& x_l, MultiVector& x_0, int level) const;

  // some helpers
  inline const Operator& P(const int i) const { return (mlapiPmod_[i]); }
  inline const Operator& R(const int i) const { return (mlapiRmod_[i]); }
//  inline const Operator& RP(const int i) const { return (mlapiRP_[i]); }
//  inline const Operator& PR(const int i) const { return (mlapiPR_[i]); }
  inline const Operator& RA(const int i) const { return (mlapiRA_[i]); }
  inline const Operator& A(const int i) const { return (mlapiA_[i]); }
  inline const Operator& Aplus(const int i) const { return (mlapiAplus_[i]); }
  inline const InverseOperator& S(const int i) const { return (mlapiS_[i]); }

  const char*                    label_;       // the label of this class
  bool                           iscomputed_;  // flag indicating whether preconditioner is ready for use
  ParameterList&                 mlparams_;    // ML parameter list
  int                            maxlevels_;   // number of levels

  RefCountPtr<Epetra_CrsMatrix>  Aplus_;  // input matrix including fine-scale term
  RefCountPtr<Epetra_CrsMatrix>  A_;      // input matrix without fine-scale term

  vector<Operator>               mlapiRmod_;   // restriction
  vector<Operator>               mlapiPmod_;   // prolongation
  vector<Operator>               mlapiRP_;     // restriction * prolongatoin
  vector<Operator>               mlapiPR_;     // prolongation * restriction
  vector<Operator>               mlapiRA_;     // restriction * A
  vector<Operator>               mlapiAplus_;  // operator Aplus on each grid (actually needed only on finest grid)
  vector<Operator>               mlapiA_;      // operator A on each grid
  vector<InverseOperator>        mlapiS_;      // smoother on each grid

}; // VM3_Solver


#endif

#endif
