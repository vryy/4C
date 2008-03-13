/*!----------------------------------------------------------------------
\file simpler_operator.H

\class LINALG::SIMPLER_Operator

\brief An approximate block factorization preconditioner based on the
       SIMPLE family of methods

This Operator implements the family of SIMPLE methods such as
SIMPLE, SIMPLER, SIMPLEC and CheapSIMPLE (a m.gee's variation of SIMPLE)

Literature:<br>

Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  
A taxonomy and comparison of parallel block multi-level              
preconditioners for the incomp. Navier-Stokes equations.             
Sandia technical report SAND2007-2761, 2007,                          
Also appeared in JCP                                                 

Pernice, M., Tocci, M.D.:                                            
A Multigrid Preconditioned Newton-Krylov method for the incomp.      
Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)  

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
#include "Epetra_Time.h"

#include "ml_MultiLevelPreconditioner.h"

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_systemmatrix.H"
#include "linalg_solver.H"
#include "linalg_utils.H"
#include "linalg_ana.H"

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

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{

/*!
\class SIMPLER_Operator

\brief A Semi-implicit Method for Pressure Linked Equations (SIMPLE) 
       preconditioner

This Operator implements the family of SIMPLE methods such as
SIMPLE, SIMPLER, SIMPLEC and CheapSIMPLE (a m.gee's variation of SIMPLE)

Literature:<br>

Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:  
A taxonomy and comparison of parallel block multi-level              
preconditioners for the incomp. Navier-Stokes equations.             
Sandia technical report SAND2007-2761, 2007,                          
Also appeared in JCP                                                 

Pernice, M., Tocci, M.D.:                                            
A Multigrid Preconditioned Newton-Krylov method for the incomp.      
Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)  

\author gee (gee@lnm.mw.tum.de)
*/
class  SIMPLER_Operator : public virtual Epetra_Operator
{
public:

  /*!
  \brief Standard Constructor
  */
  explicit SIMPLER_Operator(RCP<Epetra_Operator> A,
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
  void Setup(RCP<Epetra_Operator> A,
             const ParameterList& origvlist,
             const ParameterList& origplist);

  /*!
  \brief do one sweep of simple or simplec preconditioning
  */
  void Simple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
              LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const;

  /*!
  \brief do one sweep of simpler preconditioning
  */
  void Simpler(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
               LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const;


  /*!
  \brief do one sweep of simple or simplec preconditioning
  */
  void CheapSimple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px,
                   LINALG::ANA::Vector& vb, LINALG::ANA::Vector& pb) const;

  FILE*                                               outfile_; // error file to write output to or NULL
  ParameterList                                       vlist_;   // list for velocity solver
  ParameterList                                       plist_;   // list for pressure solver
  double                                              alpha_;   // pressure damping \in (0,1]

  MultiMapExtractor                                   mmex_;    // a  multimapetxractor to handle extracts
  RCP<BlockSparseMatrixBase>                          A_;       // 2x2 block matrix
  RCP<SparseMatrix>                                   diagAinv_;// inverse of main diagonal of A(0,0)
  RCP<SparseMatrix>                                   S_;       // Approximate Schur complement on the pressure space

  RCP<Epetra_Operator>                                Pv_;      // preconditioner for velocity subproblem
  RCP<Epetra_Operator>                                Pp_;      // preconditioner for pressure subproblem

  RCP<LINALG::ANA::Vector>                            vx_;      // velocity solution
  RCP<LINALG::ANA::Vector>                            px_;      // pressure solution
  RCP<LINALG::ANA::Vector>                            vb_;      // velocity rhs
  RCP<LINALG::ANA::Vector>                            pb_;      // pressure rhs
  RCP<LINALG::ANA::Vector>                            vwork1_;   // working vector velocity dimension
  RCP<LINALG::ANA::Vector>                            vwork2_;   // working vector velocity dimension
  RCP<LINALG::ANA::Vector>                            pwork1_;  // working vector pressure dimension
  RCP<LINALG::ANA::Vector>                            pwork2_;  // working vector pressure dimension

  RCP<LINALG::Solver>                                 vsolver_; // velocity solver
  RCP<LINALG::Solver>                                 psolver_; // pressure solver

}; // class  SIMPLER_Operator : public virtual Epetra_Operator
}  // namespace LINALG


#endif // SIMPLER_OPERATOR_H
#endif  // #ifdef CCADISCRET
