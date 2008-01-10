#ifndef VM3_SOLVER_H
#define VM3_SOLVER_H

#ifdef CCADISCRET

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_CrsMatrix.h>
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

class  VM3_Solver
{
public:

  // @{ \name Constructors and destructors


  /*!
  \brief ctor

  */
  explicit VM3_Solver( RefCountPtr<Epetra_CrsMatrix>& Sep,
                       RefCountPtr<Epetra_CrsMatrix>& Aplus,
                       RefCountPtr<Epetra_CrsMatrix>& A,
                       RefCountPtr<Epetra_Vector>& sugrvisc,
                       RefCountPtr<Epetra_Vector>& rplus,
                       RefCountPtr<Epetra_Vector>& r,
                       RefCountPtr<Epetra_Vector>& sol,
                       const RefCountPtr<Epetra_Vector> dbctoggle,
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
  \brief get Comm of this class

  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

  */
  const Epetra_Comm & Comm() const {return(A_->Comm());}


  //@}
  // @{ \name Apply methods


  /*!
  \brief routine for generating the scale-separation matrix S

   (will only be called in the first timestep) 

  */
  bool Compute();

  /*!
  \brief scale-separating routine

   precomputation of unscaled S^T*M*S
   (called in the first timestep)

   \param Sep  (In)       : scale-separation matrix
   \param Aplus (In/Out)  : (fine-scale) subgrid-viscosity matrix

  */
   void Separate(RefCountPtr<Epetra_CrsMatrix>& Sep,
                 RefCountPtr<Epetra_CrsMatrix>& Aplus);

  /*!
  \brief scaling routine

   scale precomput. matrix product by subgrid-viscosity-scaling vector
   (called in every timestep)

   \param Aplus (In)      : (fine-scale) subgrid-viscosity matrix
   \param A (In/Out)      : complete system matrix
   \param sugrvisc (In)   : subgrid-viscosity-scaling vector
   \param rplus (In)      : (fine-scale) subgrid-viscosity rhs vector
   \param r (In/Out)      : complete rhs vector
   \param sol (In)        : current solution vector
   \param increm (In)     : flag for incremental formulation

  */
   void Scale(RefCountPtr<Epetra_CrsMatrix>& Aplus,
              RefCountPtr<Epetra_CrsMatrix>& A,
              RefCountPtr<Epetra_Vector>& sugrvisc,
              RefCountPtr<Epetra_Vector>& rplus,
              RefCountPtr<Epetra_Vector>& r,
              RefCountPtr<Epetra_Vector>& sol,
              bool increm);

  /*!
  \brief Add a Epetra_CrsMatrix to another, where the first one is vector-scaled: B = B + A*V

  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Range, Row and Domain maps of A and B have to match.

  Note that this is a true parallel add!

  \param A          (in)     : Matrix to add to B (must have Filled()==true)
  \param V          (in)     : scaling vector for A
  \param B          (in/out) : Matrix to be added to (must have Filled()==false)
  */
  void AddVecScal(const Epetra_CrsMatrix& A, const Epetra_Vector& V,
                        Epetra_CrsMatrix& B);

  /*!
  \brief Add a Epetra_CrsMatrix to another, where the first one is vector-scaled: B = B + A*V and calculate right hand side by multiplying vector-scaled (A*V) with solution vector S: R = (A*V)*S

  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Range, Row and Domain maps of A and B have to match.

  Note that this is a true parallel add!

  \param A          (in)     : Matrix to add to B (must have Filled()==true)
  \param S          (in)     : solution vector
  \param V          (in)     : scaling vector for A
  \param B          (in/out) : Matrix to be added to (must have Filled()==false)
  \param R          (out)    : right-hand-side vector
  */
  void AddVecScalMult(const Epetra_CrsMatrix& A, const Epetra_Vector& S,
                      const Epetra_Vector& V, Epetra_CrsMatrix& B, Epetra_Vector& R);

  /*!
  \brief Calculate right hand side by multiplying vector-scaled (A*V) with solution vector S: R = (A*V)*S

  \param A          (in)     : Matrix A
  \param S          (in)     : solution vector
  \param V          (in)     : scaling vector for A
  \param R          (out)    : right-hand-side vector
  */
  void VecScalMult(const Epetra_CrsMatrix& A, const Epetra_Vector& S,
                   const Epetra_Vector& V, Epetra_Vector& R);


private:

  // don't want copy-ctor and = operator
  VM3_Solver(VM3_Solver& old);
  VM3_Solver operator = (const VM3_Solver& old);

  bool compute_;  // flag whether scale-separation matrix should be computed

  ParameterList& mlparams_;    // ML parameter list

  RCP<Epetra_CrsMatrix>& Sep_;       // scale-separation matrix
  RCP<Epetra_CrsMatrix>& Aplus_;     // (fine-scale) subgrid-viscosity matrix
  RCP<Epetra_CrsMatrix>& A_;         // complete system matrix
  RCP<Epetra_Vector>& sugrvisc_;     // subgrid-viscosity-scaling vector
  RCP<Epetra_Vector>& rplus_;        // (fine-scale) subgrid-viscosity rhs vector
  RCP<Epetra_Vector>& r_;            // complete rhs vector
  RCP<Epetra_Vector>& sol_;          // current solution vector

  const RCP<Epetra_Vector>& dbctoggle_;  // DBC toggle vector

}; // VM3_Solver


#endif

#endif
