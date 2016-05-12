/*!----------------------------------------------------------------------
\file solver_krylovsolver.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/

#ifndef SOLVER_KRYLOVSOLVER_H_
#define SOLVER_KRYLOVSOLVER_H_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MueLu  // TODO remove this as soon as all use Q1/2013 MueLu or newer
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>

#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_FactoryBase_fwd.hpp>
#include <MueLu_Level_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>
#endif // HAVE_MueLu

#include "solver_solvertype.H"
#include "solver_preconditionertype.H"

namespace LINALG
{
class KrylovProjector;

namespace SOLVER
{

//! krylov subspace linear solvers (aztec) with right-side preconditioning
class KrylovSolver : public SolverType
{
  public:

  KrylovSolver( const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile );

  virtual ~KrylovSolver();

  virtual void Setup(
    Teuchos::RCP<Epetra_Operator>         matrix   ,
    Teuchos::RCP<Epetra_MultiVector>      x        ,
    Teuchos::RCP<Epetra_MultiVector>      b        ,
    bool                                  refactor ,
    bool                                  reset    ,
    Teuchos::RCP<LINALG::KrylovProjector> projector
    ) = 0;

  virtual int Solve() = 0;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

  /// get plain pointer to underlaying preconditioner object (of type PreconditionerType)
  /// note: using this function outside the solver object is dangerous, since there's no
  /// guarantee, that the preconditioner object still exists (but usually it should exist
  /// as long it's needed/used within the solver object)
  virtual PreconditionerType * Preconditioner() { return preconditioner_.get(); /*&*preconditioner_;*/ }

  int Ncall() { return ncall_; }

  Teuchos::ParameterList & Params() const { return params_; }

  protected:

  /// return true, if preconditoner can be reused
  /// reuse (input): AZREUSE from parameter list
  /// reset (input): force preconditioner to be rebuilt
  bool AllowReusePreconditioner(int reuse, bool reset);

  //! function for creating preconditioner object
  void CreatePreconditioner(
    Teuchos::ParameterList & azlist,
    bool isCrsMatrix,
    Teuchos::RCP<LINALG::KrylovProjector> projector
    );

  //! a communicator
  const Epetra_Comm&               comm_;

  //! (internal) parameter lists
  Teuchos::ParameterList&          params_;

  //! file to write output to (proc 0 only, can be NULL on input)
  FILE*                            outfile_;

  //! initial guess and solution
  Teuchos::RCP<Epetra_MultiVector> x_;

  //! right hand side vector
  Teuchos::RCP<Epetra_MultiVector> b_;

  //! system of equations
  Teuchos::RCP<Epetra_Operator>    A_;

  //! counting how many times matrix was solved between resets
  int                              ncall_;

  //! preconditioner object
  Teuchos::RCP<PreconditionerType> preconditioner_;

  private:
  //! number of active DOFs in structural contact simulations.
  //! This is used to check whether preconditioner can be reused or not.
  //! Member variable only used for structural contact problems with the
  //! "contact activeDofMap" parameter in the "Linear System properties"
  //! parameter list set.
  int                              nActiveDofs_;

#ifdef HAVE_MueLu

  protected:
  //! @name Permutation methods.
  //@{

  /*! @brief extract map provided for being used with permutation factory. only rows in the given map are subject to permutations.
   *  returns Teuchos::null if no additional information is available/user-provided
   * */
  Teuchos::RCP<Epetra_Map> ExtractPermutationMap(const std::string solverSublist, const std::string mapName = "contact slaveDofMap");

  /*! @brief Member function to build permutation operators using matrix A as input for the permutation factory. epSlaveDofMap defines a row map which restricts the permutation on this map. */
  void BuildPermutationOperator(const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_Map> & epSlaveDofMap);

  /*! @brief Decide whether linear system is recommended for being permuted before solver is called.
   *
   * decide whether it is recommended to solve a permuted linear system using the given
   * linear solver. returns true, if permutation is recommended.
   * Depending on the decision it feeds the preconditioner with the following additional information:
   *     "Linear System Properties" -> "non diagonal-dominant row map" (Teuchos::RCP<Map>)
   *     "Linear System Properties" -> "near-zero diagonal row map" (Teuchos::RCP<Map>)*
   *
   * use non-permuted matrix A as input.
   * must be called after BuildPermutationOperator().
   *
   * */
  bool DecideAboutPermutation(const Teuchos::RCP<Epetra_CrsMatrix> & A);

  /*! @brief Do permutation of linear system. Routine sets permuted matrix and permuted rhs to internal solver variables. */
  void PermuteLinearSystem(const Teuchos::RCP<Epetra_CrsMatrix>& A,const Teuchos::RCP<Epetra_MultiVector>& b);

  /*! @brief Retransform solution of permuted linear system to solution of original linear system */
  void ReTransformSolution();

  // parameter xOp only needed for RowMap and DomainMap? (use original matrix A without permutation)
  void PermuteNullSpace(const Teuchos::RCP<Epetra_CrsMatrix>& A);

  /*! @brief Internal helper functions to access operators from internal data structure */
  Teuchos::RCP<const Epetra_CrsMatrix> GetOperator(const std::string name, const Teuchos::RCP<FactoryBase> & fact);
  Teuchos::RCP<Epetra_CrsMatrix>       GetOperatorNonConst(const std::string name, const Teuchos::RCP<FactoryBase> & fact);

  /*! @brief Count number of zeros on diagonal of input matrix */
  int CountZerosOnDiagonalEpetra(const Teuchos::RCP<Epetra_CrsMatrix>& A);
  /*! @brief Count number of zeros on diagonal of input matrix */
  int CountZerosOnDiagonal(const Teuchos::RCP<const Xpetra::Matrix<double, int, int> >& xOp);

  /*! @brief Returns map with global row ids of rows with absolute value < tolerance on the diagonal of matrix A */
  Teuchos::RCP<Map> FindZeroDiagonalEntries(const Teuchos::RCP<Epetra_CrsMatrix>& A, double tolerance = 1e-5);
  /*! @brief Returns map with global row ids of rows with absolute value < tolerance on the diagonal of matrix A */
  Teuchos::RCP<Map> FindZeroDiagonalEntries(const Teuchos::RCP<Xpetra::Matrix<double, int, int> >& xA, double tolerance = 1e-5);

  /*! @brief Returns map with non-diagonal dominant rows of matrix A.
   *
   *  check the ratio nof the diagonal entry and the entry with
   *  maximal absolute value in the current row (diagEntry / maxEntry ) < diagDominance (=1.0 default).
   *
   *  - if the row is diagonal dominant the ratio is 1.0
   *  - if the row is not diagonal dominant, the ratio is < 1.0
   *  - if the row has a zero on the diagonal the ratio is zero
   *  - if the row is a zero row (-> singular matrix) we divide zero by zero
   * */
  Teuchos::RCP<Map> FindNonDiagonalDominantRows(const Teuchos::RCP<Epetra_CrsMatrix>& A, double diagDominanceRatio = 1.0);
  Teuchos::RCP<Map> FindNonDiagonalDominantRows(const Teuchos::RCP<Xpetra::Matrix<double, int, int> >& xA, double diagDominanceRatio = 1.0);

  /*! @brief if true, the linear solver is allowed to decide to solve a permuted linear system (default: false)*/
  bool bAllowPermutation_;

  /*! @brief internal variable whether to solve permuted linear system or original linear system (after algorithmic decision)*/
  bool bPermuteLinearSystem_;

  /*! @brief internal variable with the permutation strategy */
  std::string permutationStrategy_;

  /*! @brief maximum allowed diagonal dominance ratio. All row GIDs with diagonalEntry/maxEntry < diagDominanceRatio are considered to be significantly non-diagonal dominant.*/
  double diagDominanceRatio_;

  /*! @brief internal flexible data container for storing permutation specific data (permutation operators...) */
  Teuchos::RCP<MueLu::Level> data_; // data container

  /*! @brief factory class which generates permutation data using input matrix A */
  Teuchos::RCP<PermutationFactory> PermFact_; // permutation factory
  //@}
#endif // HAVE_MueLu

#ifdef WRITEOUTSTATISTICS
  double dtimeprecondsetup_;
#endif
};
}
}

#endif /* SOLVER_KRYLOVSOLVER_H_ */
