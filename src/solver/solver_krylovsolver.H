/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of Baci's interface to Krylov solvers

\level 0

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

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
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#endif                                // HAVE_MueLu

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
      //! Constructor
      KrylovSolver(const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile);

      //! Destructor
      virtual ~KrylovSolver();

      /*! \brief Setup the solver object
       *
       * @param matrix Matrix of the linear system
       * @param x Solution vector of the linear system
       * @param b Right-hand side vector of the linear system
       * @param refactor Boolean flag to enforce a refactorization of the matrix
       * @param reset Boolean flag to enforce a full reset of the solver object
       * @param projector Krylov projector
       */
      virtual void Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
          Teuchos::RCP<Epetra_MultiVector> b, const bool refactor, const bool reset,
          Teuchos::RCP<LINALG::KrylovProjector> projector) = 0;

      virtual int Solve() = 0;

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y);

      /*! \brief Get plain pointer to underlying preconditioner object of type PreconditionerType
       *
       * \warning Using this function outside the solver object is dangerous, since there's no
       * guarantee, that the preconditioner object still exists (but usually it should exist
       * as long it's needed/used within the solver object)
       */
      virtual PreconditionerType* Preconditioner() { return preconditioner_.get(); }

      int Ncall() { return ncall_; }

      Teuchos::ParameterList& Params() const { return params_; }

     protected:
      /*! \brief Check whether preconditioner will be reused
       *
       * The user can control reuse/recomputation of the preconditioner by setting appropriate input
       * arguments \c reuse and \c reset. In addition, contact mechanics problems perform some
       * additional checks since they require to rebuild the preconditioner when the active set has
       * changed.
       *
       * @param[in] reuse Parameter AZREUSE from parameter list
       * @param reset Force preconditioner to be rebuilt
       * @return Boolean flag to indicate whether preconditioner is reused (\c true) or has to be
       * re-computed (\c false)
       *
       * \sa CheckReuseStatusOfActiveSet()
       */
      bool AllowReusePreconditioner(const int reuse, const bool reset);

      /*! \brief Function for creating preconditioner object
       *
       * @param azlist Aztec parameter list
       * @param isCrsMatrix Boolean flag to indicate Epetra_CrsMatrix (true) or block matrix (false)
       * @param projector Krylov projector
       */
      void CreatePreconditioner(Teuchos::ParameterList& azlist, const bool isCrsMatrix,
          Teuchos::RCP<LINALG::KrylovProjector> projector);

      //! a communicator
      const Epetra_Comm& comm_;

      //! (internal) parameter lists
      Teuchos::ParameterList& params_;

      //! file to write output to (proc 0 only, can be NULL on input)
      FILE* outfile_;

      //! initial guess and solution
      Teuchos::RCP<Epetra_MultiVector> x_;

      //! right hand side vector
      Teuchos::RCP<Epetra_MultiVector> b_;

      //! system of equations
      Teuchos::RCP<Epetra_Operator> A_;

      //! counting how many times matrix was solved between resets
      int ncall_;

      //! preconditioner object
      Teuchos::RCP<PreconditionerType> preconditioner_;

     private:
      /*! \brief Check if active set has changed. If yes, enforce to rebuild the preconditioner.
       *
       * We only can reuse the preconditioner if the active set in contact problems has not changed.
       * Therefore, we compare the current map of active DOFs with #activeDofMap_, i.e. the one
       * stored from the previous application of the preconditioner.
       *
       * The comparison is done in multiple stages:
       * 1. Check if number of active DOFs has changed:
       *   - If yes, we need to rebuild the preconditioner.
       *   - If not, then check 2.
       * 2. Compare current and stored map of active DOFs using a PointSameAs() comparison.
       *   - If map has changed, we need to rebuild the preconditioner.
       *
       * \param[in/out] bAllowReuse Boolean flag to indicate reuse (true) or rebuild
       * (false) of the preconditioner
       * \param[in] linSysParams Parameter list with some linear system information
       */
      virtual void CheckReuseStatusOfActiveSet(
          bool& bAllowReuse, const Teuchos::ParameterList* linSysParams);

      /*! \brief Map of active DOFs in structural contact simulations.
       *
       * This is used to check whether preconditioner can be reused or not.
       * We cannot reuse the preconditioner if the active set has changed.
       *
       * \note Member variable only used for structural contact problems with the
       * "contact activeDofMap" parameter in the "Linear System properties"
       * parameter list set.
       *
       * \sa CheckReuseStatusOfActiveSet()
       */
      Teuchos::RCP<Epetra_Map> activeDofMap_;

#ifdef HAVE_MueLu

     protected:
      //! @name Permutation methods.
      //!@{

      /*! @brief extract map provided for being used with permutation factory. only rows in the
       * given map are subject to permutations. returns Teuchos::null if no additional information
       * is available/user-provided
       * */
      Teuchos::RCP<Epetra_Map> ExtractPermutationMap(
          const std::string solverSublist, const std::string mapName = "contact slaveDofMap");

      /*! @brief Member function to build permutation operators using matrix A as input for the
       * permutation factory. epSlaveDofMap defines a row map which restricts the permutation on
       * this map. */
      void BuildPermutationOperator(
          const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_Map>& epSlaveDofMap);

      /*! @brief Decide whether linear system is recommended for being permuted before solver is
       * called.
       *
       * decide whether it is recommended to solve a permuted linear system using the given
       * linear solver. returns true, if permutation is recommended.
       * Depending on the decision it feeds the preconditioner with the following additional
       * information: "Linear System Properties" -> "non diagonal-dominant row map"
       * (Teuchos::RCP<Map>) "Linear System Properties" -> "near-zero diagonal row map"
       * (Teuchos::RCP<Map>)*
       *
       * use non-permuted matrix A as input.
       * must be called after BuildPermutationOperator().
       *
       * */
      bool DecideAboutPermutation(const Teuchos::RCP<Epetra_CrsMatrix>& A);

      /*! @brief Do permutation of linear system. Routine sets permuted matrix and permuted rhs to
       * internal solver variables. */
      void PermuteLinearSystem(
          const Teuchos::RCP<Epetra_CrsMatrix>& A, const Teuchos::RCP<Epetra_MultiVector>& b);

      /*! @brief Retransform solution of permuted linear system to solution of original linear
       * system */
      void ReTransformSolution();

      /*! parameter xOp only needed for RowMap and DomainMap? (use original matrix A without
       * permutation)
       */
      void PermuteNullSpace(const Teuchos::RCP<Epetra_CrsMatrix>& A);

      /*! @brief Internal helper functions to access operators from internal data structure */
      Teuchos::RCP<const Epetra_CrsMatrix> GetOperator(
          const std::string name, const Teuchos::RCP<MueLu::FactoryBase>& fact);
      Teuchos::RCP<Epetra_CrsMatrix> GetOperatorNonConst(
          const std::string name, const Teuchos::RCP<MueLu::FactoryBase>& fact);

      /*! @brief Count number of zeros on diagonal of input matrix */
      int CountZerosOnDiagonalEpetra(const Teuchos::RCP<Epetra_CrsMatrix>& A);

      /*! @brief Count number of zeros on diagonal of input matrix */
      int CountZerosOnDiagonal(const Teuchos::RCP<const Xpetra::Matrix<double, int, int>>& xOp);

      /*! @brief Returns map with global row ids of rows with absolute value < tolerance on the
       * diagonal of matrix A */
      Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> FindZeroDiagonalEntries(
          const Teuchos::RCP<Epetra_CrsMatrix>& A, double tolerance = 1e-5);

      /*! @brief Returns map with global row ids of rows with absolute value < tolerance on the
       * diagonal of matrix A */
      Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> FindZeroDiagonalEntries(
          const Teuchos::RCP<Xpetra::Matrix<double, int, int>>& xA, double tolerance = 1e-5);

      /*! @brief Returns map with non-diagonal dominant rows of matrix A.
       *
       *  check the ratio nof the diagonal entry and the entry with
       *  maximal absolute value in the current row (diagEntry / maxEntry ) < diagDominance (=1.0
       * default).
       *
       *  - if the row is diagonal dominant the ratio is 1.0
       *  - if the row is not diagonal dominant, the ratio is < 1.0
       *  - if the row has a zero on the diagonal the ratio is zero
       *  - if the row is a zero row (-> singular matrix) we divide zero by zero
       * */
      Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> FindNonDiagonalDominantRows(
          const Teuchos::RCP<Epetra_CrsMatrix>& A, double diagDominanceRatio = 1.0);
      Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> FindNonDiagonalDominantRows(
          const Teuchos::RCP<Xpetra::Matrix<double, int, int>>& xA,
          double diagDominanceRatio = 1.0);

      /*! @brief if true, the linear solver is allowed to decide to solve a permuted linear system
       * (default: false)*/
      bool bAllowPermutation_;

      /*! @brief internal variable whether to solve permuted linear system or original linear system
       * (after algorithmic decision)*/
      bool bPermuteLinearSystem_;

      /*! @brief internal variable with the permutation strategy */
      std::string permutationStrategy_;

      /*! @brief maximum allowed diagonal dominance ratio. All row GIDs with diagonalEntry/maxEntry
       * < diagDominanceRatio are considered to be significantly non-diagonal dominant.*/
      double diagDominanceRatio_;

      /*! @brief internal flexible data container for storing permutation specific data (permutation
       * operators...) */
      Teuchos::RCP<MueLu::Level> data_;

      /*! @brief factory class which generates permutation data using input matrix A */
      Teuchos::RCP<MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> PermFact_;

      //!@}

#endif  // HAVE_MueLu

#ifdef WRITEOUTSTATISTICS
      double dtimeprecondsetup_;
#endif
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_KRYLOVSOLVER_H_ */
