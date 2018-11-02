/*----------------------------------------------------------------------*/
/*!
\file inpar_solver.H

\brief Input parameters for linear solvers

\level 1

\maintainer Martin Kronbichler
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_SOLVER_H
#define INPAR_SOLVER_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace SOLVER
  {
    /// Available solvers (a few of these are likely out-dated)
    ///
    /// \author m.gee
    /// \date 7/01,11/08
    enum SolverType
    {
      aztec_msr,           ///< solver package aztec 2.1, matrix in dmsr format
      lapack_sym,          ///< symmetric lapack LU decomposition, matrix is dense
      lapack_nonsym,       ///< unsymmetric lapack LU decomposition, matrix is dense
      umfpack,             ///< solver package umfpack, matrix in compressed column format
      amesos_klu_sym,      ///< Trilinos' solver interface Amesos using KLU
      amesos_klu_nonsym,   ///< Trilinos' solver interface Amesos using KLU
      superlu,             ///< Trilinos' solver interface Amesos using SuperLU_Dist
      belos,               ///< Belos solver
      stratimikos_amesos,  ///< Stratimikos solver framework (Amesos KLU)
      stratimikos_aztec,   ///< Stratimikos solver framework (Aztec)
      stratimikos_belos,   ///< Stratimikos solver framework (Belos)
      undefined            ///< undefined solver
    };

    /// Different solvers within the Aztec2.1 library
    ///
    /// \author m.gee
    /// \date 9/01,11/08
    enum AzSolverType
    {
      azsolv_CG,        ///< cg-solver
      azsolv_GMRES,     ///< gmres-solver
      azsolv_GMRESR,    ///< restarted gmres-solver
      azsolv_CGS,       ///< cg-squared-solver (can handle unsymmetric problems)
      azsolv_BiCGSTAB,  ///< bicgstab-solver
      azsolv_LU,        ///< LU-solver (I think this doesn't work for some reasons)
      azsolv_TFQMR,     ///< quasi-minimum residual-solver (never used it)
      belos_FGMRES      ///< flexible GMRES in Belos
    };

    /// Different preconditioners within the Aztec package
    ///
    /// \author m.gee
    /// \date 9/01,11/08
    enum AzPrecType
    {
      azprec_none,                 ///< no preconditioning
      azprec_ILUT,                 ///< incomplete LU fact. with numerical drop tolerance
      azprec_ILU,                  ///< incomplete LU fact. with fill in levels
      azprec_Jacobi,               ///< Jacobi preconditioning
      azprec_Neumann,              ///< neumann polynomials
      azprec_Least_Squares,        ///< least squares something
      azprec_SymmGaussSeidel,      ///< symmetric n-step gauss-seidel preconditioner
      azprec_GaussSeidel,          ///< symmetric n-step gauss-seidel preconditioner
      azprec_DownwindGaussSeidel,  ///< symmetric n-step gauss-seidel preconditioner
      azprec_LU,                   ///< ?
      azprec_RILU,                 ///< relaxed incomplete LU
      azprec_ICC,                  ///< incomplete cholesky
      azprec_MLfluid,              ///< ML for Fluids
      azprec_MLfluid2,             ///< energy optimal unsymmetric ML for Fluids
      azprec_MLAPI,                ///< experimental mlapi preconditioner
      azprec_ML,                   ///< standard ML for structures
      azprec_BGS2x2,               ///< Block Gauss-Seidel for 2x2 system
      azprec_TekoSIMPLE,           ///< SIMPLE preconditioner for 2x2 block system (Teko)
      azprec_BGSnxn,               ///< Block Gauss-Seidel for nxn system (Teko)
      azprec_CheapSIMPLE,   ///< CheapSIMPLE (BACI implementation for contact/meshtying/constraints)
      azprec_MueLuAMG_sym,  ///< MueLu AMG preconditioner (Smoothed aggregation for symmetric
                            ///< problems)
      azprec_MueLuAMG_nonsym,      ///< MueLu AMG preconditioner (PG-AMG for nonsymmetric problems ~
                                   ///< MLFLUID2)
      azprec_MueLuAMG_contact,     ///< MueLu AMG preconditioner (special handling for condensed
                                   ///< contact formulation)
      azprec_MueLuAMG_contact2,    ///< MueLu AMG preconditioner (special enhanced handling for
                                   ///< condensed contact formulation)
      azprec_MueLuAMG_contactSP,   ///< MueLu AMG preconditioner for saddlepoint formulation
      azprec_MueLuAMG_contactPen,  ///< MueLu AMG preconditioner for penalty formulation
      azprec_AMGnxn,               ///< AMG preconditioner for a nxn block matrix
      azprec_Chebyshev             ///< Chebyshev polynomials for preconditioning
    };

    /// permutation strategies for linear solvers
    ///
    /// \author tw
    /// \date 02/13
    enum PermutationStrategy
    {
      Permutation_none,       /// < no permutation, use original linear system
      Permutation_algebraic,  /// < use algebraic permutation strategy
      Permutation_local       /// < use local permutations
    };

    /// reuse strategies for MueLu Contact preconditioner
    ///
    /// up to now these parameters are only used by the MueLu_ContactPreconditioner2
    ///
    /// \author tw
    /// \date 10/13
    enum ReuseStrategy
    {
      Reuse_nothing,  /// < do not reuse multigrid preconditioner (even if AZREUSE is > 1)
      Reuse_Ptent,  /// < try to reuse non-smoothed tentative prolongator (prerequisite AZREUSE > 1
                    /// appropriately set)
      Reuse_full    /// < try to reuse full multigrid hierarchy (this is exactly the behaviour as we
                    /// expect it from AZREUSE)
    };

    /// set the parameters for the linear solver
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace SOLVER

}  // namespace INPAR


/*----------------------------------------------------------------------*/

#endif  // #ifndef INPAR_SOLVER_H
