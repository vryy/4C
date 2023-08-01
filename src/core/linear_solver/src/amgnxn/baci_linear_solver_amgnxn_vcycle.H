/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_AMGNXN_VCYCLE_H
#define BACI_LINEAR_SOLVER_AMGNXN_VCYCLE_H

// Trilinos includes
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

// Baci includes
#include "baci_linalg_blocksparsematrix.H"
#include "baci_linear_solver_amgnxn_smoothers.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_linear_solver_preconditioner_type.H"


namespace CORE::LINEAR_SOLVER::AMGNXN
{
  class Vcycle : public GenericSmoother
  {
   public:
    Vcycle(int NumLevels, int NumSweeps, int FirstLevel);

    void SetOperators(std::vector<Teuchos::RCP<BlockedMatrix>> Avec);
    void SetProjectors(std::vector<Teuchos::RCP<BlockedMatrix>> Pvec);
    void SetRestrictors(std::vector<Teuchos::RCP<BlockedMatrix>> Rvec);
    void SetPreSmoothers(std::vector<Teuchos::RCP<GenericSmoother>> SvecPre);
    void SetPosSmoothers(std::vector<Teuchos::RCP<GenericSmoother>> SvecPos);

    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override;

   private:
    void DoVcycle(
        const BlockedVector& X, BlockedVector& Y, int level, bool InitialGuessIsZero) const;

    int NumLevels_;
    int NumSweeps_;
    int FirstLevel_;

    std::vector<Teuchos::RCP<BlockedMatrix>> Avec_;
    std::vector<Teuchos::RCP<BlockedMatrix>> Pvec_;
    std::vector<Teuchos::RCP<BlockedMatrix>> Rvec_;
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPre_;
    std::vector<Teuchos::RCP<GenericSmoother>> SvecPos_;

    bool flag_set_up_A_;
    bool flag_set_up_P_;
    bool flag_set_up_R_;
    bool flag_set_up_Pre_;
    bool flag_set_up_Pos_;
  };

  // This could be done better with templates
  class VcycleSingle : public SingleFieldSmoother
  {
   public:
    VcycleSingle(int NumLevels, int NumSweeps, int FirstLevel);

    void SetOperators(std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Avec);
    void SetProjectors(std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Pvec);
    void SetRestrictors(std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Rvec);
    void SetPreSmoothers(std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre);
    void SetPosSmoothers(std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos);

    void Apply(
        const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const override;

   private:
    void DoVcycle(const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level,
        bool InitialGuessIsZero) const;

    int NumLevels_;
    int NumSweeps_;
    int FirstLevel_;

    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Avec_;
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Pvec_;
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> Rvec_;
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre_;
    std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos_;

    bool flag_set_up_A_;
    bool flag_set_up_P_;
    bool flag_set_up_R_;
    bool flag_set_up_Pre_;
    bool flag_set_up_Pos_;
  };
}  // namespace CORE::LINEAR_SOLVER::AMGNXN

#endif  // SOLVER_AMGNXN_VCYCLE_H
