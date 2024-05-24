/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_VCYCLE_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_VCYCLE_HPP

// Trilinos includes
#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linear_solver_amgnxn_smoothers.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


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
    void do_vcycle(
        const BlockedVector& X, BlockedVector& Y, int level, bool InitialGuessIsZero) const;

    int num_levels_;
    int num_sweeps_;
    int first_level_;

    std::vector<Teuchos::RCP<BlockedMatrix>> avec_;
    std::vector<Teuchos::RCP<BlockedMatrix>> pvec_;
    std::vector<Teuchos::RCP<BlockedMatrix>> rvec_;
    std::vector<Teuchos::RCP<GenericSmoother>> svec_pre_;
    std::vector<Teuchos::RCP<GenericSmoother>> svec_pos_;

    bool flag_set_up_a_;
    bool flag_set_up_p_;
    bool flag_set_up_r_;
    bool flag_set_up_pre_;
    bool flag_set_up_pos_;
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
    void do_vcycle(const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level,
        bool InitialGuessIsZero) const;

    int num_levels_;
    int num_sweeps_;
    int first_level_;

    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> avec_;
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> pvec_;
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> rvec_;
    std::vector<Teuchos::RCP<SingleFieldSmoother>> svec_pre_;
    std::vector<Teuchos::RCP<SingleFieldSmoother>> svec_pos_;

    bool flag_set_up_a_;
    bool flag_set_up_p_;
    bool flag_set_up_r_;
    bool flag_set_up_pre_;
    bool flag_set_up_pos_;
  };
}  // namespace CORE::LINEAR_SOLVER::AMGNXN

FOUR_C_NAMESPACE_CLOSE

#endif
