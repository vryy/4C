/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
*/
/*----------------------------------------------------------------------*/
#ifndef AMGNXN_VCYCLE_H
#define AMGNXN_VCYCLE_H
#ifdef HAVE_MueLu

// Trilinos includes
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_UseDefaultTypes.hpp>

// Baci includes
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"
#include "solver_amgnxn_smoothers.H"


namespace LINALG
{
  namespace SOLVER
  {
    namespace AMGNXN
    {
      // class GenericSmoother;

      class Vcycle : public GenericSmoother
      {
       public:
        Vcycle(int NumLevels, int NumSweeps, int FirstLevel);

        void SetOperators(std::vector<Teuchos::RCP<BlockedMatrix>> Avec);
        void SetProjectors(std::vector<Teuchos::RCP<BlockedMatrix>> Pvec);
        void SetRestrictors(std::vector<Teuchos::RCP<BlockedMatrix>> Rvec);
        void SetPreSmoothers(std::vector<Teuchos::RCP<GenericSmoother>> SvecPre);
        void SetPosSmoothers(std::vector<Teuchos::RCP<GenericSmoother>> SvecPos);

        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const;

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

        void SetOperators(std::vector<Teuchos::RCP<SparseMatrix>> Avec);
        void SetProjectors(std::vector<Teuchos::RCP<SparseMatrix>> Pvec);
        void SetRestrictors(std::vector<Teuchos::RCP<SparseMatrix>> Rvec);
        void SetPreSmoothers(std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre);
        void SetPosSmoothers(std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos);

        virtual void Apply(
            const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const;

       private:
        void DoVcycle(const Epetra_MultiVector& X, Epetra_MultiVector& Y, int level,
            bool InitialGuessIsZero) const;

        int NumLevels_;
        int NumSweeps_;
        int FirstLevel_;

        std::vector<Teuchos::RCP<SparseMatrix>> Avec_;
        std::vector<Teuchos::RCP<SparseMatrix>> Pvec_;
        std::vector<Teuchos::RCP<SparseMatrix>> Rvec_;
        std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPre_;
        std::vector<Teuchos::RCP<SingleFieldSmoother>> SvecPos_;

        bool flag_set_up_A_;
        bool flag_set_up_P_;
        bool flag_set_up_R_;
        bool flag_set_up_Pre_;
        bool flag_set_up_Pos_;
      };

    }  // namespace AMGNXN
  }    // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLu
#endif  // AMGNXN_OPERATOR_H
