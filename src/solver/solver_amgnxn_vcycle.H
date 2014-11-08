/*!----------------------------------------------------------------------
\file solver_amgnxn_vcycle.H

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/
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
#include <MueLu_UseShortNames.hpp>

// Baci includes
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"


namespace LINALG
{
  namespace SOLVER
  {

    class AMGnxn_SmootherBase;

    class Richardson_Vcycle_Operator
    {
      public:

        Richardson_Vcycle_Operator
          (int NumLevels,int NumSweeps=1,double omega=1);

        void SetOperators  (std::vector< Teuchos::RCP<Epetra_Operator> > Avec);
        void SetProjectors (std::vector< Teuchos::RCP<Epetra_Operator> > Pvec);
        void SetRestrictors(std::vector< Teuchos::RCP<Epetra_Operator> > Rvec);
        void SetPreSmoothers (std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPre);
        void SetPosSmoothers (std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPos);

        void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,int start_level=0) const;

      private:

        void Richardson_Vcycle
          (const Epetra_MultiVector& X, Epetra_MultiVector& Y,int start_level) const;

        void Vcycle(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y, int level, bool InitialGuessIsZero) const;

        int NumLevels_;
        int NumSweeps_;
        double omega_;

        std::vector< Teuchos::RCP<Epetra_Operator> > Avec_;
        std::vector< Teuchos::RCP<Epetra_Operator> > Pvec_;
        std::vector< Teuchos::RCP<Epetra_Operator> > Rvec_;
        std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPre_;
        std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPos_;

        bool flag_set_up_A_;
        bool flag_set_up_P_;
        bool flag_set_up_R_;
        bool flag_set_up_Pre_;
        bool flag_set_up_Pos_;

    };


  } //namespace LINALG
} // namespace SOLVER

#endif // HAVE_MueLu
#endif // AMGNXN_OPERATOR_H
