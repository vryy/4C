/*----------------------------------------------------------------------*/
/*!
\file mfsi_overlapalgorithm.H

\brief Solve monolithic overlapping block coupled FSI problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_OVERLAPALGORITHM_H
#define MFSI_OVERLAPALGORITHM_H

#include "mfsi_algorithm.H"

namespace MFSI
{

  /// monolithic FSI algorithm with overlapping interface equations
  class OverlapAlgorithm : public Algorithm
  {
  public:

    OverlapAlgorithm(Epetra_Comm& comm);

    /// Extract initial guess from fields
    void InitialGuess(Thyra::DefaultProductVector<double>& ig);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    void SetupRHS(Thyra::DefaultProductVector<double> &f) const;

    /// setup composed system matrix from field solvers
    void SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const;

    //@}

  private:

    /// setup of NOX convergance tests
    Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp);
  };
}

#endif
#endif
