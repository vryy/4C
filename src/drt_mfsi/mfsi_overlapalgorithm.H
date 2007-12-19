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

    /// build thyra product vector from field vectors
    void SetupVector(Thyra::DefaultProductVector<double> &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale) const;

    /// setup of NOX convergance tests
    Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp);

    /// structure field inner map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > simap_;

    /// structure field interface map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > sgmap_;

    /// fluid field inner map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > fimap_;

    /// ale field inner map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > aimap_;

  };
}

#endif
#endif
