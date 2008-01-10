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

  protected:

    /// setup of NOX convergance tests
    Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp);

    /// extract the three field vectors from a given composed vector
    virtual void ExtractFieldVectors(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) const;
  private:

    /// build thyra product vector from field vectors
    void SetupVector(Thyra::DefaultProductVector<double> &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale) const;

    ///
    void AddFluidInterface(double scale,
                           Teuchos::RCP<Epetra_CrsMatrix> fgg,
                           Teuchos::RCP<Epetra_CrsMatrix> s) const;

    ///
    Teuchos::RCP<Epetra_CrsMatrix>
    ConvertFigColmap(Teuchos::RCP<Epetra_CrsMatrix> fig,
                     const std::map<int,int>& convcolmap,
                     const Epetra_Map& domainmap) const;

    ///
    Teuchos::RCP<Epetra_CrsMatrix>
    ConvertFgiRowmap(Teuchos::RCP<Epetra_CrsMatrix> fgi,
                     double scale,
                     const Epetra_Map& structrowmap) const;

    //! @name thyra vector spaces

    /// structure field map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > smap_;

    /// structure field inner map
    //Teuchos::RCP<const Thyra::VectorSpaceBase<double> > simap_;

    /// structure field interface map
    //Teuchos::RCP<const Thyra::VectorSpaceBase<double> > sgmap_;

    /// fluid field inner map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > fimap_;

    /// ale field inner map
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > aimap_;

    //@}

    //! @name global system solvers

    /// solver factory for inner structure system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > structsolverfactory_;

    /// solver factory for interface system in preconditioner
    //Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > interfacesolverfactory_;

    /// solver factory for inner fluid system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > fluidsolverfactory_;

    /// solver factory for inner ale system in preconditioner
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > alesolverfactory_;

    //@}

    /// column fluid to structure dof map
    mutable std::map<int,int> permfluidstructcolmap_;

    mutable std::map<int,int> fluidstructcolmap_;

    mutable std::map<int,int> alestructcolmap_;
  };
}

#endif
#endif
