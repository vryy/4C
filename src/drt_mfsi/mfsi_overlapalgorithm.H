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
  /*!

    This is research. Highly experimental. Subject to change.

    Combine structure, fluid and ale field in one huge block matrix. Matching
    nodes. Overlapping equations at the interface.

    The structure equations contain both internal and interface part. Fluid
    and ale blocks are reduced to their respective internal parts. The fluid
    interface equations are added to the structure interface equations. There
    are no ale equations at the interface.

    Based on Newton's method within NOX. NOX computes the sum of all Newton
    increments. The evaluation method evalModelImpl() is always called with
    the sum x. However the meaning of this sum depends on the field blocks
    used.

    - The structure block calculates the displacement increments:
      \f$ \Delta \mathbf{d}^{n+1}_{i+1} = \mathbf{d}^{n+1}_{i+1} - \mathbf{d}^{n} \f$

    - The fluid block calculates the velocity (and pressure) increments:
      \f$ \Delta \mathbf{u}^{n+1}_{i+1} = \mathbf{u}^{n+1}_{i+1} - \mathbf{u}^{n} \f$

    - The ale block calculates the absolute mesh displacement:
      \f$ \Delta \mathbf{d}^{G,n+1} = \Delta \mathbf{d}^{n+1}_{i+1} + \mathbf{d}^{n} \f$

    We assume a very simple velocity -- displacement relation at the interface
    \f$ \mathbf{u}^{n+1}_{i+1} = \frac{1}{\Delta t} \Delta \mathbf{d}^{n+1}_{i+1} \f$

    \author u.kue
    \date 01/08
   */
  class OverlapAlgorithm : public Algorithm
  {
  public:

    /// build the algorithm
    OverlapAlgorithm(Epetra_Comm& comm);

    /// Extract initial guess from fields
    void InitialGuess(Thyra::DefaultProductVector<double>& ig);

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    void SetupRHS(Thyra::DefaultProductVector<double> &f) const;

    /// setup composed system matrix from field solvers
    void SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const;

    //@}

    /// Return a simple one-line description of this object.
    std::string description() const { return "Overlapping MFSI"; }

  protected:

    /// setup of NOX convergance tests
    Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Thyra::Group> grp);

    /// extract the three field vectors from a given composed vector
    /*!
      This is called with the NOX unknown x. Thus we have the absolute value
      here, no increments.
     */
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
                     const Epetra_Map& domainmap,
                     double scale) const;

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

    //! @name linear ALE matrix

    Teuchos::RCP<Epetra_CrsMatrix> aii_;
    Teuchos::RCP<Epetra_CrsMatrix> aig_;

    //@}

    //! @name time measurement

    Teuchos::RCP<Teuchos::Time> sysmattimer_;
    Teuchos::RCP<Teuchos::Time> igtimer_;
    Teuchos::RCP<Teuchos::Time> rhstimer_;
    Teuchos::RCP<Teuchos::Time> exctracttimer_;
    Teuchos::RCP<Teuchos::Time> fluidinterfacetimer_;
    Teuchos::RCP<Teuchos::Time> figcolmaptimer_;
    Teuchos::RCP<Teuchos::Time> fgirowmaptimer_;

    //@}

    /// column fluid to structure dof map
    mutable std::map<int,int> permfluidstructcolmap_;

    mutable std::map<int,int> fluidstructcolmap_;

    std::map<int,int> alestructcolmap_;

    mutable bool havefluidstructcolmap_;
  };
}

#endif
#endif
