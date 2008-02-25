#ifdef CCADISCRET

#ifndef FSI_MONOLITHICOVERLAP_H
#define FSI_MONOLITHICOVERLAP_H

#include "fsi_monolithic.H"

namespace FSI
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
    \date 02/08
   */
  class MonolithicOverlap : Monolithic
  {
  public:
    explicit MonolithicOverlap(Epetra_Comm& comm);

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams=0);

    //@}

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Teuchos::RCP<Epetra_Vector> f) const;

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

  protected:

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergance tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp);

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
     */
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax) const;


  private:

    /// build block vector from field vectors
    void SetupVector(Teuchos::RCP<Epetra_Vector> &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale) const;

    ///
    void AddFluidInterface(double scale,
                           const LINALG::SparseMatrix& fgg,
                           LINALG::SparseMatrix& s) const;

    ///
    Teuchos::RCP<LINALG::SparseMatrix>
    ConvertFigColmap(const LINALG::SparseMatrix& fig,
                     const std::map<int,int>& convcolmap,
                     const Epetra_Map& domainmap,
                     double scale) const;

    ///
    Teuchos::RCP<LINALG::SparseMatrix>
    ConvertFgiRowmap(const LINALG::SparseMatrix& fgi,
                     double scale,
                     const Epetra_Map& structrowmap) const;

    //! @name linear ALE matrix

    Teuchos::RCP<LINALG::SparseMatrix> aii_;
    Teuchos::RCP<LINALG::SparseMatrix> aig_;

    //@}

    /// block system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;

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

    mutable bool havefluidstructcolmap_;
    mutable bool havepermfluidstructcolmap_;
  };
}

#endif
#endif
