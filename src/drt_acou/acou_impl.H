/*----------------------------------------------------------------------*/
/*!

\brief Acoustics Base Algorithm for implicit time integration

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/


#ifndef ACOU_IMPLTIMINT_H
#define ACOU_IMPLTIMINT_H

#include "acou_timeint.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | implicit time integration for acoustics               schoeder 01/14 |
 *----------------------------------------------------------------------*/
namespace ACOU
{
  // acoustics solver
  class AcouImplicitTimeInt : public AcouTimeInt
  {
   public:
    /// constructor
    AcouImplicitTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// virtual destructor
    virtual ~AcouImplicitTimeInt();

    /// set initial field to zero (needed for adjoint problem)
    virtual void SetInitialZeroField();

    /// set initial field by given function
    virtual void SetInitialField(int startfuncno);

    /// run a complete simulation
    virtual void Integrate();

    /// solve the system for trace and interior variables
    virtual void Solve();

    /// call elements to calculate system matrix
    virtual void AssembleMatAndRHS();

    /// apply Dirichlet boudnary conditions to system
    void ApplyDirichletToSystem();

    /// function to output user information
    void OutputToScreen();

    /// evaluate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    /// p-adaptivity
    virtual void UpdatePolyAndState();

    /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we
    /// can save some computing time)
    virtual void UpdateInteriorVariablesAndAssemebleRHS();

    /// Output nodal internal pressure field at new time step
    virtual void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec);

    /// Fill Touch Count Vector (needed for inverse analysis)
    // virtual void FillTouchCountVec(Teuchos::RCP<Epetra_Vector> touchcount);

    /// For result tests
    Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }

    /// return the name of the time integration scheme
    virtual std::string Name() = 0;

    // get values for speed of sound gradient
    virtual double GetSoSGradient(int rowid)
    {
      dserror("inverse analysis no longer supported with implicit time integration");
      return -1.0;
    }

    // get values for density gradient
    virtual double GetDensityGradient(int rowid)
    {
      dserror("inverse analysis no longer supported with implicit time integration");
      return -1.0;
    }

   protected:
    int sourcefuncno_;  /// function number of source term for wave equation

    /// cpu-time measures
    double dtele_;    /// element evaluation time
    double dtsolve_;  /// solver time

    bool acouopt_;       /// is optimization also with respect to acoustical parameters?
    int outputcount_;    /// output count in inverse analysis
    bool errormaps_;     /// flag indicating whether error maps are calculated
    double padapttol_;   /// p-adaptivity error tolerance
    bool calcerr_;       /// flag for error calculation
    bool allelesequal_;  /// flag true if all elements have same shape and same material (for faster
                         /// evaluation for uniform meshes)

    /// (standard) system matrix
    Teuchos::RCP<LINALG::SparseOperator> sysmat_;

    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector> residual_;

    /// element based error vector, only if desired by user
    Teuchos::RCP<Epetra_Vector> error_;
  };

  class AcouTimeIntAderTriTet : public AcouImplicitTimeInt
  {
   public:
    /// constructor
    AcouTimeIntAderTriTet(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// return the name of the time integration scheme
    virtual std::string Name() { return "AderTriTet"; }

    /// run a complete simulation
    virtual void Integrate();

    /// solve the system for trace and interior variables
    virtual void Solve();

    /// call elements to calculate system matrix
    virtual void AssembleMatAndRHS();

    /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we
    /// can save some computing time)
    virtual void UpdateInteriorVariablesAndAssemebleRHS();

    /// evaluate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

   private:
    Teuchos::RCP<Epetra_Vector> tempsrc_;
  };

}  // namespace ACOU

#endif /* ACOU_IMPLTIMINT_H */
