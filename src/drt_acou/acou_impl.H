/*----------------------------------------------------------------------*/
/*!
\file acou_impl.H

\brief Acoustics Base Algorithm for implicit time integration

       Derived classes are implicit Euler, trapezoidal rule, BDF2,
       BDF3, BDF4 and several diagonally implicit RK-schemes

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ACOU_IMPLTIMINT_H
#define ACOU_IMPLTIMINT_H

#include "acou_timeint.H"

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
  AcouImplicitTimeInt(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
  );

  /// virtual destructor
  virtual ~AcouImplicitTimeInt();

  /// set initial field to zero (needed for adjoint problem)
  virtual void SetInitialZeroField();

  /// set initial field by given function
  virtual void SetInitialField(int startfuncno);

  /// set initial field by given light distribution
  virtual void SetInitialPhotoAcousticField(Teuchos::RCP<Epetra_Vector> light, Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform);

  /// run a complete simulation
  virtual void Integrate(Teuchos::RCP<Epetra_MultiVector> history = Teuchos::null, Teuchos::RCP<LINALG::MapExtractor> splitter = Teuchos::null);

  /// solve the system for trace and interior variables
  virtual void Solve();

  /// call elements to calculate system matrix
  virtual void AssembleMatAndRHS();

  /// apply Dirichlet boudnary conditions to system
  void ApplyDirichletToSystem();

  /// update: current solution becomes old solution of next timestep
  virtual void TimeUpdate();

  /// output
  void Output(Teuchos::RCP<Epetra_MultiVector> history = Teuchos::null, Teuchos::RCP<LINALG::MapExtractor> splitter = Teuchos::null);

  /// function to output user information
  void OutputToScreen();

  /// evaluate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /// p-adaptivity
  virtual void UpdatePolyAndState();

  /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we can save some computing time)
  virtual void UpdateInteriorVariablesAndAssemebleRHS();

  /// Output nodal internal pressure field at new time step
  void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec);

  /// Output nodal internal psi field at new time step
  void NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec);

  /// Fill Touch Count Vector (needed for inverse analysis)
  void FillTouchCountVec(Teuchos::RCP<Epetra_Vector> touchcount);

  /// For result tests
  Teuchos::RCP<Epetra_Vector> Velnp()   { return velnp_; }

  /// init a monitor file
  void InitMonitorFile();

  /// fill monitor file
  void FillMonitorFile(Teuchos::RCP<Epetra_Vector> ip);

  /// return the name of the time integration scheme
  virtual std::string Name() = 0;

protected:

  int     sourcefuncno_;                /// function number of source term for wave equation

  /// cpu-time measures
  double  dtele_;                       /// element evaluation time
  double  dtsolve_;                     /// solver time

  bool writemonitor_;                   /// flag indicating if a monitor file is written
  bool writestress_;                   /// flag indicating if a stresses or -if no- the displacment gradient is written to the output
  bool errormaps_;                      /// flag indicating whether error maps are calculated
  double padapttol_;                    /// p-adaptivity error tolerance
  bool calcerr_;                        /// flag for error calculation
  bool allelesequal_;                   /// flag true if all elements have same shape and same material (for faster evaluation for uniform meshes)

  /// (standard) system matrix
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// (standard) residual vector (rhs for the incremental form),
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// needed in case of adjoint problem solve
  Teuchos::RCP<Epetra_MultiVector>    adjoint_rhs_;

  /// element based error vector, only if desired by user
  Teuchos::RCP<Epetra_Vector>    error_;
};

} // namespace ACOU

#endif /* ACOU_IMPLTIMINT_H */
