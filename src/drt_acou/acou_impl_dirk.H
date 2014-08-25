/*!----------------------------------------------------------------------
\file acou_impl_dirk.H
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_DIRK_H
#define ACOU_IMPL_DIRK_H

#include "acou_impl.H"

namespace ACOU
{

class TimIntImplDIRK : public AcouImplicitTimeInt
{

public:
  /// constructor
  TimIntImplDIRK(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      );

  /// run a complete simulation - overrides base class, because we need to solve several stages each time step
  virtual void Integrate(Teuchos::RCP<Epetra_MultiVector> history = Teuchos::null, Teuchos::RCP<LINALG::MapExtractor> splitter = Teuchos::null);

  /// call elements to calculate system matrix
  void AssembleMatAndRHS(int stage);

  /// p-adaptivity
  virtual void UpdatePolyAndState();

  /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we can save some computing time)
  void UpdateInteriorVariablesAndAssemebleRHS(int stage);

  /// Routine to loop over all stages in one time step
  virtual void Solve();

  /// apply Dirichlet boudnary conditions to system
  void ApplyDirichletToSystem(int stage);

  /// return the name of the time integration scheme
  virtual std::string Name();

private:

  /// coefficients of given dirk scheme, filled in during construction
  double dirk_a_[6][6];
  double dirk_b_[6];
  int dirk_q_;

  // for dirk as class variable
  bool resonly_;

  /// vectors needed for interstage states, nomenclatur is in correspondence with Nguyen, 2011, "High-order implicit hybridizable
  /// discontinuous Galerkin methods for acoustics and elastodynamics".
  Teuchos::RCP<Epetra_Vector> t_;
};
} // namespace ACOU


#endif /* ACOU_IMPL_DIRK_H */
