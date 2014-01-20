/*!----------------------------------------------------------------------
\file acou_impl_bdf.H
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_BDF_H
#define ACOU_IMPL_BDF_H

#include "acou_impl.H"

namespace ACOU
{

class TimIntImplBDF : public AcouImplicitTimeInt
{

public:
  /// constructor
  TimIntImplBDF(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      );

  /// read restart information for given time step
  virtual void ReadRestart(int step);

  /// set initial field to zero (needed for adjoint problem)
  virtual void SetInitialZeroField();

  /// set initial field by given function
  virtual void SetInitialField(int startfuncno, double pulse);

  /// set initial field by given light distribution
  virtual void SetInitialPhotoAcousticField(double pulse, Teuchos::RCP<Epetra_Vector> light, Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform);

  /// call elements to calculate system matrix
  virtual void AssembleMatAndRHS();

  /// update: current solution becomes old solution of next timestep
  virtual void TimeUpdate();

  /// write restart
  virtual void WriteRestart();

  /// Call elements to update interior variabels and to calculate the RHS vector (in this way, we can save some computing time)
  virtual void UpdateInteriorVariablesAndAssemebleRHS();

  /// return the name of the time integration scheme
  virtual std::string Name();

private:
  /// order of the BDF scheme
  int order_;

  // in contrast to the base class, BDF schemes require additional history variables
  Teuchos::RCP<Epetra_Vector>    velnmm_;
  Teuchos::RCP<Epetra_Vector>    velnmmm_;
  Teuchos::RCP<Epetra_Vector>    intvelnmm_;
  Teuchos::RCP<Epetra_Vector>    intvelnmmm_;
};

} // namespace ACOU

#endif /* ACOU_IMPL_BDF_H */
