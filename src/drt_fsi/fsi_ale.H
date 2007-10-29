
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_ALE_H
#define FSI_ALE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../io/io_drt.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>

#include <set>


namespace FSI {


/// General ALE interface
class Ale
{
public:
  virtual ~Ale() {}

  /// announce a new time step
  //virtual void NewStep(int step, double time) {}

  virtual void PrepareTimeStep() = 0;
  virtual void Solve() = 0;
  virtual void Update() = 0;
  virtual void Output() = 0;

  virtual void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im) = 0;

  virtual void ApplyInterfaceDisplacements(Teuchos::RefCountPtr<Epetra_Vector> idisp) = 0;

  virtual const DRT::Discretization& Discretization() const = 0;

  virtual Teuchos::RefCountPtr<Epetra_Vector> ExtractDisplacement() = 0;
};


/// linear ALE algorithm
class AleLinear : public Ale,
                  public NOX::Epetra::Interface::Required,
                  public NOX::Epetra::Interface::Jacobian
{
#ifdef RESULTTEST
  friend class AleResultTest;
#endif

public:

  AleLinear(Teuchos::RefCountPtr<DRT::Discretization> actdis,
            Teuchos::RefCountPtr<LINALG::Solver> solver,
            Teuchos::RefCountPtr<ParameterList> params,
            Teuchos::RefCountPtr<IO::DiscretizationWriter> output);

  void Integrate();

  virtual void PrepareTimeStep();
  virtual void Solve();
  virtual void Update();
  virtual void Output();

  virtual bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  virtual bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

  void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RefCountPtr<Epetra_Vector> idisp);

  const DRT::Discretization& Discretization() const { return *discret_; };

  virtual Teuchos::RefCountPtr<Epetra_Vector> ExtractDisplacement();

private:

  Teuchos::RefCountPtr<DRT::Discretization> discret_;
  Teuchos::RefCountPtr<LINALG::Solver> solver_;
  Teuchos::RefCountPtr<ParameterList> params_;
  Teuchos::RefCountPtr<IO::DiscretizationWriter> output_;

  int step_;
  int nstep_;
  double time_;
  double maxtime_;
  double dt_;

  int maxentriesperrow_;

  Teuchos::RefCountPtr<Epetra_CrsMatrix> sysmat_;

  Teuchos::RefCountPtr<Epetra_Vector>    dirichtoggle_;
  //Teuchos::RefCountPtr<Epetra_Vector>    zeros_;

  Teuchos::RefCountPtr<Epetra_Vector>    dispnp_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispn_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispnm_;

  Teuchos::RefCountPtr<Epetra_Vector>    residual_;

  Teuchos::RefCountPtr<Epetra_Map> imeshmap_;
  Teuchos::RefCountPtr<Epetra_Import> extractor_;

private:

  void EvaluateElements();

  //Epetra_Vector &x_;
  bool haveF_;
  bool haveJacobian_;
};
}

#endif
#endif
#endif
