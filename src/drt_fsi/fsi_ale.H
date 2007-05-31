
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_ALE_H
#define FSI_ALE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../io/io_drt.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Required.H>

namespace FSI {


/// General ALE interface
class Ale
{
public:
  virtual ~Ale() {}
};


/// linear ALE algorithm
class AleLinear : public NOX::Epetra::Interface::Required,
                  public NOX::Epetra::Interface::Jacobian
{
#ifdef RESULTTEST
  friend class AleResultTest;
#endif

public:

  AleLinear(Teuchos::RefCountPtr<DRT::Discretization> actdis,
            LINALG::Solver&       solver,
            ParameterList&        params,
            DiscretizationWriter& output);

  void Integrate();

  virtual bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  virtual bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

private:

  Teuchos::RefCountPtr<DRT::Discretization> discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  DiscretizationWriter&         output_;

  int step_;
  int nstep_;
  double time_;
  double maxtime_;
  double dt_;

  int maxentriesperrow_;

  Teuchos::RefCountPtr<Epetra_CrsMatrix> sysmat_;

  Teuchos::RefCountPtr<Epetra_Vector>    dirichtoggle_;
  Teuchos::RefCountPtr<Epetra_Vector>    zeros_;

  Teuchos::RefCountPtr<Epetra_Vector>    dispnp_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispn_;
  Teuchos::RefCountPtr<Epetra_Vector>    dispnm_;

  Teuchos::RefCountPtr<Epetra_Vector>    residual_;

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
