
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_DIRICHLETNEUMANN_H
#define FSI_DIRICHLETNEUMANN_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <AztecOO.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include <vector>

#include "fsi_structure.H"
#include "fsi_fluid.H"
#include "fsi_ale.H"
#include "fsi_coupling.H"


namespace FSI {


/// Dirichlet Neumann coupled interface system
class DirichletNeumannCoupling : public NOX::Epetra::Interface::Required
{
public:

  DirichletNeumannCoupling(Epetra_Comm& comm);

  /// setup FSI algorithm
  void Setup();

  /// outer level FSI time loop
  void Timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface);

  /// compute FSI interface residuum d - S^-1(F(d))
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

private:

  Teuchos::RefCountPtr<Epetra_Vector> FluidOp(Teuchos::RefCountPtr<Epetra_Vector> idisp);

  Teuchos::RefCountPtr<Epetra_Vector> StructOp(Teuchos::RefCountPtr<Epetra_Vector> iforce);

  Teuchos::RefCountPtr<Epetra_Vector> InterfaceVelocity(Teuchos::RefCountPtr<Epetra_Vector> idispnp);

  Teuchos::RefCountPtr<Epetra_Vector> StructToAle(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> StructToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> FluidToStruct(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> AleToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv);

  //! current interface displacements
  Teuchos::RefCountPtr<Epetra_Vector> InterfaceDisp();

  /// setup structure algorithm
  void SetupStructure();

  /// setup fluid algorithm
  void SetupFluid();

  /// setup ale algorithm
  void SetupAle();

  void PrepareTimeStep();

  void Update();

  void Output();

  Epetra_Comm& comm_;

  RefCountPtr<Structure> structure_;
  RefCountPtr<Fluid>     fluid_;
  RefCountPtr<Ale>       ale_;

  /// coupling of structure and fluid at the interface
  Coupling               coupsf_;

  /// coupling of structure and ale at the interface
  Coupling               coupsa_;

  /// coupling of fluid and ale
  Coupling               coupfa_;

  //! main nonlinear solver object that contains everything else
  Teuchos::RefCountPtr<NOX::Solver::Manager> solver_;

  Teuchos::RefCountPtr<Epetra_Vector> idispn_;

  //! counters on how many times the residuum was called in a time step
  /*!
    NOX knows different types of residuum calls depending on
    circumstances (normal, finite difference, matrix free
    jacobi). It is possible to do approximations depending on the
    type.
   */
  std::vector<int> counter_;

  //! number of residuum calculations per nonlinear solve in one time step
  std::vector<int> linsolvcount_;

  //! we couple based on displacements
  bool displacementcoupling_;

  int step_;
  int nstep_;
  double time_;
  double maxtime_;
  double dt_;
};

}

#endif
#endif
#endif
