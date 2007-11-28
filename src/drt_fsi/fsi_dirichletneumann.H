/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann.H

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

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
#include "fsi_coupling_mortar.H"


/// FSI: Fluid Structure Interaction
/*!
  This is the place where we bring fluid and structure together.
 */
namespace FSI {


/// Dirichlet Neumann coupled interface system
/*!
  This is the algorithm class of partitioned FSI problems. Here we
  do the time loop and the coupling between fields. The fields
  themselves are solved using appropiate field algorithms (that are
  used for standalone solvers as well.) The FSI interface problem is
  solved using NOX.

  Many different things come together in this class. First and
  foremost there is the Timeloop() method, that contains the overall
  FSI time stepping scheme. The time loop is build around the FSI
  interface problem, that is the problem of finding the interface
  displacements that satisfy the coupled problem. The time loop knows
  nothing about the fields themselves.

  Inside the time loop the interface problem is solved using NOX. To
  do so Timeloop() needs to know an object which defines the nonlinear
  residual of the FSI problem. This is always a RefCountPtr to the
  DirichletNeumannCoupling object itself!

  So the second part of this class consists of the interface residual
  evaluation computeF(). This method does one FSI cycle, that is one
  solve of all three participating fields. But for sake of clarity
  this cycle is expressed via the two operator methods FluidOp() and
  StructOp().

  FluidOp() takes a interface displacement, applies it to the ale
  field, solves the ale field, calculates the interface velocities,
  applies them to the fluid field, solves the fluid field on the
  newly deformed fluid mesh and returns the interface forces.

  StructOp() takes interface forces, applies them to the structural
  field, solves the field and returns the interface displacements.

  This coupling process build on the available field solvers. However,
  the independant parallel distribution of the fields complicates the
  exchange of coupling information. Therefore three instances of the
  Coupling class are used that couple those fields. On top of these
  there are helper methods StructToAle(), StructToFluid(),
  FluidToStruct() and AleToFluid() to easily exchange distributed
  interface vectors between fields.

  The FSI algorithm requires repeated evaluations of the interface
  residual via computeF(). So the field solvers themselves must be
  clean, subsequent calls with the same interface input must yield the
  same results. The time stepping therefore needs the further methods
  PrepareTimeStep() to start a new time step as wellas Update() and
  Output() to finish the current step, save the result and write the
  files.

  \author u.kue
  \date 06/07
 */
class DirichletNeumannCoupling : public NOX::Epetra::Interface::Required
{
public:

  /// setup FSI algorithm
  DirichletNeumannCoupling(Epetra_Comm& comm);

  /// read restart data
  void ReadRestart(int step);

  /// outer level FSI time loop
  void Timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface);

  /// compute FSI interface residuum S^{-1}(F(d)) - d
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// structural field solver
  Structure& StructureField() { return *structure_; }

  /// fluid field solver
  Fluid& FluidField() { return *fluid_; }

  /// ale field solver
  Ale& AleField() { return *ale_; }

private:

  ///
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>
  CreateLinearSystem(ParameterList& nlParams,
                     const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface,
                     NOX::Epetra::Vector& noxSoln,
                     Teuchos::RefCountPtr<NOX::Utils> utils);

  ///
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> CreateStatusTest(ParameterList& nlParams,
                                                                Teuchos::RefCountPtr<NOX::Epetra::Group> grp);

  /// interface fluid operator
  Teuchos::RefCountPtr<Epetra_Vector> FluidOp(Teuchos::RefCountPtr<Epetra_Vector> idisp,
                                              const FillType fillFlag);

  /// interface structural operator
  Teuchos::RefCountPtr<Epetra_Vector> StructOp(Teuchos::RefCountPtr<Epetra_Vector> iforce,
                                               const FillType fillFlag);

  /// interface velocity calculation given the interface displacements
  Teuchos::RefCountPtr<Epetra_Vector> InterfaceVelocity(Teuchos::RefCountPtr<Epetra_Vector> idispnp);

  Teuchos::RefCountPtr<Epetra_Vector> StructToAle(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> StructToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> FluidToStruct(Teuchos::RefCountPtr<Epetra_Vector> iv);
  Teuchos::RefCountPtr<Epetra_Vector> AleToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv);

  /// current interface displacements
  Teuchos::RefCountPtr<Epetra_Vector> InterfaceDisp();

  /// current interface forces
  Teuchos::RefCountPtr<Epetra_Vector> InterfaceForce();

  /// setup structure algorithm
  void SetupStructure();

  /// setup fluid algorithm
  void SetupFluid();

  /// setup ale algorithm
  void SetupAle();

  /// start a new time step
  void PrepareTimeStep();

  /// take current results for converged and save for next time step
  void Update();

  /// write output
  void Output();

  /// comunication (mainly for screen output)
  Epetra_Comm& comm_;

  /// structural field solver
  RefCountPtr<Structure> structure_;

  /// fluid field solver
  RefCountPtr<Fluid>     fluid_;

  /// ale field solver
  RefCountPtr<Ale>       ale_;

  /// coupling of structure and fluid at the interface
  Coupling               coupsf_;

  /// coupling of structure and ale at the interface
  Coupling               coupsa_;

  /// coupling of fluid and ale
  Coupling               coupfa_;

  /// coupling of structure and fluid at the interface, with mortar. Experimental!
  CouplingMortar         coupsfm_;

  /// interface displacement from time step begin
  Teuchos::RefCountPtr<Epetra_Vector> idispn_;

  //! connection of interface dofs for finite differences
  Teuchos::RefCountPtr<Epetra_CrsGraph> rawGraph_;

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

  /// print parameters and stuff
  /*!
    \warning This variable is only valid during when the time loop runs.
   */
  Teuchos::RefCountPtr<NOX::Utils> utils_;

  int step_;
  int nstep_;
  double time_;
  double maxtime_;
  double dt_;
  int mfresitemax_;

  /// nodes at the fluid-structure interface match
  bool matchingnodes_;
};

}

#endif
#endif
