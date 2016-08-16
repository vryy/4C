/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann.H

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362

\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_H
#define FSI_DIRICHLETNEUMANN_H

#include "fsi_partitioned.H"

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
  residual of the FSI problem. This is always a Teuchos::RCP to the
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
  the independent parallel distribution of the fields complicates the
  exchange of coupling information. Therefore three instances of the
  Coupling class are used that couple those fields. On top of these
  there are helper methods StructToAle(), StructToFluid(),
  FluidToStruct() and AleToFluid() to easily exchange distributed
  interface vectors between fields.

  The FSI algorithm requires repeated evaluations of the interface
  residual via computeF(). So the field solvers themselves must be
  clean, subsequent calls with the same interface input must yield the
  same results. The time stepping therefore needs the further methods
  PrepareTimeStep() to start a new time step as well as Update() and
  Output() to finish the current step, save the result and write the
  files.

  \author u.kue
  \date 06/07
 */
class DirichletNeumann : public Partitioned
{
public:

  /// constructor
  explicit DirichletNeumann(const Epetra_Comm& comm);

  /// setup this object
  virtual void Setup();

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag);

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

private:

  //! we couple based on displacements
  bool displacementcoupling_;

};

}

#endif
