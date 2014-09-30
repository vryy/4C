/*----------------------------------------------------------------------*/
/*!
 \file timint_adjoint.H
 \brief Adjoint time integration for a quasi static hyperelasticity

 <pre>
 Maintainer: Sebastian Kehl
 kehl@mhpc.mw.tum.de
 089 - 289-10361
 </pre>
 */

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef TIMINT_ADJOINT_H
#define TIMINT_ADJOINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"

// forward declarations
namespace DRT
{
class Discretization;

namespace UTILS
{
class LocsysManager;
}
}

namespace LINALG
{
class Solver;
class MapExtractor;
class SparseMatrix;
class SparseOperator;
class BlockSparseMatrixBase;
}

namespace IO
{
class DiscretizationWriter;
}

/*----------------------------------------------------------------------*/
namespace STR
{

class TimIntAdjoint
{

public:

  //! constructor
  TimIntAdjoint(Teuchos::RCP<DRT::Discretization> dis, std::vector<double> time);

  //! destructor
  virtual ~TimIntAdjoint() {}

  // interface to the calling optimization routine to set up the rhs and the primal variables
  void SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs, Teuchos::RCP<Epetra_MultiVector> dis);

  //! do reverse time integration
  void Integrate();

  // result vector
  Teuchos::RCP<Epetra_MultiVector> ExtractSolution() {return disdual_;}

private:

  //! intro
  void PrintLogo();

  //! solve the linear system
  void Solve();

  // evaluate stiffness matrix
  void EvaluateStiff();

  // set up things needed for this timestep
  void PrepareStep();

  // update step and time
  void UpdateStep();

  // create a linear solver according to sdyn section
  void CreateLinearSolver(const Teuchos::ParameterList& sdyn);

  // DBC Map
  void GetDBCMap();

  //! return (rotatory) transformation matrix of local co-ordinate systems
  Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

  //! print information of time steps to screen
  int printtoscreen_;

  //! the discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! the dofrowmap of the discretization
  const Epetra_Map* dofrowmap_;

  //! a solver
  Teuchos::RCP<LINALG::Solver> solver_;

  //! to apply Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

  //! map to the dirichlet boundary
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  //! DBC toggle
  Teuchos::RCP<Epetra_Vector> dbctoggle_;

  //! primal variables (displacements)
  Teuchos::RCP<Epetra_MultiVector> dis_;

  //! dual variables (lagrange multiplier)
  Teuchos::RCP<Epetra_MultiVector> disdual_;

  //! rhs for the adjoints
  Teuchos::RCP<Epetra_MultiVector> rhs_;

  //! dis_ at step n
  Teuchos::RCP<Epetra_Vector> disn_;

  //! disdualn_ at step n
  Teuchos::RCP<Epetra_Vector> disdualn_;

  //! disdualn_ at step n
  Teuchos::RCP<Epetra_Vector> rhsn_;

  //! stiffness matrix
  Teuchos::RCP<LINALG::SparseOperator> stiff_;

  //! stiffness matrix with neumann terms
  Teuchos::RCP<LINALG::SparseOperator> stiffn_;

  //! a vector with zeros
  Teuchos::RCP<Epetra_Vector> zeros_;

  //! number of simulation steps in the primal problem
  int msteps_;

  //! time step at time n
  int stepn_;

  //! timesteps coming from the primal problem
  std::vector<double> time_;

  //! current time
  double timen_;

  //! rhs and primal variables set?
  bool isinit_;

};

}

/*----------------------------------------------------------------------*/
#endif
