/*----------------------------------------------------------------------*/
/*!
\file timint_adjoint.H
\brief Time integration for an hyperleastic quasi static adjoint problem

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef TIMINT_ADJOINT_H
#define TIMINT_ADJOINT_H

/*----------------------------------------------------------------------*/
/* headers */
#include <string>
#include <iostream>
#include <fstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_Time.h>


// forward declarations
namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    class LocsysManager;
  }
}

namespace UTILS
{
  class ConstrManager;
  class ConstraintSolver;
  class SurfStressManager;
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
//! Structural Dynamics
namespace STR
{


  class TimIntAdjoint
  {

  public:

    //! Constructor
    TimIntAdjoint(Teuchos::RCP<DRT::Discretization> dis);

    //! Empty constructor
    TimIntAdjoint() { ; }

    // interface to the calling optimization routine to set up the rhs and the primal variables
    void SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs,
                      Teuchos::RCP<Epetra_MultiVector> dis);

    //! Do time integration of multiple steps
    void Integrate();

    // result vector
    Teuchos::RCP<Epetra_MultiVector> ExtractSolution() {return disdual_;};


  private:

    //! Intro
    void PrintLogo();

    //! Solve the linear system
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

    //! print information of time steps to screen
    int printtoscreen_;

    //! some funny stuff
    Teuchos::RCP<DRT::Discretization> discret_;

    //! a dofrowmap of the discretization
    const Epetra_Map* dofrowmap_;

    //! a solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //! Map to the dirichlet boundary
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

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

    //! a vector with zeros
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! number of simulation steps in the primal problem
    int msteps_;

    //! time step at time n
    int stepn_;

    //! rhs and primal variables set?
    bool isinit_;


  };

}

/*----------------------------------------------------------------------*/
#endif
