/*----------------------------------------------------------------------*/
/*!
\brief Adjoint time integration for a quasi static hyperelasticity

\level 3

\maintainer Sebastian Brandstaeter
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
}  // namespace DRT

namespace LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
  class SparseOperator;
  class BlockSparseMatrixBase;
}  // namespace LINALG

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
    TimIntAdjoint(Teuchos::RCP<DRT::Discretization> dis);

    //! destructor
    virtual ~TimIntAdjoint() {}

    //! interface to the calling optimization routine to set up the rhs and the primal variables
    /*! \param In
        rhs - rigth hand side for the adjoint equation.
        \param In
        mtime - vector with time for each entry in rhs.
        \param In
        dis - primal solution to the adjoint problem.
        \param In
        mtime - vector with time for each entry in the primal solution.
        \param In
    */
    virtual void SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs, std::vector<double> mtime,
        Teuchos::RCP<Epetra_MultiVector> dis, std::vector<double> time);

    //! do reverse time integration
    void Integrate();

    //! dual variables corresponding to the standard adjoint equation
    Teuchos::RCP<Epetra_MultiVector> ExtractSolution() { return disdual_; }

    //! dual variables corresponding to the "prestressed adjoint equation"
    virtual Teuchos::RCP<Epetra_MultiVector> ExtractPrestressSolution() { return Teuchos::null; }

   protected:
    //! intro
    void PrintLogo();

    //! solve the linear system
    virtual void Solve();

    //! evaluate stiffness matrix
    void EvaluateStiff();

    //! set up things needed for this timestep
    void PrepareStep();

    //! update step and time
    void Update();

    //! update state
    virtual void UpdateStepState();

    //! set time step history
    void SetTimeStepHistory(int timestep);

    //! create a linear solver according to sdyn section
    void CreateLinearSolver(const Teuchos::ParameterList& sdyn);

    //! create the writer
    void CreateWriter();

    //! write output
    void OutputStep();

    // DBC Map
    void GetDBCMap();

    //! return (rotatory) transformation matrix of local co-ordinate systems
    Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

    //! find step in the rhs vector according to this time
    int FindStep(double time);

    //! print information of time steps to screen
    int printtoscreen_;

    //! the discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! a writer for the output
    Teuchos::RCP<IO::DiscretizationWriter> writer_;


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

    //! to stepn_ corresponding step of the rhs
    int steprhsn_;

    //! timesteps coming from the primal problem
    std::vector<double> time_;

    //! dt (from forward problem)
    double dt_;

    //! timesteps for which there is a rhs not zero
    std::vector<double> mtime_;

    //! current time
    double timen_;

    //! "write step" counter
    int writestep_;

    //! "write time"
    double writetime_;

    //! rhs and primal variables set?
    bool isinit_;
  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif
