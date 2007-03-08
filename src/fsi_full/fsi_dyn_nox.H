
#ifndef FSI_DYN_NOX_H
#define FSI_DYN_NOX_H

#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Import.h>
#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

#include <vector>

extern "C"
{
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//! FSI interface equation ready to solve with NOX
/*!

This is the huge fat NOX interface class for FSI.  The purpose of this
class is to implement the NOX::Epetra::Interface::Required interface
and thus implement the computeF method. The idea behind all this is to
use a nonlinear solver (NOX) on the FSI interface problem. That is the
FSI problem is cast into the form
\f[
R(d_\Gamma) = d_\Gamma - S(F(d_\Gamma))
\f]
where \f$ d_\Gamma \f$ represents the displacements at the FSI
interface \f$ \Gamma \f$. The solvers \f$ S() \f$ and \f$ F() \f$ are
therefore the structural and fluid solver condensed to interface
degrees of freedom. (The Schur complements. Of course, internally the
whole field has to be solved. And it is!)

The above equation can be subjected to a nonlinear solver. NOX gives
us plenty to choose from. Some are supported. You select an solver by
appropriate parameter list settings in the ccarat input file. That is
NOX reads its parameters from a Teuchos::ParameterList object. Such a
object is read from the dat file. The following list describes how the
final parameter list has to look.

<b>Input Parameters</b>

The most important thing is to chose a Jacobian approximation.

<li> "Jacobian" = "None"
<li> "Jacobian" = "Matrix Free"
<li> "Jacobian" = "Finite Difference"
<li> "Jacobian" = "Finite Difference 1"

Secondly you can select a preconditioner.

<li> "Preconditioner" = "None"
<li> "Preconditioner" = "Finite Difference"
<li> "Preconditioner" = "Finite Difference 1"

Solvers are selected using:

"Direction"

<li> "Method" = "Newton"
<li> "Method" = "Modified-Newton"
<li> "Method" = "Quasi-Newton"
<li> "Method" = "Broyden"
<li> "Method" = "NonlinearCG"

After the solve a line search (relaxation) step might be in order:

"Line Search"

<li> "Method" = "Full Step"
<li> "Method" = "NonlinearCG"
<li> "Method" = "Aitken"

Please note that NOX specifies a bunch of (possible) options you can
chose. Both for its nonlinear solvers and for the inner linear solver
and preconditioner. These settings are interrelated.

<b>References</b>

  <ul>

  <li>Aitken and steepest descent relaxation are the traditional way
  to do thing. Due to Wall and Mok.

\author u.kue
\date 03/07
*/
class FSI_InterfaceProblem : public NOX::Epetra::Interface::Required
{
public:

//! Constructor
  FSI_InterfaceProblem(Epetra_Comm& Comm,
                       FSI_STRUCT_WORK* struct_work,
                       FSI_FLUID_WORK* fluid_work,
                       FSI_ALE_WORK* ale_work);

  //! Compute the function, F, given the specified input vector x.  Returns true if computation was successful.
  /*!
    This is the place where we do one FSI cycle (ale, fluid,
    structure).

    \param x         (in): prescribed interface displacements
    \param F        (out): new residuum
    \param fillFlag  (in): situation in which residuum is needed
   */
  virtual bool computeF(const Epetra_Vector& x,
                        Epetra_Vector& F,
			const FillType fillFlag);

  //! current interface displacements extracted from the nodes
  Teuchos::RefCountPtr<Epetra_Vector> soln();

  //! the huge time loop
  void timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface);

private:

  //! FSI control variables
  FSI_DYNAMIC *fsidyn_;

  //! distribution of interface dofs extracted from the structural mesh distribution
  Teuchos::RefCountPtr<Epetra_Map> StandardMap_;

  //! redundant distribution of interface dofs
  Teuchos::RefCountPtr<Epetra_Map> redundantmap_;

  //! interface displacements distributed
  Teuchos::RefCountPtr<Epetra_Vector> soln_;

  //! redundant version of interface displacements
  Teuchos::RefCountPtr<Epetra_Vector> redundantsol_;

  //! redundant version of interface residuum
  Teuchos::RefCountPtr<Epetra_Vector> redundantf_;

  //! connection of interface dofs for colored finite differences
  Teuchos::RefCountPtr<Epetra_CrsGraph> rawGraph_;

  //! importer to make interface displacements redundant
  Teuchos::RefCountPtr<Epetra_Import> importredundant_;

  //! traditional structure variables
  FSI_STRUCT_WORK* struct_work_;

  //! traditional fluid variables
  FSI_FLUID_WORK* fluid_work_;

  //! traditional ale variables
  FSI_ALE_WORK* ale_work_;

  //! MPI_COMM_WORLD wrapped for trilinos
  Epetra_Comm& Comm_;

  //! counters on how many times the residuum was called in a time step
  /*!
    NOX knows different types of residuum calls depending on
    circumstances (normal, finite difference, matrix free
    jacobi). It is possible to do approximations depending on the
    type.
   */
  std::vector<int> counter_;
};


#endif
#endif
#endif
