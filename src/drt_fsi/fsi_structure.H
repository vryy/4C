/*----------------------------------------------------------------------*/
/*!
\file

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
#ifdef TRILINOS_PACKAGE

#ifndef FSI_STRUCTURE_H
#define FSI_STRUCTURE_H

#include "../drt_structure/strugenalpha.H"

#include <Epetra_Map.h>
#include <Teuchos_RefCountPtr.hpp>


namespace FSI
{

/// Structural dynamics algorithm enhanced for FSI
/*!
  Dirichlet-Neumann FSI needs to apply interface forces to the
  structure and takes interface displacements after the structural
  solve. The coupling algorithm is implemented in
  DirichletNeumannCoupling.

  \note the coupling of this class to its base StruGenAlpha is really
  close. You have to understand StruGenAlpha before you look here.

  \author u.kue
  \date 10/07
 */
class Structure : public StruGenAlpha
{
public:
  Structure(Teuchos::RefCountPtr<ParameterList> params,
            Teuchos::RefCountPtr<DRT::Discretization> dis,
            Teuchos::RefCountPtr<LINALG::Solver> solver,
            Teuchos::RefCountPtr<IO::DiscretizationWriter> output);

  /// announce a new time step to the solver
  /*!
    Start the time step by a predictor call.
   */
  void PrepareTimeStep();

  //! @name Interface Map
  //! the map of all interface dofs will be set from the outside

  Teuchos::RefCountPtr<Epetra_Map> InterfaceMap() { return idispmap_; }
  void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);

  //@}

  //! @name Extract interface values

  /// extract displacements
  Teuchos::RefCountPtr<Epetra_Vector> ExtractInterfaceDisplacement();

  /// extract velocities
  Teuchos::RefCountPtr<Epetra_Vector> ExtractInterfaceVel();

  /// extract accelerations
  Teuchos::RefCountPtr<Epetra_Vector> ExtractInterfaceAcc();

  //@}

  //! @name Apply interface forces

  /// apply interface forces to structural solver
  /*!
    This prepares a new solve of the structural field within one time
    step. The middle values are newly created.

    \note This is not yet the most efficient implementation.
   */
  void ApplyInterfaceForces(Teuchos::RefCountPtr<Epetra_Vector> iforce);

  //@}

  /// discretization access
  const DRT::Discretization& Discretization() const { return discret_; };

  /// linear structure solve with just a interface load
  /*!
    The very special solve done in steepest descent relaxation
    calculation (and matrix free Newton Krylov).

    \note Can only be called after a valid structural solve.
   */
  Teuchos::RefCountPtr<Epetra_Vector> RelaxationSolve(Teuchos::RefCountPtr<Epetra_Vector> iforce);

private:

  /// element call and effective stiffness calculation
  void CalculateStiffness();

  /// interface dof map
  Teuchos::RefCountPtr<Epetra_Map> idispmap_;

  /// communicator between interface dof map and field dof map
  Teuchos::RefCountPtr<Epetra_Import> extractor_;

  /// copy of external forces needed to reapply the interface forces
  Teuchos::RefCountPtr<Epetra_Vector> fextncopy_;

  //! @name local copies of input parameters
  Teuchos::RefCountPtr<ParameterList> params_;
  Teuchos::RefCountPtr<LINALG::Solver> solver_;
  Teuchos::RefCountPtr<IO::DiscretizationWriter> output_;
  //@}
};

}

#endif
#endif
#endif
