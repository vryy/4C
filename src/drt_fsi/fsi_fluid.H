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

#ifndef FSI_FLUID_H
#define FSI_FLUID_H

#include "../drt_fluid/fluidimplicitintegration.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include "fsi_utils.H"


namespace FSI
{

/// Fluid algorithm enhanced for FSI
/*!
  Dirichlet-Neumann FSI needs to apply interface velocities to the
  fluid and takes interface forces after the fluid solve. At the same
  time interface the displacements lead to a mesh deformation and thus
  to a grid velocity. The coupling algorithm is implemented in
  DirichletNeumannCoupling.

*/
class Fluid : public FluidImplicitTimeInt
{
public:
  Fluid(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output);

  //! @name Interface Map
  //! the map of all interface dofs will be set from the outside

  /// the interface map contains all velocity dofs at the interface
  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  /// the mesh map contains all velocity dofs that are covered by an ALE node
  void SetMeshMap(Teuchos::RCP<Epetra_Map> mm);

  //@}

  //! @name Extract interface forces

  /// After the fluid solve we need the forces at the FSI interface.
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

  //@}

  //! @name Extract interface values

  /// at the interface the velocity is prescribed as a Dirichlet condition
  void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

  //@}

  //! @name Extract mesh values

  /// tell the mesh displacement to the fluid solver
  void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp);

  /// tell the mesh velocity to the fluid solver
  void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel);

  //@}

  /// discretization access
  const DRT::Discretization& Discretization() const { return *discret_; };

  /// linear fluid solve with just a interface load
  /*!
    The very special solve done in steepest descent relaxation
    calculation (and matrix free Newton Krylov).

    \note Can only be called after a valid fluid solve.
   */
  Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

  /// integrate FSI interface shape functions
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();


  //! @name Number of Newton iterations
  //! For simplified FD MFNK solve we want to temporally limit the
  /// number of Newton steps inside the fluid solver

  int Itemax() const;
  void SetItemax(int itemax);

  //@}

  //! @name Access to current time step for restart

  double time() const { return time_; }
  int step() const { return step_; }

  //@}

  /// communication object at the interface
  const LINALG::MapExtractor& Interface() const { return interface_; }

private:

  /// the interface map setup for interface <-> full translation
  LINALG::MapExtractor interface_;

  /// ALE dof map
  LINALG::MapExtractor meshmap_;

  /// local variable for relaxation calculation
  Teuchos::RCP<Epetra_Vector> relax_;

  /// grid displacement
  Teuchos::RCP<Epetra_Vector> griddisp_;

  //! @name local copies of input parameters
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<ParameterList> params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;
  //@}
};

}

#endif
#endif
