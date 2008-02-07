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

#ifndef FSI_STRUCTURE_H
#define FSI_STRUCTURE_H

#include "../drt_structure/strugenalpha.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include "fsi_utils.H"


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
  Structure(Teuchos::RCP<ParameterList> params,
            Teuchos::RCP<DRT::Discretization> dis,
            Teuchos::RCP<LINALG::Solver> solver,
            Teuchos::RCP<IO::DiscretizationWriter> output);

  /// announce a new time step to the solver
  /*!
    Start the time step by a predictor call.
   */
  void PrepareTimeStep();

  //! @name Interface Map
  //! the map of all interface dofs will be set from the outside

  Teuchos::RCP<const Epetra_Map> InterfaceMap() { return interface_.CondMap(); }
  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  //@}

  //! @name Extract interface values

  /// extract displacements
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

  /// extract displacements
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceDisplacement();

  /// Predictor for interface displacements
  Teuchos::RCP<Epetra_Vector> PredictInterfaceDisplacement();

  //@}

  //! @name Apply interface forces

  /// apply interface forces to structural solver
  /*!
    This prepares a new solve of the structural field within one time
    step. The middle values are newly created.

    \note This is not yet the most efficient implementation.
   */
  void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);

  //@}

  /// discretization access
  const DRT::Discretization& Discretization() const { return discret_; };

  /// solve the current system
  void Solve();

  /// linear structure solve with just a interface load
  /*!
    The very special solve done in steepest descent relaxation
    calculation (and matrix free Newton Krylov).

    \note Can only be called after a valid structural solve.
   */
  Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

  /// communication object at the interface
  const LINALG::MapExtractor& Interface() const { return interface_; }

private:

  /// element call and effective stiffness calculation
  void CalculateStiffness();

  /// the interface map setup for interface <-> full translation
  LINALG::MapExtractor interface_;

  /// copy of external forces needed to reapply the interface forces
  Teuchos::RCP<Epetra_Vector> fextncopy_;

  //! @name local copies of input parameters
  Teuchos::RCP<ParameterList> params_;
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;
  //@}
};

}

#endif
#endif
