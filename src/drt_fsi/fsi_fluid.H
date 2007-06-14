
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_FLUID_H
#define FSI_FLUID_H

#include "../drt_fluid/fluidimplicitintegration.H"

#include <Epetra_Map.h>
#include <Teuchos_RefCountPtr.hpp>


namespace FSI
{

/// Fluid algorithm enhanced for FSI
class Fluid : public FluidImplicitTimeInt
{
public:
  Fluid(RefCountPtr<DRT::Discretization> dis,
        RefCountPtr<LINALG::Solver> solver,
        RefCountPtr<ParameterList> params,
        RefCountPtr<DiscretizationWriter> output);

  /// announce a new time step
  //virtual void NewStep(int step, double time);

  void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);

  Teuchos::RefCountPtr<Epetra_Vector> ExtractInterfaceForces();

  void ApplyInterfaceVelocities(Teuchos::RefCountPtr<Epetra_Vector> ivel);

  void SetMeshMap(Teuchos::RefCountPtr<Epetra_Map> mm);

  void ApplyMeshDisplacement(Teuchos::RefCountPtr<Epetra_Vector> fluiddisp);

  const DRT::Discretization& Discretization() const { return *discret_; };

  /// linear fluid solve with just a interface load
  Teuchos::RefCountPtr<Epetra_Vector> RelaxationSolve(Teuchos::RefCountPtr<Epetra_Vector> ivel);

private:

  Teuchos::RefCountPtr<Epetra_Map> ivelmap_;
  Teuchos::RefCountPtr<Epetra_Export> extractor_;

  Teuchos::RefCountPtr<Epetra_Map> meshmap_;
  Teuchos::RefCountPtr<Epetra_Export> meshextractor_;

  Teuchos::RefCountPtr<Epetra_Vector> relax_;

  RefCountPtr<LINALG::Solver> solver_;
  RefCountPtr<ParameterList> params_;
  RefCountPtr<DiscretizationWriter> output_;
};

}

#endif
#endif
#endif
