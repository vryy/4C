
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_STRUCTURE_H
#define FSI_STRUCTURE_H

#include "../drt_lib/strugenalpha.H"

#include <Epetra_Map.h>
#include <Teuchos_RefCountPtr.hpp>


namespace FSI
{

/// Structural dynamics algorithm enhanced for FSI
class Structure : public StruGenAlpha
{
public:
  Structure(Teuchos::RefCountPtr<ParameterList> params,
            Teuchos::RefCountPtr<DRT::Discretization> dis,
            Teuchos::RefCountPtr<LINALG::Solver> solver,
            Teuchos::RefCountPtr<DiscretizationWriter> output);

  /// announce a new time step
  //virtual void NewStep(int step, double time);

  void PrepareTimeStep();

  void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);

  Teuchos::RefCountPtr<Epetra_Vector> ExtractInterfaceDisplacement();

  void ApplyInterfaceForces(Teuchos::RefCountPtr<Epetra_Vector> iforce);

  const DRT::Discretization& Discretization() const { return discret_; };

private:

  void CalculateStiffness();

  Teuchos::RefCountPtr<Epetra_Map> idispmap_;
  Teuchos::RefCountPtr<Epetra_Export> extractor_;
  //Teuchos::RefCountPtr<Epetra_Import> inserter_;

  Teuchos::RefCountPtr<Epetra_Vector> fextncopy_;

  Teuchos::RefCountPtr<ParameterList> params_;
  //Teuchos::RefCountPtr<DRT::Discretization> dis_;
  Teuchos::RefCountPtr<LINALG::Solver> solver_;
  Teuchos::RefCountPtr<DiscretizationWriter> output_;
};

}

#endif
#endif
#endif
