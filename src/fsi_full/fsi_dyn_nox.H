
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
#include <Epetra_Export.h>
#include <AztecOO.h>
#include "Teuchos_ParameterList.hpp"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FSI_InterfaceProblem : public NOX::Epetra::Interface::Required
{
public:
  FSI_InterfaceProblem(Epetra_Comm& Comm,
                       FSI_STRUCT_WORK* struct_work,
                       FSI_FLUID_WORK* fluid_work,
                       FSI_ALE_WORK* ale_work);

  //! Compute the function, F, given the specified input vector x.  Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x,
                        Epetra_Vector& F,
			const FillType fillFlag);

  Teuchos::RefCountPtr<Epetra_Vector> soln();

  void timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface);

private:

  FSI_DYNAMIC *fsidyn_;

  Teuchos::RefCountPtr<Epetra_Map> StandardMap_;
  Teuchos::RefCountPtr<Epetra_Map> redundantmap_;
  Teuchos::RefCountPtr<Epetra_Vector> soln_;
  Teuchos::RefCountPtr<Epetra_Vector> redundantsol_;
  Teuchos::RefCountPtr<Epetra_CrsGraph> rawGraph_;

  Teuchos::RefCountPtr<Epetra_Export> exportredundant_;

  FSI_STRUCT_WORK* struct_work_;
  FSI_FLUID_WORK* fluid_work_;
  FSI_ALE_WORK* ale_work_;

  Epetra_Comm& Comm_;
};


#endif
#endif
#endif
