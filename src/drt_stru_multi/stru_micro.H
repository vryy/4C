#if 0

#ifndef STRU_MICRO_H
#define STRU_MICRO_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "../drt_lib/strugenalpha.H"

#include <Epetra_Map.h>
#include <Teuchos_RefCountPtr.hpp>

using namespace Teuchos;

namespace StructureMulti
{


/// Structural dynamics algorithm enhanced for microstructure simulations
/// in case of multiscale analyses

class Micro : public StruGenAlpha
{
public:
  Micro(RefCountPtr<ParameterList> params,
        RefCountPtr<DRT::Discretization> dis,
        RefCountPtr<LINALG::Solver> solver,
        RefCountPtr<IO::DiscretizationWriter> output);

  void DetermineToggle();

  void EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd);

  void SetOldDisp(RefCountPtr<Epetra_Vector> disp);

  void ReturnNewDisp();

  void ClearDisp();



private:

  RefCountPtr<ParameterList> params_;
  RefCountPtr<LINALG::Solver> solver_;
  RefCountPtr<IO::DiscretizationWriter> output_;

};

}

#endif
#endif
#endif

#endif
