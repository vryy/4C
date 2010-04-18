/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_base_algorithm.H

\brief Fluid Base Algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_BASE_ALGORITHM_H
#define ADAPTER_FLUID_BASE_ALGORITHM_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "adapter_fluid.H"


// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// fluid field solver
  class FluidBaseAlgorithm
  {
  public:

    /// constructor
    explicit FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn, bool isale);

    /// virtual destructor to support polymorph destruction
    virtual ~FluidBaseAlgorithm();

    /// fluid field solver
    Fluid& FluidField() { return *fluid_; }
    /// return instance of fluid time integration class
    const Fluid& FluidField() const { return *fluid_; }

  private:

    /// setup fluid algorithm (overriding some fluid parameters with
    /// values specified in given problem-dependent ParameterList)
    void SetupFluid(const Teuchos::ParameterList& prbdyn, bool& isale);

    /// fluid field solver
    Teuchos::RCP<Fluid> fluid_;
  };

}

#endif
#endif
