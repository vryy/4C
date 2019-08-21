/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class defining the fluid boundary element

\level 1

\maintainer  Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_BOUNDARY_INTERFACE_H
#define FLUID_ELE_BOUNDARY_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
  class Discretization;
  class Condition;

  namespace ELEMENTS
  {
    class FluidBoundary;

    /// Interface base class for FluidEleBoundaryCalc
    /*!
      This class exists to provide a common interface for all template
      versions of FluidBoundaryCalc.
     */
    //  class FluidBoundaryImplInterface : DRT::SingletonDestruction
    class FluidBoundaryInterface : DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      FluidBoundaryInterface(){};

      /// Empty destructor
      virtual ~FluidBoundaryInterface(){};

      virtual void EvaluateAction(DRT::ELEMENTS::FluidBoundary* ele1,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3) = 0;

      virtual int EvaluateNeumann(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseMatrix* elemat1) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
