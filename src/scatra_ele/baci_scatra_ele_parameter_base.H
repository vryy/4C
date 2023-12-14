/*----------------------------------------------------------------------*/
/*! \file

\brief abstract interface for parameter classes required for scalar transport element evaluation

This abstract, purely virtual singleton class provides a common inheritance interface for several
derived singleton classes holding static parameters required for scalar transport element
evaluation, e.g., general parameters, turbulence parameters, and problem specific parameters. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation.


\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_SCATRA_ELE_PARAMETER_BASE_H
#define BACI_SCATRA_ELE_PARAMETER_BASE_H


#include <memory>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterBase
    {
     public:
      /// Virtual destructor.
      virtual ~ScaTraEleParameterBase() = default;

      //! set parameters
      virtual void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
