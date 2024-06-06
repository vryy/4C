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
#ifndef FOUR_C_SCATRA_ELE_PARAMETER_BASE_HPP
#define FOUR_C_SCATRA_ELE_PARAMETER_BASE_HPP


#include "4C_config.hpp"

#include <memory>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

namespace Discret
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
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
