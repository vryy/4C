/*----------------------------------------------------------------------*/
/*!

\brief abstract interface for parameter classes required for scalar transport element evaluation

This abstract, purely virtual singleton class provides a common inheritance interface for several
derived singleton classes holding static parameters required for scalar transport element
evaluation, e.g., general parameters, turbulence parameters, and problem specific parameters. All
parameters are usually set only once at the beginning of a simulation, namely during initialization
of the global time integrator, and then never touched again throughout the simulation.

\maintainer Anh-Tu Vuong

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_PARAMETER_BASE_H
#define SCATRA_ELE_PARAMETER_BASE_H

#include "../drt_lib/drt_singletondestruction.H"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleParameterBase : public DRT::SingletonDestruction
    {
     public:
      //! singleton destruction
      virtual void Done() = 0;

      //! set parameters
      virtual void SetParameters(Teuchos::ParameterList& parameters  //!< parameter list
          ) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
