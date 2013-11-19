/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_parameter_timint.cpp

\brief Common interface for fluid element parameter

<pre>
Maintainers: Ursula Rasthofer & Volker Gravemeier
             {rasthofer,vgravem}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef FLUID_ELE_PARAMETER_INTERFACE_H
#define FLUID_ELE_PARAMETER_INTERFACE_H

#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
namespace ELEMENTS
{

/// Interface base class for FluidEleParameter

class FluidEleParameterInterface : public DRT::SingletonDestruction
{
public:

  //! empty constructor
  FluidEleParameterInterface() {}

  //! empty destructor
  virtual ~FluidEleParameterInterface() {}

};

}
}

#endif
