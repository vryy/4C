#ifndef FLUID_ELE_FACTORY_H
#define FLUID_ELE_FACTORY_H

#include "fluid_ele_interface.H"

namespace DRT
{

namespace ELEMENTS
{
class FluidFactory
{
public:
  FluidFactory();

  virtual ~FluidFactory(){}

  static FluidEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype, string problem);

private:
  template<DRT::Element::DiscretizationType distype>
  static FluidEleInterface* DefineProblemType(string problem);

};

}
}

#endif
