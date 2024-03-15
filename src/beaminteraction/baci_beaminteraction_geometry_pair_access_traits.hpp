/*----------------------------------------------------------------------*/
/*! \file

\brief Specialize structures to correctly initialize the element data containers for Hermite
elements

\level 1

*/


#ifndef BACI_BEAMINTERACTION_GEOMETRY_PAIR_ACCESS_TRAITS_HPP
#define BACI_BEAMINTERACTION_GEOMETRY_PAIR_ACCESS_TRAITS_HPP

#include "baci_config.hpp"

#include "baci_beam3_base.hpp"
#include "baci_geometry_pair_element.hpp"

BACI_NAMESPACE_OPEN


namespace GEOMETRYPAIR
{
  template <>
  struct SetShapeFunctionData<t_hermite>
  {
    static void Set(ShapeFunctionData<t_hermite>& shape_function_data, const DRT::Element* element)
    {
      // Get the reference length of the beam element
      const auto* beam_element = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element);
      if (beam_element == nullptr)
        dserror(
            "The element pointer has to point to a valid beam element when evaluating the shape "
            "function data of a hermite beam, as we need to get RefLength()!");
      shape_function_data.ref_length_ = beam_element->RefLength();
    }
  };
}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE

#endif
