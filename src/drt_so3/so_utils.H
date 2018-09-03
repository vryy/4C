/*!----------------------------------------------------------------------
\file so_utils.H

\brief A collection of helper methods for solid elements

<pre>
\level 1
\maintainer Martin Pfaller
</pre>
*-----------------------------------------------------------------------*/

#ifndef SO_UTILS_H
#define SO_UTILS_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace DRT
{
  namespace ELEMENTS
  {
    namespace UTILS
    {
      template <DRT::Element::DiscretizationType distype>
      void CalcR(const DRT::Element* ele, const std::vector<double>& disp,
          LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim,
              DRT::UTILS::DisTypeToDim<distype>::dim>& R);
    }
  }  // namespace ELEMENTS
}  // namespace DRT

#endif
