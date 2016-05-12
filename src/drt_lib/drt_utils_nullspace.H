/*!----------------------------------------------------------------------
\file drt_utils_nullspace.H
\brief Declaration A collection of helper methods for namespace DRT
<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef UTILS_NULLSPACE_H
#define UTILS_NULLSPACE_H

#include <vector>

namespace DRT
{
// forward declarations
class Discretization;

/// collection of useful and important but less central stuff
namespace UTILS
{


  void ComputeStructure3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );
  void ComputeStructure2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );
  void ComputeBeam2DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );
  void ComputeBeam3DNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );
  void ComputeXFluidDNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );
  void ComputeFluidDNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_NULLSPACE_H




