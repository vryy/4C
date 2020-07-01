/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/
#ifndef UTILS_NULLSPACE_H
#define UTILS_NULLSPACE_H

#include <vector>

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace UTILS
  {
    void ComputeStructure3DNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
    void ComputeStructure2DNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
    void ComputeBeam2DNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
    void ComputeBeam3DNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
    void ComputeXFluidDNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
    void ComputeFluidDNullSpace(
        DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);

  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_NULLSPACE_H
