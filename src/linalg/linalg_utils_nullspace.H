/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0

*/
/*---------------------------------------------------------------------*/
#ifndef LINALG_UTILS_NULLSPACE_H
#define LINALG_UTILS_NULLSPACE_H

#include <vector>

#include "../drt_lib/drt_discret.H"

namespace LINALG
{
  // forward declarations
  class Discretization;

  void ComputeStructure3DNullSpace(
      DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
  void ComputeStructure2DNullSpace(
      DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
  void ComputeShell3DNullSpace(
      DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
  void ComputeXFluidDNullSpace(
      DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
  void ComputeFluidDNullSpace(
      DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns);
}  // namespace LINALG

#endif
