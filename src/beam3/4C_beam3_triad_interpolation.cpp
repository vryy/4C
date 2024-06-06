/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief base class for a triad interpolation scheme


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_beam3_triad_interpolation.hpp"

#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_utils_exceptions.hpp"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
LargeRotations::TriadInterpolation<T>::TriadInterpolation()
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
Teuchos::RCP<LargeRotations::TriadInterpolation<T>> LargeRotations::TriadInterpolation<T>::Create(
    unsigned int numnodes)
{
  // so far, the only implemented variant is the one based on local rotation vectors

  switch (numnodes)
  {
    case 2:
    {
      return Teuchos::rcp(new LargeRotations::TriadInterpolationLocalRotationVectors<2, T>());
    }
    case 3:
    {
      return Teuchos::rcp(new LargeRotations::TriadInterpolationLocalRotationVectors<3, T>());
    }
    case 4:
    {
      return Teuchos::rcp(new LargeRotations::TriadInterpolationLocalRotationVectors<4, T>());
    }
    case 5:
    {
      return Teuchos::rcp(new LargeRotations::TriadInterpolationLocalRotationVectors<5, T>());
    }
    default:
    {
      FOUR_C_THROW("%d is no valid number of nodes used for triad interpolation! choose 2,3,4 or 5",
          numnodes);
      break;
    }
  }

  return Teuchos::null;
}

// explicit template instantiations
template class LargeRotations::TriadInterpolation<double>;
template class LargeRotations::TriadInterpolation<Sacado::Fad::DFad<double>>;

FOUR_C_NAMESPACE_CLOSE
