/*-----------------------------------------------------------------------------------------------*/
/*!

\brief base class for a triad interpolation scheme

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "triad_interpolation.H"

#include "triad_interpolation_local_rotation_vectors.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>
#include <Sacado.hpp>

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
LARGEROTATIONS::TriadInterpolation<T>::TriadInterpolation()
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
Teuchos::RCP<LARGEROTATIONS::TriadInterpolation<T>> LARGEROTATIONS::TriadInterpolation<T>::Create(
    unsigned int numnodes)
{
  // so far, the only implemented variant is the one based on local rotation vectors

  switch (numnodes)
  {
    case 2:
    {
      return Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<2, T>());
    }
    case 3:
    {
      return Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, T>());
    }
    case 4:
    {
      return Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<4, T>());
    }
    case 5:
    {
      return Teuchos::rcp(new LARGEROTATIONS::TriadInterpolationLocalRotationVectors<5, T>());
    }
    default:
    {
      dserror("%d is no valid number of nodes used for triad interpolation! choose 2,3,4 or 5",
          numnodes);
      break;
    }
  }

  return Teuchos::null;
}

// explicit template instantiations
template class LARGEROTATIONS::TriadInterpolation<double>;
template class LARGEROTATIONS::TriadInterpolation<Sacado::Fad::DFad<double>>;
