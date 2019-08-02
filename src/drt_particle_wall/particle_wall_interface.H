/*---------------------------------------------------------------------------*/
/*!
\brief pure virtual interface to provide access to particle wall handler

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_INTERFACE_H
#define PARTICLE_WALL_INTERFACE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}

namespace PARTICLEWALL
{
  class WallDataState;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  class WallHandlerInterface
  {
   public:
    //! virtual destructor
    virtual ~WallHandlerInterface(){};

    //! get wall discretization
    virtual Teuchos::RCP<const DRT::Discretization> GetWallDiscretization() const = 0;

    //! get wall data state container
    virtual std::shared_ptr<PARTICLEWALL::WallDataState> GetWallDataState() const = 0;

    //! get reference to potential wall neighbors
    virtual const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors() const = 0;

    //! determine nodal positions of column wall element
    virtual void DetermineColWallEleNodalPos(
        DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) const = 0;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
