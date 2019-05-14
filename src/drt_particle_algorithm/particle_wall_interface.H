/*---------------------------------------------------------------------------*/
/*!
\file particle_wall_interface.H

\brief pure virtual interface to provide access to particle wall handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class WallHandlerInterface
  {
   public:
    //! virtual destructor
    virtual ~WallHandlerInterface(){};

    //! get wall discretization
    virtual Teuchos::RCP<const DRT::Discretization> GetWallDiscretization() const = 0;

    //! get wall displacements
    virtual Teuchos::RCP<const Epetra_Vector> GetDisp() const = 0;

    //! get wall velocities
    virtual Teuchos::RCP<const Epetra_Vector> GetVel() const = 0;

    //! get reference to potential wall neighbors
    virtual const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors() const = 0;

    //! determine nodal positions of column wall element
    virtual void DetermineColWallEleNodalPos(
        DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) const = 0;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
