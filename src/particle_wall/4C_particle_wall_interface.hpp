/*---------------------------------------------------------------------------*/
/*! \file
\brief interface to provide restricted access to particle wall handler
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_WALL_INTERFACE_HPP
#define FOUR_C_PARTICLE_WALL_INTERFACE_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_particle_engine_typedefs.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Discret
{
  class Discretization;
}

namespace PARTICLEWALL
{
  class WallDataState;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  /*!
   * \brief interface to provide restricted access to particle wall handler
   *
   * The particle algorithm holds an instance of the particle wall handler, thus having full access
   * and control over it. This abstract interface class to the particle wall handler provides
   * restricted access to be used in all other classes.
   *
   * \note Methods in this class are documented briefly. Refer to the full documentation of the
   *       particle wall handler class!
   *
   * \author Sebastian Fuchs \date 10/2018
   */
  class WallHandlerInterface
  {
   public:
    //! virtual destructor
    virtual ~WallHandlerInterface() = default;

    /*!
     * \brief get wall discretization
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return wall discretization
     */
    virtual Teuchos::RCP<const Discret::Discretization> get_wall_discretization() const = 0;

    /*!
     * \brief get wall data state container
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return wall data state container
     */
    virtual std::shared_ptr<PARTICLEWALL::WallDataState> GetWallDataState() const = 0;

    /*!
     * \brief get reference to potential wall neighbors
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return potential particle wall neighbor pairs
     */
    virtual const PARTICLEENGINE::PotentialWallNeighbors& get_potential_wall_neighbors() const = 0;

    /*!
     * \brief determine nodal positions of column wall element
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param ele[in]             column wall element
     * \param colelenodalpos[out] current nodal position
     */
    virtual void determine_col_wall_ele_nodal_pos(Core::Elements::Element* ele,
        std::map<int, Core::LinAlg::Matrix<3, 1>>& colelenodalpos) const = 0;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
