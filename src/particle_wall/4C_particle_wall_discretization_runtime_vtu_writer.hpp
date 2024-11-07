// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP
#define FOUR_C_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::IO
{
  class DiscretizationVisualizationWriterMesh;
  class DiscretizationReader;
}  // namespace Core::IO

namespace PARTICLEWALL
{
  class WallDataState;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  /*!
   * \brief particle wall discretization runtime vtu writer class
   *
   * A class that writes visualization output for particle wall discretization in vtk/vtp format at
   * runtime.
   *
   * \author Sebastian Fuchs \date 08/2019
   */
  class WallDiscretizationRuntimeVtuWriter final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] walldiscretization wall discretization
     * \param[in] walldatastate      wall data state container
     * \param[in] restart_time       restart time of the simulation
     *
     */
    explicit WallDiscretizationRuntimeVtuWriter(
        const std::shared_ptr<Core::FE::Discretization> walldiscretization,
        const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate, double restart_time);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~WallDiscretizationRuntimeVtuWriter() = default;

    /*!
     * \brief write wall discretization runtime output
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] step output step
     * \param[in] time output time
     */
    void write_wall_discretization_runtime_output(const int step, const double time) const;

   private:
    //! wall discretization
    std::shared_ptr<Core::FE::Discretization> walldiscretization_;

    //! wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate_;

    //! vtu writer object
    std::unique_ptr<Core::IO::DiscretizationVisualizationWriterMesh> runtime_vtuwriter_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
