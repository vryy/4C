/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP
#define FOUR_C_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include <Teuchos_RCP.hpp>

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationVisualizationWriterMesh;
  class DiscretizationReader;
}  // namespace IO

namespace PARTICLEWALL
{
  class WallDataState;
}

namespace DRT
{
  class Discretization;
}

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
        const Teuchos::RCP<DRT::Discretization> walldiscretization,
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
    void WriteWallDiscretizationRuntimeOutput(const int step, const double time) const;

   private:
    //! wall discretization
    Teuchos::RCP<DRT::Discretization> walldiscretization_;

    //! wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate_;

    //! vtu writer object
    std::unique_ptr<IO::DiscretizationVisualizationWriterMesh> runtime_vtuwriter_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
