/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
class DiscretizationRuntimeVtuWriter;

namespace IO
{
  class DiscretizationReader;
}

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
     */
    explicit WallDiscretizationRuntimeVtuWriter();

    //! destructor
    ~WallDiscretizationRuntimeVtuWriter();

    /*!
     * \brief init wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] walldiscretization wall discretization
     * \param[in] walldatastate      wall data state container
     */
    void Init(const Teuchos::RCP<DRT::Discretization> walldiscretization,
        const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate);

    /*!
     * \brief setup wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] write_binary_output     flag for output in binary format
     */
    void Setup(bool write_binary_output);

    /*!
     * \brief write restart of wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

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
    //! setup time of runtime vtu writer
    double setuptime_;

    //! wall discretization
    Teuchos::RCP<DRT::Discretization> walldiscretization_;

    //! wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate_;

    //! vtu writer object
    std::unique_ptr<DiscretizationRuntimeVtuWriter> runtime_vtuwriter_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
