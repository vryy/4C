/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2019 |
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
 | class definitions                                          sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  class WallDiscretizationRuntimeVtuWriter final
  {
   public:
    //! constructor
    explicit WallDiscretizationRuntimeVtuWriter();

    //! destructor
    ~WallDiscretizationRuntimeVtuWriter();

    //! init wall discretization runtime vtu writer
    void Init(const Teuchos::RCP<DRT::Discretization> walldiscretization,
        const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate);

    //! setup wall discretization runtime vtu writer
    void Setup(bool write_binary_output);

    //! write restart of wall discretization runtime vtu writer
    void WriteRestart(const int step, const double time) const;

    //! read restart of wall discretization runtime vtu writer
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! write wall discretization runtime output
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
