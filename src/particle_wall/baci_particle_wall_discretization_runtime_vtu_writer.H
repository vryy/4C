/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define BACI_PARTICLE_WALL_DISCRETIZATION_RUNTIME_VTU_WRITER_H

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

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
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
    std::unique_ptr<DiscretizationRuntimeVtuWriter> runtime_vtuwriter_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
