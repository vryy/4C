/*---------------------------------------------------------------------------*/
/*! \file
\brief gravity acceleration handler for particle simulations

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_GRAVITY_H
#define PARTICLE_GRAVITY_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief gravity acceleration handler for particle simulations
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class GravityHandler
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] params particle simulation parameter list
     */
    explicit GravityHandler(const Teuchos::ParameterList& params);

    //! destructor
    ~GravityHandler(){};

    /*!
     * \brief init gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] gravity gravity acceleration
     */
    void Init(const std::vector<double>& gravity);

    /*!
     * \brief setup gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void Setup();

    /*!
     * \brief write restart of gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief get gravity acceleration
     *
     * Evaluate the gravity ramp function at the given time to get the scaled gravity acceleration.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in]  time           evaluation time
     * \param[out] scaled_gravity scaled gravity acceleration
     */
    void GetGravityAcceleration(const double time, std::vector<double>& scaled_gravity);

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! gravity acceleration vector
    std::vector<double> gravity_;

    //! gravity ramp function number
    const int gravityrampfctnumber_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
