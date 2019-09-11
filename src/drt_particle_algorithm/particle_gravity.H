/*---------------------------------------------------------------------------*/
/*! \file
\brief gravity acceleration handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_GRAVITY_H
#define PARTICLE_GRAVITY_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class GravityHandler
  {
   public:
    //! constructor
    explicit GravityHandler(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~GravityHandler(){};

    //! init gravity handler
    void Init(const std::vector<double>& gravity);

    //! setup gravity handler
    void Setup();

    //! write restart of gravity handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of gravity handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get gravity acceleration
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
