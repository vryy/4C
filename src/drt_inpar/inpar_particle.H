/*---------------------------------------------------------------------------*/
/*!
\file inpar_particle.H

\brief input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef INPAR_PARTICLE_H
#define INPAR_PARTICLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | input parameters for particle problems                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace INPAR
{
  namespace PARTICLE
  {
    /*---------------------------------------------------------------------------*
     | general control parameters for particle simulations                       |
     *---------------------------------------------------------------------------*/

    //! type of particle time integration
    enum DynamicType
    {
      dyna_semiimpliciteuler,  //! semi implicit euler scheme (explicit)
      dyna_velocityverlet      //! velocity verlet scheme (explicit)
    };

    //! type of particle interaction
    enum InteractionType
    {
      interaction_none,  //! no particle interaction
    };

    //! data format for written numeric data via vtp
    enum OutputDataFormat
    {
      binary,
      ascii
    };


    //! set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PARTICLE

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
#endif
