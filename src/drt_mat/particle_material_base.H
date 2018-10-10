/*---------------------------------------------------------------------------*/
/*!
\file particle_material_base.H

\brief particle material base

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_MATERIAL_BASE_H
#define PARTICLE_MATERIAL_BASE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialBase : virtual public Parameter
    {
     public:
      //! constructor
      ParticleMaterialBase(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! virtual desctructor
      virtual ~ParticleMaterialBase(){};

      //! @name material parameters
      //@{

      //! initial radius
      const double initRadius_;

      //! initial density
      const double initDensity_;

      //@}

      //! create material instance of matching type with parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() = 0;
    };

  }  // namespace PAR

}  // namespace MAT

/*---------------------------------------------------------------------------*/
#endif
