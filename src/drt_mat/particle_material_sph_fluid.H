/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH fluid

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_MATERIAL_SPH_H
#define PARTICLE_MATERIAL_SPH_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_material_base.H"
#include "particle_material_thermo.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialSPHFluid : public ParticleMaterialBase, public ParticleMaterialThermo
    {
     public:
      //! constructor
      ParticleMaterialSPHFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! virtual desctructor
      virtual ~ParticleMaterialSPHFluid(){};

      //! speed of sound
      virtual double SpeedOfSound() const { return std::sqrt(bulkModulus_ / initDensity_); };

      //! @name material parameters
      //@{

      //! reference density factor in equation of state
      const double refDensFac_;

      //! exponent in equation of state
      const double exponent_;

      //! background pressure
      const double backgroundPressure_;

      //! bulk modulus
      const double bulkModulus_;

      //! dynamic shear viscosity
      const double dynamicViscosity_;

      //! bulk viscosity
      const double bulkViscosity_;

      //! artificial viscosity
      const double artificialViscosity_;

      //@}

      //! create material instance of matching type with parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
    };

  }  // namespace PAR

  class ParticleMaterialSPHFluidType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleMaterialSPHType"; };

    static ParticleMaterialSPHFluidType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleMaterialSPHFluidType instance_;
  };

  class ParticleMaterialSPHFluid : public Material
  {
   public:
    //! constructor (empty material object)
    ParticleMaterialSPHFluid();

    //! constructor (with given material parameters)
    explicit ParticleMaterialSPHFluid(MAT::PAR::ParticleMaterialSPHFluid* params);

    //! @name Packing and Unpacking

    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return ParticleMaterialSPHFluidType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_particle_sph_fluid;
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ParticleMaterialSPHFluid(*this));
    }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    //! my material parameters
    MAT::PAR::ParticleMaterialSPHFluid* params_;
  };

}  // namespace MAT

/*---------------------------------------------------------------------------*/
#endif
