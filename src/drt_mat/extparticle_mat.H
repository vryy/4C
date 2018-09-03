/*----------------------------------------------------------------------*/
/*!
\file extparticle_mat.H

\brief extended particle material

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/
#ifndef ExtParticleMat_H
#define ExtParticleMat_H

#include "particle_mat.H"

namespace MAT
{
  namespace PAR
  {
    //! parameters for ExtParticleMat material
    class ExtParticleMat : public ParticleMat
    {
     public:
      //! constructor
      ExtParticleMat(Teuchos::RCP<MAT::PAR::Material> matdata);
      //! destructor
      virtual ~ExtParticleMat() { return; };

      //! return the nominal speed of sound in case of small density variations and liquid phase
      virtual double SpeedOfSoundL() const { return std::pow(bulkModulus_ / initDensity_, 0.5); }

      //! create instance of particle additive manufacturing material
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// factor defining the equation of state: pressure=c_0^2(density-refdensfac_*initDensity_)
      /// with the speed of sound c_0^2=BULK_MODULUS/initDensity_;
      const double refdensfac_;
      /// exponent occurring in equation of state. For exponents !=1, we have:
      /// pressure=p0*((density/initDensity_)^exponent_-refdensfac_) with
      /// p0=c_0^2*Density_/exponent_;
      const double exponent_;
      /// bulk modulus. 1/compressibility for the liquid phase
      const double bulkModulus_;
      /// dynamic (shear) viscosity for the liquid phase
      const double dynamicViscosity_;
      /// bulk viscosity for the liquid phase
      const double bulkViscosity_;
      /// artificial viscosity for the liquid phase
      const double artificialViscosity_;

    };  // class MAT::PAR::ExtParticleMat
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ExtParticleMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ExtParticleMatType"; };

    static ExtParticleMatType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ExtParticleMatType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for ExtParticleMat material
  class ExtParticleMat : public ParticleMat
  {
   public:
    //! construct empty ExtParticleMat material
    ExtParticleMat();

    //! construct ExtParticleMat material with specific material parameters
    explicit ExtParticleMat(MAT::PAR::ExtParticleMat* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const { return ExtParticleMatType::Instance().UniqueParObjectId(); };

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void Unpack(const std::vector<char>& data);
    //@}

    //! return material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_extparticlemat; };

    //! clone ExtParticleMat material
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ExtParticleMat(*this)); };

   private:
    //! return material parameters
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! material parameters
    MAT::PAR::ExtParticleMat* params_;
  };
}  // namespace MAT
#endif
