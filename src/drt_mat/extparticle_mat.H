/*!----------------------------------------------------------------------*/
/*!
\file extparticle_mat.H

\brief Multi-purpose Particle-MeshFree material with support for thermodynamics

\level 3

\maintainer Alessandro Cattabiani
*/
/*----------------------------------------------------------------------*/
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
      ExtParticleMat(
            Teuchos::RCP<MAT::PAR::Material> matdata
            );
      //! destructor
      virtual ~ExtParticleMat(){return;};

      //! return the threshold specific enthalpy (solid <-> transition)
      virtual double SpecEnthalpyST() const { return transitionTemperature_ * CPS_; }

      //! return the threshold specific enthalpy (transition <-> liquid)
      virtual double SpecEnthalpyTL() const { return SpecEnthalpyST() + latentHeat_; }

      //! return the nominal speed of sound in case of small density variations and solid phase
      virtual double SpeedOfSoundS() const { return std::pow(youngModulus_/initDensity_,0.5); }

      //! return the nominal speed of sound in case of small density variations and liquid phase
      virtual double SpeedOfSoundL() const { return std::pow(bulkModulus_/initDensity_,0.5); }

      //! return the nominal speed of sound in case of small density variations and transition phase
      virtual double SpeedOfSoundT(const double specEnthalpy) const
      {
        // safety check
        if (specEnthalpy < SpecEnthalpyST() || specEnthalpy > SpecEnthalpyTL())
          dserror("This particle is not in transition phase");

        // linear interpolation between the Young and compressibility moduli
        const double alpha = ( bulkModulus_ - youngModulus_)/latentHeat_;
        const double qq = youngModulus_ - SpecEnthalpyST() * alpha;
        const double xx = alpha * specEnthalpy + qq;

        return std::pow(xx/initDensity_,0.5);
      }

      //! create instance of particle additive manufacturing material
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// dismembering radius
      const double dismemberRadius_;
      /// initial temperature
      const double initTemperature_;
      /// specific heat - constant pressure - solid state
      const double CPS_;
      /// specific heat - constant pressure - liquid state
      const double CPL_;
      /// specific latent heat - solid <-> liquid
      const double latentHeat_;
      /// transition temperature - solid <-> liquid
      const double transitionTemperature_;
      /// volumetric thermal expansion coefficient - solid state - usually 3 times the linear coefficient
      const double thermalExpansionS_;
      /// volumetric thermal expansion coefficient - liquid state
      const double thermalExpansionL_;
      /// volumetric thermal expansion coefficient - solid <-> liquid - Beware! It is linked to specEnthalpy instead of temperature!
      const double thermalExpansionT_;
      /// bulk modulus. 1/compressibility for the liquid phase
      const double bulkModulus_;
      /// dynamic (shear) viscosity for the liquid phase
      const double dynamicViscosity_;
      /// thermal conductivity of the material (it is considered the same for all the phases) - W/(m * K)
      const double thermalConductivity_;


    };  // class MAT::PAR::ExtParticleMat
  } // namespace PAR


  /*----------------------------------------------------------------------*/
  class ExtParticleMatType : public DRT::ParObjectType
  {
    public:

      std::string Name() const {return "ExtParticleMatType";};

      static ExtParticleMatType& Instance() {return instance_;};

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
      int UniqueParObjectId() const {return ExtParticleMatType::Instance().UniqueParObjectId();};

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
      INPAR::MAT::MaterialType MaterialType() const {return INPAR::MAT::m_extparticlemat;};

      //! clone ExtParticleMat material
      Teuchos::RCP<Material> Clone() const {return Teuchos::rcp(new ExtParticleMat(*this));};

    private:

      //! return material parameters
      virtual MAT::PAR::Parameter* Parameter() const {return params_;}

      //! material parameters
      MAT::PAR::ExtParticleMat* params_;
  };
}
#endif
