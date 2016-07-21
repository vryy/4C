/*!----------------------------------------------------------------------*/
/*!
\file particleAMmat.H

\brief Particle material with support for thermodynamics

\level 3

\maintainer Alessandro Cattabiani
*/
/*----------------------------------------------------------------------*/
#ifndef PARTICLEAMMAT_H
#define PARTICLEAMMAT_H

#include "particle_mat.H"

namespace MAT
{
  namespace PAR
  {
    //! parameters for ParticleAMmat material
    class ParticleAMmat : public ParticleMat
    {
      public:

        //! constructor
      ParticleAMmat(
            Teuchos::RCP<MAT::PAR::Material> matdata
            );

        //! destructor
        virtual ~ParticleAMmat(){return;};

        //! create instance of particle additive manufacturing material
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        /// initial temperature
        const double temperature_;
        /// specific heat - constant pressure - solid state
        const double CPS_;
        /// specific heat - constant pressure - liquid state
        const double CPL_;
        /// specific latent heat - solid <-> liquid
        const double SL_latent_heat_max_;
        /// transition temperature - solid <-> liquid
        const double SL_transitionTemperature_;
        /// volumetric thermal expansion coefficient - solid state - usually 3 times the linear coefficient
        const double S_thermalExpansion_;
        /// volumetric thermal expansion coefficient - liquid state
        const double L_thermalExpansion_;
        /// volumetric thermal expansion coefficient - solid <-> liquid - Beware! It is linked to latent heat instead of temperature!
        const double SL_thermalExpansion_;

    };  // class MAT::PAR::ParticleAMmat
  } // namespace PAR


  /*----------------------------------------------------------------------*/
  class ParticleAMmatType : public DRT::ParObjectType
  {
    public:

      std::string Name() const {return "ParticleAMmatType";};

      static ParticleAMmatType& Instance() {return instance_;};

      DRT::ParObject* Create(const std::vector<char>& data);

    private:

      static ParticleAMmatType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for ParticleAMmat material
  class ParticleAMmat : public ParticleMat
  {
    public:

      //! construct empty ParticleAMmat material
    ParticleAMmat();

      //! construct ParticleAMmat material with specific material parameters
      explicit ParticleAMmat(MAT::PAR::ParticleAMmat* params);

      //! @name packing and unpacking
      /*!
        \brief Return unique ParObject id

        Every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
      */
      int UniqueParObjectId() const {return ParticleAMmatType::Instance().UniqueParObjectId();};

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
      INPAR::MAT::MaterialType MaterialType() const {return INPAR::MAT::m_particleAMmat;};

      //! clone ParticleAMmat material
      Teuchos::RCP<Material> Clone() const {return Teuchos::rcp(new ParticleAMmat(*this));};

    private:

      //! return material parameters
      virtual MAT::PAR::Parameter* Parameter() const {return params_;}

      //! material parameters
      MAT::PAR::ParticleAMmat* params_;
  };
}
#endif
