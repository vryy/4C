/*----------------------------------------------------------------------*/
/*! \file
\brief material properties of particles

\level 2

\maintainer  Christoph Meier

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_MAT_H
#define PARTICLE_MAT_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for particles
    class ParticleMat : public Parameter
    {
     public:
      /// standard constructor
      ParticleMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ParticleMat() { ; }

      /// @name material parameters
      //@{

      /// density
      const double initDensity_;
      /// initial radius of particle
      const double initRadius_;
      /// poisson ratio
      const double poissonRatio_;
      /// youngs modulus
      const double youngModulus_;
      /// yield strength
      const double yieldStrength_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ParticleMat

  }  // namespace PAR

  class ParticleMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleMatType"; }

    static ParticleMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// particle material
  class ParticleMat : public Material
  {
   public:
    /// construct empty material object
    ParticleMat();

    /// construct the material object given material parameters
    explicit ParticleMat(MAT::PAR::ParticleMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return ParticleMatType::Instance().UniqueParObjectId();
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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_particlemat; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ParticleMat(*this)); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ParticleMat* params_;
  };

}  // namespace MAT

#endif
