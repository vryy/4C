/*!----------------------------------------------------------------------*/
/*!
\file particle_mat_ellipsoids.H

\brief material properties of ellipsoidal particles

\level 2

<pre>
\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15262
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PARTICLE_MAT_ELLIPSOIDS_H
#define PARTICLE_MAT_ELLIPSOIDS_H

#include "particle_mat.H"

namespace MAT
{
  namespace PAR
  {
    //! material parameters
    class ParticleMatEllipsoids : public ParticleMat
    {
     public:
      //! constructor
      ParticleMatEllipsoids(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~ParticleMatEllipsoids() { return; };

      //! create instance of particle material
      Teuchos::RCP<MAT::Material> CreateMaterial();

      //! return semi-axes of ellipsoidal particles
      const LINALG::Matrix<3, 1>& SemiAxes() const { return semiaxes_; };

     private:
      //! semi-axes of ellipsoidal particles
      LINALG::Matrix<3, 1> semiaxes_;
    };  // class ParticleMatEllipsoids
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  class ParticleMatEllipsoidsType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleMatEllipsoidsType"; };

    static ParticleMatEllipsoidsType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleMatEllipsoidsType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for particle material
  class ParticleMatEllipsoids : public ParticleMat
  {
   public:
    //! construct empty particle material
    ParticleMatEllipsoids();

    //! construct particle material with specific material parameters
    explicit ParticleMatEllipsoids(MAT::PAR::ParticleMatEllipsoids* params);

    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const
    {
      return ParticleMatEllipsoidsType::Instance().UniqueParObjectId();
    };

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
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_particlemat_ellipsoids; };

    //! clone particle material
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ParticleMatEllipsoids(*this)); };

    //! return semi-axes of ellipsoidal particles
    const LINALG::Matrix<3, 1>& SemiAxes() const { return params_->SemiAxes(); };

   private:
    //! material parameters
    MAT::PAR::ParticleMatEllipsoids* params_;
  };  // wrapper for particle material
}  // namespace MAT
#endif
