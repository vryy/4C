/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH boundary

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_PARTICLE_SPH_BOUNDARY_HPP
#define FOUR_C_MAT_PARTICLE_SPH_BOUNDARY_HPP

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_mat_particle_base.hpp"
#include "4C_mat_particle_thermo.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialSPHBoundary : public ParticleMaterialBase, public ParticleMaterialThermo
    {
     public:
      //! constructor
      ParticleMaterialSPHBoundary(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      //! create material instance of matching type with parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;
    };

  }  // namespace PAR

  class ParticleMaterialSPHBoundaryType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ParticleMaterialSPHBoundaryType"; };

    static ParticleMaterialSPHBoundaryType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ParticleMaterialSPHBoundaryType instance_;
  };

  class ParticleMaterialSPHBoundary : public CORE::MAT::Material
  {
   public:
    //! constructor (empty material object)
    ParticleMaterialSPHBoundary();

    //! constructor (with given material parameters)
    explicit ParticleMaterialSPHBoundary(MAT::PAR::ParticleMaterialSPHBoundary* params);

    //! @name Packing and Unpacking

    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ParticleMaterialSPHBoundaryType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_particle_sph_boundary;
    }

    //! return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new ParticleMaterialSPHBoundary(*this));
    }

    //! return quick accessible material parameter data
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    //! my material parameters
    MAT::PAR::ParticleMaterialSPHBoundary* params_;
  };

}  // namespace MAT

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
