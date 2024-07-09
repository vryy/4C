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
namespace Mat
{
  namespace PAR
  {
    class ParticleMaterialSPHBoundary : public ParticleMaterialBase, public ParticleMaterialThermo
    {
     public:
      //! constructor
      ParticleMaterialSPHBoundary(const Core::Mat::PAR::Parameter::Data& matdata);

      //! create material instance of matching type with parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };

  }  // namespace PAR

  class ParticleMaterialSPHBoundaryType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "ParticleMaterialSPHBoundaryType"; };

    static ParticleMaterialSPHBoundaryType& instance() { return instance_; };

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    static ParticleMaterialSPHBoundaryType instance_;
  };

  class ParticleMaterialSPHBoundary : public Core::Mat::Material
  {
   public:
    //! constructor (empty material object)
    ParticleMaterialSPHBoundary();

    //! constructor (with given material parameters)
    explicit ParticleMaterialSPHBoundary(Mat::PAR::ParticleMaterialSPHBoundary* params);

    //! @name Packing and Unpacking

    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int unique_par_object_id() const override
    {
      return ParticleMaterialSPHBoundaryType::instance().unique_par_object_id();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by unique_par_object_id() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      unique_par_object_id().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void unpack(const std::vector<char>& data) override;

    //@}

    //! material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_particle_sph_boundary;
    }

    //! return copy of this material object
    Teuchos::RCP<Core::Mat::Material> clone() const override
    {
      return Teuchos::rcp(new ParticleMaterialSPHBoundary(*this));
    }

    //! return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

   private:
    //! my material parameters
    Mat::PAR::ParticleMaterialSPHBoundary* params_;
  };

}  // namespace Mat

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
