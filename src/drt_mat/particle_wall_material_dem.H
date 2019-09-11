/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall material for DEM

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_MATERIAL_BASE_H
#define PARTICLE_WALL_MATERIAL_BASE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    class ParticleWallMaterialDEM : public Parameter
    {
     public:
      //! constructor
      ParticleWallMaterialDEM(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! desctructor
      ~ParticleWallMaterialDEM() override{};

      //! create material instance of matching type with parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! @name material parameters
      //@{

      //! tangential contact friction coefficient
      const double frictionTang_;

      //! rolling contact friction coefficient
      const double frictionRoll_;

      //! adhesion surface energy
      const double adhesionSurfaceEnergy_;

      //@}
    };

  }  // namespace PAR

  class ParticleWallMaterialDEMType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ParticleWallMaterialDEMType"; };

    static ParticleWallMaterialDEMType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ParticleWallMaterialDEMType instance_;
  };

  class ParticleWallMaterialDEM : public Material
  {
   public:
    //! constructor (empty material object)
    ParticleWallMaterialDEM();

    //! constructor (with given material parameters)
    explicit ParticleWallMaterialDEM(MAT::PAR::ParticleWallMaterialDEM* params);

    //! @name Packing and Unpacking

    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ParticleWallMaterialDEMType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const override;

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

    //! @name Access methods

    //@{

    //! material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_particle_wall_dem;
    }

    //! return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ParticleWallMaterialDEM(*this));
    }

    //! return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    //! return tangential contact friction coefficient
    double MuTangential() const { return params_->frictionTang_; }

    //! return rolling contact friction coefficient
    double MuRolling() const { return params_->frictionRoll_; }

    //! return adhesion surface energy
    double AdhesionSurfaceEnergy() const { return params_->adhesionSurfaceEnergy_; }

    //@}

   private:
    //! my material parameters
    MAT::PAR::ParticleWallMaterialDEM* params_;
  };

}  // namespace MAT

/*---------------------------------------------------------------------------*/
#endif
