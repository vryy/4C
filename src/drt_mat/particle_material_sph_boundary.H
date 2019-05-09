/*---------------------------------------------------------------------------*/
/*!
\brief particle material for SPH boundary

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_MATERIAL_SPH_BOUNDARY_H
#define PARTICLE_MATERIAL_SPH_BOUNDARY_H

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
    class ParticleMaterialSPHBoundary : public ParticleMaterialBase, public ParticleMaterialThermo
    {
     public:
      //! constructor
      ParticleMaterialSPHBoundary(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! virtual desctructor
      virtual ~ParticleMaterialSPHBoundary(){};

      //! create material instance of matching type with parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
    };

  }  // namespace PAR

  class ParticleMaterialSPHBoundaryType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ParticleMaterialSPHBoundaryType"; };

    static ParticleMaterialSPHBoundaryType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ParticleMaterialSPHBoundaryType instance_;
  };

  class ParticleMaterialSPHBoundary : public Material
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
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
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
      return INPAR::MAT::m_particle_sph_boundary;
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ParticleMaterialSPHBoundary(*this));
    }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    //! my material parameters
    MAT::PAR::ParticleMaterialSPHBoundary* params_;
  };

}  // namespace MAT

/*---------------------------------------------------------------------------*/
#endif
