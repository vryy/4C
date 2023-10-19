/*----------------------------------------------------------------------*/
/*! \file
\brief material according to mixture-fraction approach

\level 2

*----------------------------------------------------------------------*/
#ifndef BACI_MAT_MIXFRAC_H
#define BACI_MAT_MIXFRAC_H



#include "baci_lib_parobjectfactory.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for material according to mixture-fraction approach
    class MixFrac : public Parameter
    {
     public:
      /// standard constructor
      MixFrac(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      ~MixFrac() override { ; }

      /// @name material parameters
      //@{

      /// kinematic viscosity
      const double kinvisc_;

      /// kinematic diffusivity
      const double kindiff_;

      /// equation-of-state factor a
      const double eosfaca_;

      /// equation-of-state factor b
      const double eosfacb_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class MixFrac

  }  // namespace PAR

  class MixFracType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "MixFracType"; }

    static MixFracType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MixFracType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for material according to mixture-fraction approach
  class MixFrac : public Material
  {
   public:
    /// construct empty material object
    MixFrac();

    /// construct the material object given material parameters
    explicit MixFrac(MAT::PAR::MixFrac* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return MixFracType::Instance().UniqueParObjectId(); }

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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_mixfrac; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new MixFrac(*this)); }

    /// compute dynamic viscosity
    double ComputeViscosity(const double temp) const;

    /// compute dynamic diffusivity
    double ComputeDiffusivity(const double temp) const;

    /// compute density
    double ComputeDensity(const double mixfrac) const;

    /// kinematic viscosity
    double KinVisc() const { return params_->kinvisc_; }
    /// kinematic diffusivity
    double KinDiff() const { return params_->kindiff_; }
    /// equation-of-state factor a
    double EosFacA() const { return params_->eosfaca_; }
    /// equation-of-state factor b
    double EosFacB() const { return params_->eosfacb_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::MixFrac* params_;
  };

}  // namespace MAT

#endif
