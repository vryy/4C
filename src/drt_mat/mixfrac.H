/*----------------------------------------------------------------------*/
/*! \file
\brief material according to mixture-fraction approach

\level 2

*----------------------------------------------------------------------*/
#ifndef MIXFRAC_H
#define MIXFRAC_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

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
      virtual ~MixFrac() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class MixFrac

  }  // namespace PAR

  class MixFracType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MixFracType"; }

    static MixFracType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return MixFracType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_mixfrac; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new MixFrac(*this)); }

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
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::MixFrac* params_;
  };

}  // namespace MAT

#endif
