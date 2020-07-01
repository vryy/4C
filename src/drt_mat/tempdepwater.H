/*----------------------------------------------------------------------*/
/*! \file
\brief temperature-dependent water according to "VDI Waermeatlas"

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef TEMPDEPWATER_H
#define TEMPDEPWATER_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for temperature-dependent water
    class TempDepWater : public Parameter
    {
     public:
      /// standard constructor
      TempDepWater(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~TempDepWater() { ; }

      /// @name material parameters
      //@{

      /// critical density (kg/m^3)
      const double critdens_;
      /// critical temperature (K)
      const double crittemp_;
      /// specific heat capacity
      const double shc_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class TempDepWater

  }  // namespace PAR

  class TempDepWaterType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "TempDepWaterType"; }

    static TempDepWaterType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static TempDepWaterType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// temperature-dependent water
  class TempDepWater : public Material
  {
   public:
    /// construct empty material object
    TempDepWater();

    /// construct the material object given material parameters
    explicit TempDepWater(MAT::PAR::TempDepWater* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const
    {
      return TempDepWaterType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_tempdepwater; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new TempDepWater(*this)); }

    /// compute viscosity
    double ComputeViscosity(const double temp) const;

    /// compute diffusivity
    double ComputeDiffusivity(const double temp) const;

    /// compute density
    double ComputeDensity(const double temp) const;

    /// return material parameters for element calculation
    //@{

    /// critical density (kg/m^3)
    double CritDens() const { return params_->critdens_; }
    /// reference temperature (K)
    double CritTemp() const { return params_->crittemp_; }
    /// specific heat capacity
    double Shc() const { return params_->shc_; }

    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::TempDepWater* params_;
  };

}  // namespace MAT

#endif
