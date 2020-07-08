/*----------------------------------------------------------------------*/
/*! \file
\brief temperature-dependent gas according to Sutherland law

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef SUTHERLAND_H
#define SUTHERLAND_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for material with temperature dependence according to Sutherland law
    class Sutherland : public Parameter
    {
     public:
      /// standard constructor
      Sutherland(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Sutherland() { ; }

      /// @name material parameters
      //@{

      /// reference dynamic viscosity (kg/(m*s))
      const double refvisc_;
      /// reference temperature (K)
      const double reftemp_;
      /// Sutherland temperature (K)
      const double suthtemp_;
      /// specific heat capacity
      const double shc_;
      /// Prandtl number
      const double pranum_;
      /// (initial) thermodynamic pressure
      const double thermpress_;
      /// specific gas constant R
      const double gasconst_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Sutherland

  }  // namespace PAR

  class SutherlandType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "SutherlandType"; }

    static SutherlandType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static SutherlandType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// material with temperature dependence according to Sutherland law
  class Sutherland : public Material
  {
   public:
    /// construct empty material object
    Sutherland();

    /// construct the material object given material parameters
    explicit Sutherland(MAT::PAR::Sutherland* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return SutherlandType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_sutherland; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Sutherland(*this)); }

    /// compute viscosity
    double ComputeViscosity(const double temp) const;

    /// compute diffusivity
    double ComputeDiffusivity(const double temp) const;

    /// compute density
    double ComputeDensity(const double temp, const double thermpress) const;

    /// return material parameters for element calculation
    //@{

    /// reference dynamic viscosity (kg/(m*s))
    double RefVisc() const { return params_->refvisc_; }
    /// reference temperature (K)
    double RefTemp() const { return params_->reftemp_; }
    /// Sutherland temperature (K)
    double SuthTemp() const { return params_->suthtemp_; }
    /// specific heat capacity
    double Shc() const { return params_->shc_; }
    /// Prandtl number
    double PraNum() const { return params_->pranum_; }
    /// (initial) thermodynamic pressure
    double ThermPress() const { return params_->thermpress_; }
    /// specific gas constant R
    double GasConst() const { return params_->gasconst_; }

    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::Sutherland* params_;
  };

}  // namespace MAT

#endif
