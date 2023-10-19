/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material according to Sutherland law with
       Arrhenius-type chemical kinetics (progress variable)

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_MAT_ARRHENIUS_PV_H
#define BACI_MAT_ARRHENIUS_PV_H



#include "baci_lib_parobjectfactory.H"
#include "baci_mat_material.H"
#include "baci_mat_par_parameter.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material with Arrhenius-type chemical kinetics (progress
    /// variable)
    class ArrheniusPV : public Parameter
    {
     public:
      /// standard constructor
      ArrheniusPV(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      ~ArrheniusPV() override { ; }

      /// @name material parameters
      //@{

      /// reference dynamic viscosity (kg/(m*s))
      const double refvisc_;
      /// reference temperature (K)
      const double reftemp_;
      /// Sutherland temperature (K)
      const double suthtemp_;
      /// Prandtl number
      const double pranum_;
      /// pre-exponential constant
      const double preexcon_;
      /// exponent of temperature dependence
      const double tempexp_;
      /// activation temperature
      const double actemp_;
      /// specific heat capacity of unburnt phase
      const double unbshc_;
      /// specific heat capacity of burnt phase
      const double burshc_;
      /// temperature of unburnt phase
      const double unbtemp_;
      /// temperature of burnt phase
      const double burtemp_;
      /// density of unburnt phase
      const double unbdens_;
      /// density of burnt phase
      const double burdens_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class ArrheniusPV

  }  // namespace PAR

  class ArrheniusPVType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ArrheniusPVType"; }

    static ArrheniusPVType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ArrheniusPVType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material with Arrhenius-type chemical kinetics (progress
  /// variable)
  class ArrheniusPV : public Material
  {
   public:
    /// construct empty material object
    ArrheniusPV();

    /// construct the material object given material parameters
    explicit ArrheniusPV(MAT::PAR::ArrheniusPV* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */

    int UniqueParObjectId() const override
    {
      return ArrheniusPVType::Instance().UniqueParObjectId();
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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_arrhenius_pv; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new ArrheniusPV(*this)); }

    /// compute temperature
    double ComputeTemperature(const double provar) const;

    /// compute density
    double ComputeDensity(const double provar) const;

    /// compute factor for scalar time derivative and convective scalar term
    double ComputeFactor(const double provar) const;

    /// compute specific heat capacity at constant pressure
    double ComputeShc(const double provar) const;

    /// compute viscosity
    double ComputeViscosity(const double temp) const;

    /// compute diffusivity
    double ComputeDiffusivity(const double temp) const;

    /// compute reaction coefficient
    double ComputeReactionCoeff(const double temp) const;

    /// return material parameters for element calculation
    //@{

    /// reference dynamic viscosity (kg/(m*s))
    double RefVisc() const { return params_->refvisc_; }
    /// reference temperature (K)
    double RefTemp() const { return params_->reftemp_; }
    /// Sutherland temperature (K)
    double SuthTemp() const { return params_->suthtemp_; }
    /// Prandtl number
    double PraNum() const { return params_->pranum_; }
    /// pre-exponential constant
    double PreExCon() const { return params_->preexcon_; }
    /// exponent of temperature dependence
    double TempExp() const { return params_->tempexp_; }
    /// activation temperature
    double AcTemp() const { return params_->actemp_; }
    /// specific heat capacity of unburnt phase
    double UnbShc() const { return params_->unbshc_; }
    /// specific heat capacity of burnt phase
    double BurShc() const { return params_->burshc_; }
    /// temperature of unburnt phase
    double UnbTemp() const { return params_->unbtemp_; }
    /// temperature of burnt phase
    double BurTemp() const { return params_->burtemp_; }
    /// density of unburnt phase
    double UnbDens() const { return params_->unbdens_; }
    /// density of burnt phase
    double BurDens() const { return params_->burdens_; }

    //@}

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::ArrheniusPV* params_;
  };

}  // namespace MAT

#endif
