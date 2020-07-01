/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material with simplified chemical kinetics

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FERECH_PV_H
#define FERECH_PV_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material with simplified chemical
    /// kinetics due to Ferziger and Echekki (1993) (original version and
    /// modification by Poinsot and Veynante (2005)) (progress variable)
    class FerEchPV : public Parameter
    {
     public:
      /// standard constructor
      FerEchPV(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FerEchPV() { ; }

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
      /// reaction-rate constant
      const double reacratecon_;
      /// critical value of progress variable
      const double pvcrit_;
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
      /// modification factor (0.0=original, 1.0=modified)
      const double mod_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class FerEchPV

  }  // namespace PAR


  class FerEchPVType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "FerEchPVType"; }

    static FerEchPVType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static FerEchPVType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material with simplified chemical
  /// kinetics due to Ferziger and Echekki (1993) (original version and
  /// modification by Poinsot and Veynante (2005)) (progress variable)
  class FerEchPV : public Material
  {
   public:
    /// construct empty material object
    FerEchPV();

    /// construct the material object given material parameters
    explicit FerEchPV(MAT::PAR::FerEchPV* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return FerEchPVType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_ferech_pv; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new FerEchPV(*this)); }

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
    /// reaction-rate constant
    double ReacRateCon() const { return params_->reacratecon_; }
    /// critical value of progress variable
    double PvCrit() const { return params_->pvcrit_; }
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
    /// modification factor (0.0=original, 1.0=modified)
    double Mod() const { return params_->mod_; }

    //@}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::FerEchPV* params_;
  };

}  // namespace MAT

#endif
