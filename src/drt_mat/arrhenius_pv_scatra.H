/*----------------------------------------------------------------------*/
/*!
\file arrhenius_pv_scatra.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ARRHENIUS_PV_SCATRA_H
#define ARRHENIUS_PV_SCATRA_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material with Arrhenius-type chemical kinetics (progress variable)
    class ArrheniusPVScatra : public Parameter
    {
    public:

      /// standard constructor
      ArrheniusPVScatra(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ArrheniusPVScatra() { ; }

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
      /// flame speed
      const double flamespeed_;
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

    };  // class ArrheniusPV

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material with Arrhenius-type chemical kinetics (progress variable)
  class ArrheniusPVScatra : public Material
  {
  public:

    /// construct empty material object
    ArrheniusPVScatra();

    /// construct the material object given material parameters
    explicit ArrheniusPVScatra(MAT::PAR::ArrheniusPVScatra* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_ArrheniusPVScatra; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_arrhenius_pv_scatra; }

    /// compute rate parameter alpha
    double ComputeAlpha() const;

    /// compute rate parameter beta
    double ComputeBeta() const;

    /// compute diffusive flame thickness
    double ComputeDiffFlameThickness() const;

    /// compute flame thickness due to Blint (1986)
    double ComputeBlintFlameThickness() const;

    /// compute specific heat capacity at constant pressure
    double ComputeShc(const double provar) const;

    /// compute temperature
    double ComputeTemperature(const double provar) const;

    /// compute density
    double ComputeDensity(const double provar) const;

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
    /// flame speed
    double FlameSpeed() const { return params_->flamespeed_; }
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

  private:

    /// my material parameters
    MAT::PAR::ArrheniusPVScatra* params_;
  };

}

#endif

#endif
