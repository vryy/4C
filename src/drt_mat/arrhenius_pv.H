/*----------------------------------------------------------------------*/
/*!
\file arrhenius_pv.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ARRHENIUS_PV_H
#define ARRHENIUS_PV_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion material according to Sutherland law with Arrhenius-type chemical kinetics (progress variable)
    ///
    /// This object exists only once
    class ArrheniusPV : public Parameter
    {
    public:

      /// standard constructor
      ArrheniusPV(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ArrheniusPV() { ; }

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
      /// pre-exponential constant
      const double preexcon_;
      /// exponent of temperature dependence
      const double tempexp_;
      /// activation temperature
      const double actemp_;
      /// temperature of unburnt phase
      const double unbtemp_;
      /// temperature of burnt phase
      const double burtemp_;

      //@}

    };  // class ArrheniusPV

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// convection-diffusion material according to Sutherland law with Arrhenius-type chemical kinetics (progress variable)
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
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_ArrheniusPV; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_arrhenius_pv; }

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
    /// pre-exponential constant
    double PreExCon() const { return params_->preexcon_; }
    /// exponent of temperature dependence
    double TempExp() const { return params_->tempexp_; }
    /// activation temperature
    double AcTemp() const { return params_->actemp_; }
    /// temperature of unburnt phase
    double UnbTemp() const { return params_->unbtemp_; }
    /// temperature of burnt phase
    double BurTemp() const { return params_->burtemp_; }

    //@}

  private:

    /// my material parameters
    MAT::PAR::ArrheniusPV* params_;
  };

}

#endif

#endif
