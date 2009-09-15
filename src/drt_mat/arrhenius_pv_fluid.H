/*----------------------------------------------------------------------*/
/*!
\file arrhenius_pv_fluid.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ARRHENIUS_PV_FLUID_H
#define ARRHENIUS_PV_FLUID_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for fluid material coupled to Arrhenius-type chemical kinetics (progress variable)
    class ArrheniusPVFluid : public Parameter
    {
    public:

      /// standard constructor
      ArrheniusPVFluid(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ArrheniusPVFluid() { ; }

      /// @name material parameters
      //@{

      /// reference dynamic viscosity (kg/(m*s))
      const double refvisc_;
      /// reference temperature (K)
      const double reftemp_;
      /// Sutherland temperature (K)
      const double suthtemp_;
      /// temperature of unburnt phase
      const double unbtemp_;
      /// temperature of burnt phase
      const double burtemp_;
      /// density of unburnt phase
      const double unbdens_;
      /// density of burnt phase
      const double burdens_;

      //@}

    };  // class ArrheniusPVFluid

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// wrapper for fluid material coupled to Arrhenius-type chemical kinetics (progress variable)
  class ArrheniusPVFluid : public Material
  {
  public:

    /// construct empty material object
    ArrheniusPVFluid();

    /// construct the material object given material parameters
    explicit ArrheniusPVFluid(MAT::PAR::ArrheniusPVFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_ArrheniusPVFluid; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_arrhenius_pv_fluid; }

    /// compute temperature
    double ComputeTemperature(const double provar) const;

    /// compute density
    double ComputeDensity(const double provar) const;

    /// compute viscosity
    double ComputeViscosity(const double temp) const;

    /// return material parameters for element calculation
    //@{

    /// reference dynamic viscosity (kg/(m*s))
    double RefVisc() const { return params_->refvisc_; }
    /// reference temperature (K)
    double RefTemp() const { return params_->reftemp_; }
    /// Sutherland temperature (K)
    double SuthTemp() const { return params_->suthtemp_; }
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
    MAT::PAR::ArrheniusPVFluid* params_;
  };

}

#endif

#endif
