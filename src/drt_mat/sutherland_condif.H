/*----------------------------------------------------------------------*/
/*!
\file sutherland_condif.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SUTHERLAND_CONDIF_H
#define SUTHERLAND_CONDIF_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for temperature-dependent viscosity according to Sutherland law
    ///
    /// This object exists only once for each read Newton fluid.
    class SutherlandCondif : public Parameter
    {
    public:

      /// standard constructor
      SutherlandCondif(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~SutherlandCondif() { ; }

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

      //@}

    };  // class SutherlandCondif

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// temperature-dependent viscosity according to Sutherland law
  class SutherlandCondif : public Material
  {
  public:

    /// construct empty material object
    SutherlandCondif();

    /// construct the material object given material parameters
    explicit SutherlandCondif(MAT::PAR::SutherlandCondif* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_SutherlandCondif; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_sutherland_condif; }

    /// compute diffusivity
    double ComputeDiffusivity(const double temp) const;

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

    //@}

  private:

    /// my material parameters
    MAT::PAR::SutherlandCondif* params_;
  };

}

#endif

#endif
