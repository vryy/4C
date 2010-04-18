/*----------------------------------------------------------------------*/
/*!
\file modpowerlaw.H

<pre>
Maintainer: Ursula Mayer
			mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef MODPOWERLAW_H
#define MODPOWERLAW_H

#ifdef CCADISCRET

#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class ModPowerLaw : public Parameter
    {
    public:

      /// standard constructor
      ModPowerLaw(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ModPowerLaw() { ; }

      /// @name material parameters
      //@{
     const double m_cons_; 		    /* consistency */
     const double delta_;       	/* safety factor */
     const double a_exp_;			/* exponent */
     const double density_;     	/* density */
      //@}

    };  // class ModPowerLaw

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Nonlinear viscosity according to a modified power law
  class ModPowerLaw : public Material
  {
  public:

    /// construct empty material object
    ModPowerLaw();

    /// construct the material object given material parameters
    explicit ModPowerLaw(MAT::PAR::ModPowerLaw* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_ModPowerLaw; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_modpowerlaw; }

    /// return material parameters for element calculation

    /// consistency constant
    double MCons() const { return params_->m_cons_; }
    /// safety factor
    double Delta() const { return params_->delta_; }
    /// exponent
    double AExp() const { return params_->a_exp_; }
    /// density
    double Density() const { return params_->density_; }

  private:

    /// my material parameters
    MAT::PAR::ModPowerLaw* params_;

  };

}

#endif

#endif
