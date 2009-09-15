/*!----------------------------------------------------------------------
\file mixfrac_fluid.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*----------------------------------------------------------------------*/
#ifndef MIXFRAC_FLUID_H
#define MIXFRAC_FLUID_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for fluid material coupled to mixture-fraction approach
    class MixFracFluid : public Parameter
    {
    public:

      /// standard constructor
      MixFracFluid(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ MixFracFluid() { ; }

      /// @name material parameters
      //@{

      /// dynamic viscosity
      const double viscosity_;

      /// equation-of-state factor a
      const double eosfaca_;

      /// equation-of-state factor b
      const double eosfacb_;

      //@}

    };  // class MixFracFluid

  } // namespace PAR


  /*----------------------------------------------------------------------*/
  /// Wrapper for fluid material coupled to mixture-fraction approach
  class  MixFracFluid : public Material
  {
  public:

    /// construct empty material object
    MixFracFluid();

    /// construct the material object given material parameters
    explicit MixFracFluid(MAT::PAR::MixFracFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_MixFracFluid; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_mixfrac_fluid; }

    /// compute density
    double ComputeDensity(const double mixfrac) const;

    /// dynamic viscosity
    double Viscosity() const { return params_->viscosity_; }
    /// equation-of-state factor a
    double EosFacA() const { return params_->eosfaca_; }
    /// equation-of-state factor b
    double EosFacB() const { return params_->eosfacb_; }

  private:

    /// my material parameters
    MAT::PAR::MixFracFluid* params_;

  };

}

#endif

#endif
