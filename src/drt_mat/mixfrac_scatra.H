/*!----------------------------------------------------------------------
\file mixfrac_scatra.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*----------------------------------------------------------------------*/
#ifndef MIXFRAC_SCATRA_H
#define MIXFRAC_SCATRA_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material according to mixture-fraction approach
    class MixFracScatra : public Parameter
    {
    public:

      /// standard constructor
      MixFracScatra(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ MixFracScatra() { ; }

      /// @name material parameters
      //@{

      /// dynamic diffusivity
      const double diffusivity_;

      /// equation-of-state factor a
      const double eosfaca_;

      /// equation-of-state factor b
      const double eosfacb_;

      //@}

    };  // class MixFracScatra

  } // namespace PAR


  /*----------------------------------------------------------------------*/
  /// Wrapper for scalar transport material according to mixture-fraction approach
  class  MixFracScatra : public Material
  {
  public:

    /// construct empty material object
    MixFracScatra();

    /// construct the material object given material parameters
    explicit MixFracScatra(MAT::PAR::MixFracScatra* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_MixFracScatra; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_mixfrac_scatra; }

    /// compute density
    double ComputeDensity(const double mixfrac) const;

    /// dynamic diffusivity
    double Diffusivity() const { return params_->diffusivity_; }
    /// equation-of-state factor a
    double EosFacA() const { return params_->eosfaca_; }
    /// equation-of-state factor b
    double EosFacB() const { return params_->eosfacb_; }

  private:

    /// my material parameters
    MAT::PAR::MixFracScatra* params_;

  };

}

#endif

#endif
