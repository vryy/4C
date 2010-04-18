/*!----------------------------------------------------------------------
\file scatra.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_H
#define SCATRA_MAT_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// parameters for scalar transport material
    class ScatraMat : public Parameter
    {
    public:

      /// standard constructor
      ScatraMat(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ScatraMat() { ; }

      /// @name material parameters
      //@{

      /// kinematic diffusivity
      const double diffusivity_;
      /// reaction coefficient
      const double reacoeff_;

      //@}

    };  // class Scatra

  } // namespace PAR


  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class ScatraMat : public Material
  {
  public:

    /// construct empty material object
    ScatraMat();

    /// construct the material object given material parameters
    explicit ScatraMat(MAT::PAR::ScatraMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_ScatraMat; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_scatra; }

    /// diffusivity
    double Diffusivity() const { return params_->diffusivity_; }
    /// reaction coefficient
    double ReaCoeff() const { return params_->reacoeff_; }

  private:

    /// my material parameters
    MAT::PAR::ScatraMat* params_;

  };

}

#endif

#endif
