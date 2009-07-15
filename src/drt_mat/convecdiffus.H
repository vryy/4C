/*!----------------------------------------------------------------------
\file convecdiffus.H

<pre>
Maintainer: ???
</pre>
*----------------------------------------------------------------------*/
#ifndef CONVECDIFFUS_H
#define CONVECDIFFUS_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class ConvecDiffus : public Parameter
    {
    public:

      /// standard constructor
      ConvecDiffus(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ConvecDiffus() { ; }

      /// @name material parameters
      //@{

      /// kinematic diffusivity
      const double diffusivity_;
      /// reaction coefficient
      const double reacoeff_;
      /// specific heat capacity
      const double shc_;

      //@}

    };  // class ConvecDiffus

  } // namespace PAR


  /*----------------------------------------------------------------------*/
  /// Wrapper for Newtonian fluid material
  class ConvecDiffus : public Material
  {
  public:

    /// construct empty material object
    ConvecDiffus();

    /// construct the material object given material parameters
    explicit ConvecDiffus(MAT::PAR::ConvecDiffus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_ConvecDiffus; }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_condif; }

    /// diffusivity
    double Diffusivity() const { return params_->diffusivity_; }
    /// reaction coefficient
    double ReaCoeff() const { return params_->reacoeff_; }
    /// specific heat capacity
    double Shc() const { return params_->shc_; }

  private:

    /// my material parameters
    MAT::PAR::ConvecDiffus* params_;

  };

}

#endif

#endif
