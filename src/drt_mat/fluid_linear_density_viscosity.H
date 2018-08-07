/*----------------------------------------------------------------------*/
/*!
\file fluid_linear_density_viscosity.H

\brief Linear law (pressure-dependent) for the density and the viscosity

<pre>
\level 1

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef LINEARDENSITYVISCOSITY_H
#define LINEARDENSITYVISCOSITY_H



#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for fluid with linear law (pressure-dependent) for
    /// the density and the viscosity
    ///
    /// This object exists only once for each read fluid.
    class LinearDensityViscosity : public Parameter
    {
    public:

      /// standard constructor
      LinearDensityViscosity(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~LinearDensityViscosity() { ; }

      /// @name material parameters
      //@{

      /// reference density
      const double refdensity_;
      /// reference viscosity
      const double refviscosity_;
      /// reference pressure
      const double refpressure_;
      /// density-pressure coefficient
      const double coeffdensity_;
      /// viscosity-pressure coefficient
      const double coeffviscosity_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class LinearDensityViscosity

  }  // namespace PAR

class LinearDensityViscosityType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "LinearDensityViscosityType"; }

  static LinearDensityViscosityType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static LinearDensityViscosityType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for fluid with linear law (pressure-dependent) for
  /// the density and the viscosity
  ///
  /// This object exists (several times) at every element
  class LinearDensityViscosity : public Material
  {
  public:

    /// construct empty material object
    LinearDensityViscosity();

    /// construct the material object given material parameters
    explicit LinearDensityViscosity(MAT::PAR::LinearDensityViscosity* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return LinearDensityViscosityType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_fluid_linear_density_viscosity; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new LinearDensityViscosity( *this ) ); }

    /// compute density
    double ComputeDensity(const double press) const;

    /// compute viscosity
    double ComputeViscosity(const double press) const;

    /// return material parameters for element calculation
    //@{

    /// return reference density
    double RefDensity() const { return params_->refdensity_; }

    /// return reference viscosity
    double RefViscosity() const { return params_->refviscosity_; }

    /// return reference pressure
    double RefPressure() const { return params_->refpressure_; }

    /// return density-pressure coefficient
    double CoeffDensity() const { return params_->coeffdensity_; }

    /// return viscosity-pressure coefficient
    double CoeffViscosity() const { return params_->coeffviscosity_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::LinearDensityViscosity* params_;
  };

}

#endif

