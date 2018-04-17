/*----------------------------------------------------------------------*/
/*!
\file newtonianfluid.H

\brief Newtonian fluid material

\level 1

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef NEWTONIANFLUID_H
#define NEWTONIANFLUID_H



#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for Newtonian fluid
    ///
    /// This object exists only once for each read Newton fluid.
    class NewtonianFluid : public Parameter
    {
    public:

      /// standard constructor
      NewtonianFluid(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~NewtonianFluid() { ; }

      /// @name material parameters
      //@{

      /// kinematic or dynamic viscosity
      const double viscosity_;
      /// density
      const double density_;
      /// surface tension coefficient
      const double gamma_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class NewtonianFluid

  }  // namespace PAR

class NewtonianFluidType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "NewtonianFluidType"; }

  static NewtonianFluidType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static NewtonianFluidType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for Newtonian fluid material
  ///
  /// This object exists (several times) at every element
  class NewtonianFluid : public Material
  {
  public:

    /// construct empty material object
    NewtonianFluid();

    /// construct the material object given material parameters
    explicit NewtonianFluid(MAT::PAR::NewtonianFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return NewtonianFluidType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_fluid; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new NewtonianFluid( *this ) ); }

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return density
    virtual double Density() const { return params_->density_; }

    /// return surface tension coefficient
    double Gamma() const { return params_->gamma_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::NewtonianFluid* params_;
  };

}

#endif

