/*!----------------------------------------------------------------------
\file meshFreeAMmat.H

\brief material specifications of the meshfree method

\level 2

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/
#ifndef SRC_DRT_MAT_MESHFREEAMMAT_H_
#define SRC_DRT_MAT_MESHFREEAMMAT_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
namespace PAR
{

  /*----------------------------------------------------------------------*/
  /// material parameters for particles
  class MeshFreeAMMat : public Parameter
  {
  public:

    /// standard constructor
    MeshFreeAMMat(
      Teuchos::RCP<MAT::PAR::Material> matdata
      );

    /// destructor
    virtual ~MeshFreeAMMat() { ; }

    /// @name material parameters
    //@{

    /// initial mass
    const double initMass_;
    /// initial radius
    const double initRadius_;
    /// initial temperature
    const double initTemperature_;
    /// specific heat - constant pressure - solid state
    const double CPS_;
    /// specific heat - constant pressure - liquid state
    const double CPL_;
    /// specific latent heat - solid <-> liquid
    const double SL_latent_heat_max_;
    /// transition temperature - solid <-> liquid
    const double SL_transitionTemperature_;
    /// volumetric thermal expansion coefficient - solid state - usually 3 times the linear coefficient
    const double S_thermalExpansion_;
    /// volumetric thermal expansion coefficient - liquid state
    const double L_thermalExpansion_;
    /// volumetric thermal expansion coefficient - solid <-> liquid - Beware! It is linked to latent heat instead of temperature!
    const double SL_thermalExpansion_;

    //@}

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  };  // class MeshFreeAMMat

  } // namespace PAR

class MeshFreeAMMatType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "MeshFreeAMMatType"; }

  static MeshFreeAMMatType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static MeshFreeAMMatType instance_;
};

/*----------------------------------------------------------------------*/
/// mesh free material
class MeshFreeAMMat : public Material
{
public:

  /// construct empty material object
  MeshFreeAMMat();

  /// construct the material object given material parameters
  explicit MeshFreeAMMat(MAT::PAR::MeshFreeAMMat* params);

  //! @name Packing and Unpacking

  /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
  */
  virtual int UniqueParObjectId() const { return MeshFreeAMMatType::Instance().UniqueParObjectId(); }

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
  virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_meshFreeAMmat; }

  /// return copy of this material object
  virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new MeshFreeAMMat( *this ) ); }

  /// return initial mass
  double InitMass() const { return params_->initMass_; }

  /// return initial radius
  double InitRadius() const { return params_->initRadius_; }

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const { return params_; }

private:

  /// my material parameters
  MAT::PAR::MeshFreeAMMat* params_;
};

}

#endif /* SRC_DRT_MAT_MESHFREEAMMAT_H_ */
