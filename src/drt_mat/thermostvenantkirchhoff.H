/*----------------------------------------------------------------------*/
/*!
\file thermostvenantkirchhoff.H

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 02/10 |
 *----------------------------------------------------------------------*/
#ifndef THERMOSTVENANTKIRCHHOFF_H_
#define THERMOSTVENANTKIRCHHOFF_H_

#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

namespace PAR
{

  /*----------------------------------------------------------------------*/
  /// material parameters for de St. Venant--Kirchhoff with temperature
  /// dependent term
  ///
  /// <h3>Input line</h3>
  /// MAT 1 MAT_Struct_ThrStVenantK YOUNG 400 NUE 0.3 DENS 1 THEXPANS 1 INITTEMP 20
  class ThermoStVenantKirchhoff : public Parameter
  {
  public:

    /// standard constructor
    ThermoStVenantKirchhoff(Teuchos::RCP<MAT::PAR::Material> matdata);

    /// destructor
    virtual ~ThermoStVenantKirchhoff() { ; }

    /// @name material parameters
    //@{

    /// Young's modulus
    const double youngs_;
    /// Possion's ratio
    const double poissonratio_;
    /// mass density
    const double density_;
    /// linear coefficient of thermal expansion
    const double thermexpans_;
    /// initial temperature (constant) \f \Theta_0 \f
    const double thetainit_;

    //@}

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  };  // class ThermoStVenantKirchhoff

}  // namespace PAR

class ThermoStVenantKirchhoffType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ThermoStVenantKirchhoffType"; }

  static ThermoStVenantKirchhoffType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static ThermoStVenantKirchhoffType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for St.-Venant-Kirchhoff material with temperature term
  class ThermoStVenantKirchhoff : public Material
  {
  public:

    /// construct empty material object
    ThermoStVenantKirchhoff();

    /// construct the material object given material parameters
    explicit ThermoStVenantKirchhoff(MAT::PAR::ThermoStVenantKirchhoff* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return ThermoStVenantKirchhoffType::Instance().UniqueParObjectId();
    }

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
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_thermostvenant;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ThermoStVenantKirchhoff( *this ) ); }

    // evaluates stresses for 3d
    void Evaluate(
      const Epetra_SerialDenseVector* glstrain_e,
      Epetra_SerialDenseMatrix* cmat_e,
      Epetra_SerialDenseVector* stress_e
      );

    // evaluates stresses for 3d
    void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,
      LINALG::Matrix<6,6>& cmat,
      LINALG::Matrix<6,1>& stress
      );

    /// Main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(
      const LINALG::Matrix<1,1>& Ntemp,  // temperature of element
      LINALG::Matrix<6,1>& ctemp,
      LINALG::Matrix<6,1>& stresstemp
      );

    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6,6>& cmat);

    // computes isotropic elasticity tensor in matrix notion for 2d
    // plane strain, rotational symmetry
    void SetupCmat2d(Epetra_SerialDenseMatrix* cmat);

    // computes temperature dependent isotropic elasticity tensor in matrix
    // notion for 3d
    void SetupCthermo(LINALG::Matrix<6,1>& ctemp);

    // computes temperature dependent isotropic elasticity tensor in matrix
    // notion for 2d plane strain, rotational symmetry
    void SetupCthermo2d(Epetra_SerialDenseMatrix* ctemp);

    // computes the constant temperature fraction of the stress
    void Stempconst(
      LINALG::Matrix<6,1>& ctemp,
      LINALG::Matrix<6,1>& stempconst
      );

    // calculates stress-temperature modulus
    double STModulus();

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio \f \nu \f
    double PoissonRatio() const { return params_->poissonratio_; }

    /// density
    double Density() const { return params_->density_; }

    /// Initial temperature
    double InitTemp() const { return params_->thetainit_; }

    /// shear modulus \f \mu \f
    double ShearMod() const
    {
      return 0.5*params_->youngs_/(1.0+params_->poissonratio_);
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ThermoStVenantKirchhoff* params_;

  };  // ThermoStVenantKirchhoff
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
#endif  // THERMOSTVENANTKIRCHHOFF_H
