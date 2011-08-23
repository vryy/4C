/*----------------------------------------------------------------------*/
/*!
\file thermoplasticlinelast.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following ideal von Mises plasticity and a linear thermo-elastic material
       law (Thermo St.Venant Kirchhoff).

       Combination of the two material model
       - PlastiLinElastic
       - ThermoStVenantKirchhoff

       plastic materail with kinematic hardening
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress \sigma_y = const.

       geometric linear, for small strains including temperature dependency

       example input line:
       MAT 1 MAT_Struct_ThrPlasticLinElast YOUNG 206.9 NUE 0.29 DENS 0.0
             THEXPANS 1.72e-5 INITTEMP 293.15 YIELD 0.45 KINHARD 0.0 TOL 1.0e-6

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef THERMOPLASTICLINELAST_H_
#define THERMOPLASTICLINELAST_H_

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

using namespace Teuchos;

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for plastic lin elastic Material
    class ThermoPlasticLinElast : public Parameter
    {
    public:

      /// standard constructor
      ThermoPlasticLinElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ThermoPlasticLinElast() { ; }

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
      /// yield stress (constant for perfect plasticity)
      const double yield_;
      /// kinematic hardening
      const double kinhard_;
      // tolerance for local Newton iteration
      const double abstol_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ThermoThermoPlasticLinElast

  }  // namespace PAR


class ThermoPlasticLinElastType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ThermoPlasticLinElastType"; }

  static ThermoPlasticLinElastType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static ThermoPlasticLinElastType instance_;

};  // ThermoThermoPlasticLinElastType


  /*----------------------------------------------------------------------*/
  /// wrapper for linear thermo-elasto-plastic material
  class ThermoPlasticLinElast : public Material
  {
  public:

    /// construct empty material object
    ThermoPlasticLinElast();

    /// construct the material object given material parameters
    explicit ThermoPlasticLinElast(MAT::PAR::ThermoPlasticLinElast* params);

    //! @name Packing and Unpacking

    //!  \brief Return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return ThermoPlasticLinElastType::Instance().UniqueParObjectId(); }

    //!  \brief Pack this class so it can be communicated
    //!
    //!  Resizes the vector data and stores all information of a class in it.
    //!  The first information to be stored in data has to be the
    //!  unique parobject id delivered by UniqueParObjectId() which will then
    //!  identify the exact class on the receiving processor.
    //!
    //!  \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    //!  \brief Unpack data from a char vector into this class
    //!
    //!  The vector data contains all information to rebuild the
    //!  exact copy of an instance of a class on a different processor.
    //!  The first entry in data has to be an integer which is the unique
    //!  parobject id defined at the top of this file and delivered by
    //!  UniqueParObjectId().
    //!
    //!  \param data (in) : vector storing all data to be unpacked into this
    //!  instance.
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_thermopllinelast; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ThermoPlasticLinElast( *this ) ); }

    /// Initialize internal stress variables
    void Setup(
      const int numgp ///< number of Gauss points
      );

    /// Update internal stress variables
    void Update();

    /// Reset internal stress variables
    void Reset();

    /// Evaluate material
    void Evaluate(
      const LINALG::Matrix<6,1>& linstrain,
      LINALG::Matrix<6,1>& plstrain, ///< return the plastic strains to the element
      const int gp,  ///< current Gaus point
      Teuchos::ParameterList& params,  ///< parameter list for communication
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>& stress  ///< 2nd PK-stress
      );

    // computes stress
    void Stress(
      const double p,  //!< volumetric stress tensor
      const LINALG::Matrix<6,1>& devstress,  //!< deviatoric stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
      );

    void RelStress(
      const LINALG::Matrix<6,1>& devstrain,  //!< deviatoric stress tensor
      const LINALG::Matrix<6,1>& beta,  //!< back stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& eta //!< relative stress
      );

    /// computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6,6>& cmat);

    /// computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic(
      LINALG::Matrix<6,6>& cmat,
      double Dgamma,
      double G,
      double q,
      const LINALG::Matrix<6,1> flowvector,
      const LINALG::Matrix<6,1> eta,
      double heaviside,
      double Hkin
      );

    /// computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic2(
      LINALG::Matrix<6,6>& cmat,
      double Dgamma,
      double q,
      LINALG::Matrix<6,1> unitflow
      );

    void SetupContinuumCmatElastoPlastic(
      LINALG::Matrix<6,6>& cmat,
      double Dgamma,
      double q,
      LINALG::Matrix<6,1> unitflow
      );

    /// Return density
    virtual double Density() const { return params_->density_; }

    /// Check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ && (strainplcurr_!=Teuchos::null) );
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Finite difference check for debugging purposes
    void FDCheck(
      LINALG::Matrix<6,1>& stress,  //!< updated stress sigma_n+1
      LINALG::Matrix<6,6>& cmatFD, //!< material tangent calculated with FD of stresses
      LINALG::Matrix<6,1>& beta,  //!< updated back stresses
      double p,  //!< volumetric stress
      const LINALG::Matrix<6,1>& strain,  //!< elastic strain vector
      double Dgamma,  //!< plastic multiplier
      double G,  //!< shear modulus
      double qbar,  //!< elastic trial von Mises effective stress
      double kappa, //!< bulk modulus
      LINALG::Matrix<6,1>& N, // flow vector
      double heaviside  //!< Heaviside function
      );

    /// Check implementation of material tangent for debugging purposes
    void CheckImpl();

    //! Flag plastic step was called
    bool plastic_step;

    /// @name temperature specific methods
    //@{

    /// main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(
      const LINALG::Matrix<1,1>& Ntemp,  // temperature of element
      LINALG::Matrix<6,1>& ctemp,
      LINALG::Matrix<6,1>& stresstemp
      );

    /// computes the constant temperature fraction of the stress
    void Stempconst(
      LINALG::Matrix<6,1>& ctemp,
      LINALG::Matrix<6,1>& stempconst
      );

    /// computes temperature dependent isotropic elasticity tensor in matrix
    /// notion for 3d
    void SetupCthermo(LINALG::Matrix<6,1>& ctemp);

    /// calculates stress-temperature modulus
    double STModulus();

    /// Initial temperature
    double InitTemp() const { return params_->thetainit_; }

    //@}


  private:

    /// my material parameters
    MAT::PAR::ThermoPlasticLinElast* params_;

    /// plastic history vector
    // old plastic strain at t_n
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    // current plastic strain at t_n+1
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    // old accumulated plastic strain at t_n
    Teuchos::RCP< vector< LINALG::Matrix<1,1> > > strainbarpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    // current accumulated plastic strain at t_n+1
    Teuchos::RCP< vector< LINALG::Matrix<1,1> > > strainbarplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    // old back stress at t_n
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresslast_;  //!< \f${\beta}_{n}\f$
    // current back stress at t_n+1
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresscurr_;  //!< \f${\beta}_{n+1}\f$


    /// indicator if #Initialize routine has been called
    bool isinit_;

  };  // class ThermoThermoPlasticLinElast : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
#endif  // THERMOPLASTICLINELAST_H_


