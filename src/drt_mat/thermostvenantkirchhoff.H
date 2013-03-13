/*----------------------------------------------------------------------*/
/*!
\file thermostvenantkirchhoff.H
\brief St.Venant Kirchhoff with an additional temperature dependent term
       describing heat expansion

       example input line:
       MAT 1   MAT_Struct_ThrStVenantK   YOUNG 1.48e8 NUE 0.3 DENS 9.130e-6
       THEXPANS 1.72e-5 INITTEMP 293.15

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


/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{

namespace PAR
{

  /*----------------------------------------------------------------------*/
  //! material parameters for de St. Venant--Kirchhoff with temperature
  //! dependent term
  //!
  //! <h3>Input line</h3>
  //! MAT 1 MAT_Struct_ThrStVenantK YOUNG 400 NUE 0.3 DENS 1 THEXPANS 1 INITTEMP 20
  class ThermoStVenantKirchhoff : public Parameter
  {
  public:

    //! standard constructor
    ThermoStVenantKirchhoff(Teuchos::RCP<MAT::PAR::Material> matdata);

    //! destructor
    virtual ~ThermoStVenantKirchhoff() { ; }

    //! @name material parameters
    //@{

    //! Young's modulus (temperature dependent --> polynomial expression)
    const std::vector<double> youngs_;
    //! Possion's ratio \f$ \nu \f$
    const double poissonratio_;
    //! mass density \f$ \rho \f$
    const double density_;
    //! linear coefficient of thermal expansion  \f$ \alpha_T \f$
    const double thermexpans_;
    //! heat capacity \f$ C_V \f$
    const double capa_;
    //! heat conductivity \f$ k \f$
    const double conduct_;
    //! initial temperature (constant) \f$ \theta_0 \f$
    const double thetainit_;

    //@}

    //! create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  };  // class ThermoStVenantKirchhoff

}  // namespace PAR

class ThermoStVenantKirchhoffType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ThermoStVenantKirchhoffType"; }

  static ThermoStVenantKirchhoffType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ThermoStVenantKirchhoffType instance_;
};

  /*----------------------------------------------------------------------*/
  //! Wrapper for St.-Venant-Kirchhoff material with temperature term
  class ThermoStVenantKirchhoff : public Material
  {
  public:

    //! construct empty material object
    ThermoStVenantKirchhoff();

    //! construct the material object given material parameters
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
    */
    virtual void Pack(
      DRT::PackBuffer& data  //!< (i/o): char vector to store class information
      ) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
    */
    virtual void Unpack(
      const std::vector<char>& data  //!< (i) : vector storing all data to be unpacked into this instance.
      );

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_thermostvenant;
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp( new ThermoStVenantKirchhoff( *this ) );
    }

    //! evaluates stresses for 3d
    void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,  //!< Green-Lagrange strain
      LINALG::Matrix<6,6>& cmat,  //!< elastic material tangent
      LINALG::Matrix<6,1>& stress,  //!< stress
      Teuchos::ParameterList& params  //!< parameter list
      );

    //! main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(
      const LINALG::Matrix<1,1>& Ntemp,  //!< temperature of element (N . T)
      LINALG::Matrix<6,1>& ctemp,  //!< temperature dependent material tangent
      LINALG::Matrix<6,1>& stresstemp,  //!< temperature dependent stress term
      Teuchos::ParameterList& params  //!< parameter list
      );

    //! computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(
      LINALG::Matrix<6,6>& cmat,
      Teuchos::ParameterList params
      );

    //! computes temperature dependent isotropic elasticity tensor in matrix
    //! notion for 3d
    void SetupCthermo(
      LINALG::Matrix<6,1>& ctemp,
      Teuchos::ParameterList& params  //!< parameter list
      );

    //! calculates stress-temperature modulus
    double STModulus(
      Teuchos::ParameterList& params  //!< parameter list
      ) const;

    //! Young's modulus, here we suppose a constant value for E
    double Youngs() const
    {
      dserror("return value 'E' not available due to possibly temperature-dependent Young's modulus.");
      return (params_->youngs_)[0];
    }

    //! Poisson's ratio \f$ \nu \f$
    double PoissonRatio() const { return params_->poissonratio_; }

    //! density \f$ \rho \f$
    virtual double Density() const { return params_->density_; }

    //! conductivity \f$ k \f$
    double Conductivity() const { return params_->conduct_; }

    //! material capacity \f$ C_V \f$
    double Capacity() const { return params_->capa_; }

    //! Initial temperature \f$ \theta_0 \f$
    double InitTemp() const { return params_->thetainit_; }

    //! shear modulus \f$ \mu \f$, here we suppose a constant value for E
    double ShearMod() const
    {
      dserror("return value 'mu' not available due to possibly temperature-dependent Young's modulus");
      return 0.5*(params_->youngs_)[0]/(1.0+params_->poissonratio_);
    }

    //! calculate temperature dependent material parameter and return value
    double GetMatParameterAtTempnp(
      const std::vector<double>* paramvector,  //!< (i) given parameter is a vector
      const double& tempnp  // tmpr (i) current temperature
      ) const ;

    //! Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    //! my material parameters
    MAT::PAR::ThermoStVenantKirchhoff* params_;

  };  // ThermoStVenantKirchhoff
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // THERMOSTVENANTKIRCHHOFF_H
