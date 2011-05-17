/*----------------------------------------------------------------------*/
/*!
\file plasticlinelast.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following ideal von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       perfect plasticity:
        - no hardening allowed
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress \sigma_y = const.

       geometric linear, for small strains

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PLASTICLINELAST_H_
#define PLASTICLINELAST_H_

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
    class PlasticLinElast : public Parameter
    {
    public:

      /// standard constructor
      PlasticLinElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticLinElast() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      /// yield stress (constant for perfect plasticity)
      const double yield_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PlasticLinElast

  }  // namespace PAR


class PlasticLinElastType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "PlasticLinElastType"; }

  static PlasticLinElastType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static PlasticLinElastType instance_;
};  // PlasticLinElastType


  /*----------------------------------------------------------------------*/
  /// Wrapper for Plastic linear elastic (St.Venant Kirchhoff) Material
  class PlasticLinElast : public Material
  {
  public:

    /// construct empty material object
    PlasticLinElast();

    /// construct the material object given material parameters
    explicit PlasticLinElast(MAT::PAR::PlasticLinElast* params);

    //! @name Packing and Unpacking

    //!  \brief Return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return PlasticLinElastType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_pllinelast; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new PlasticLinElast( *this ) ); }

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

    /// computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6,6>& cmat);

    /// computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic(
      LINALG::Matrix<6,6>& cmat,
      double Dgamma,
      double G,
      double q,
      const LINALG::Matrix<6,1> flowvector,
      double heaviside
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
    double Density() const { return params_->density_; }

    /// Check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ && (strainplcurr_!=Teuchos::null) );
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Finite difference check for debugging purposes
    void FDCheck(
      const LINALG::Matrix<6,1>& strain,  //!< elastic strain vector
      double Dgamma,  //!< plastic multiplier
      double G,  //!< shear modulus
      double kappa, //!< bulk modulus
      const LINALG::Matrix<6,1>& N, // flow vector
      double heaviside,  //!< Heaviside function
      double p,  //!< volumetric stress
      LINALG::Matrix<6,1>& stress  //!< updated stress sigma_n+1
      );

    /// Check implementation of material tangent for debugging purposes
    void CheckImpl();


  private:

    /// my material parameters
    MAT::PAR::PlasticLinElast* params_;

    /// plastic history deformation
    // old plastic strain at t_n
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    // current plastic strain at t_n+1
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$

    /// indicator if #Initialize routine has been called
    bool isinit_;

  };  // class PlasticLinElast : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
#endif  // PLASTICLINELAST_H_


