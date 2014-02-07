/*----------------------------------------------------------------------*/
/*!
\file plasticneohooke.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following rate-independent von Mises plasticity with combined
       (nonlinear) isotropic and kinematic hardening and a nearly
       incompressible Neo Hookean material law.

       Refer also to the master thesis of A.A. van der Stelt (2009).

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PLASTICNEOHOOKE_H_
#define PLASTICNEOHOOKE_H_


#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for plastic neo-Hooke
    class PlasticNeoHooke : public Parameter
    {
    public:

      /// standard constructor
      PlasticNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      /// isotropic hardening
      const double isohard_;
      /// yield stress
      const double yield_;
      /// inf yield stress for nonlinear iso hardening
      const double infyield_;
      /// exponent for nonlinear iso hardening
      const double exp_;
      /// kinematic hardening
      const double kinhard_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PlasticNeoHooke

  }  // namespace PAR

class PlasticNeoHookeType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "PlasticNeoHookeType"; }

  static PlasticNeoHookeType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static PlasticNeoHookeType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for Plastic Neo-Hookean Material
  class PlasticNeoHooke : public So3Material
  {
  public:

    /// construct empty material object
    PlasticNeoHooke();

    /// construct the material object given material parameters
    explicit PlasticNeoHooke(MAT::PAR::PlasticNeoHooke* params);

    //! @name Packing and Unpacking

    //!  \brief Return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return PlasticNeoHookeType::Instance().UniqueParObjectId(); }

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plneohooke; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new PlasticNeoHooke( *this ) ); }

    /// Initialize internal stress variables
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Update internal stress variables
    virtual void Update();

    /// Evaluate material
    virtual void Evaluate
    (
      const LINALG::Matrix<3,3>* defgrd,                  ///< deformation gradient
      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,    ///< Green-Lagrange strains
      Teuchos::ParameterList& params,                     ///< parameter list for communication
      LINALG::Matrix<NUM_STRESS_3D,1>* stress,            ///< 2nd PK-stress
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,   ///< material stiffness matrix
      const int eleGID
    );

    /// Return density
    virtual double Density() const { return params_->density_; }

    /// Check if history variables are already initialized
    bool Initialized() const { return (isinit_ && (histplasticrcgcurr_!=Teuchos::null)); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::PlasticNeoHooke* params_;

    /// plastic history deformation
    // current plastic strain (right cauchy-green tensor)
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histplasticrcgcurr_;
    // last converged state of plastic right cauchy-green tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histplasticrcglast_;
    // current equivalent plastic strain
    Teuchos::RCP<std::vector<LINALG::Matrix<1,1> > > histeplasticscurr_;
    // last converged state of equivalent plastic strain
    Teuchos::RCP<std::vector<LINALG::Matrix<1,1> > > histeplasticslast_;
    // current backward stress /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > histplasticbstresscurr_;
    // backward stress of last converged state  /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > histplasticbstresslast_;
    // current eigenvalues /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > lamdacurr_;
    // eigenvalues of the last converged state  /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > lamdalast_;
    // matrix of the current three eigenvectors in 3D /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > nsetcurr_;
    // matrix of the three eigenvectors in 3D of last converged state /*added*/
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > nsetlast_;

    /// indicator if #Initialize routine has been called
    bool isinit_;

  };  // class PlasticNeoHooke : public Material
}  // namespace MAT

#endif  // PLASTICNEOHOOKE_H_
